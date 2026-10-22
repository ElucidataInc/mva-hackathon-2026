#!/usr/bin/env python3
"""
Single-pass streaming profiler for a VCF/BCF-decoded VCF.

Computes (without needing any external tool):
  * header / provenance / build / annotation detection
  * record & allele counts, variant classes, bi- vs multi-allelic sites
  * Ts/Tv, indel length distribution, per-chromosome breakdowns
  * sort order, duplicate records, REF concordance against a FASTA
  * FILTER, QUAL, depth and GATK INFO-annotation distributions
  * known vs novel (rsID / DB flag / gnomAD AF)
  * VEP / SnpEff / ClinVar annotation summaries
  * per-sample genotype counts, Ts/Tv, Het/HomAlt, missingness, DP, GQ,
    singletons, private variants, phasing, allele balance, chrX/chrY metrics
  * Mendelian errors (PED), somatic metrics (TLOD/NLOD/VAF/F1R2/F2R1/SB)
  * variant density + mean depth in genomic bins (for Manhattan-style maps)
  * a down-sampled genotype dosage matrix for relatedness / PCA / ROH
"""

import argparse
import os
import random
import re
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcf_common import (Fasta, TRANSITIONS, chrom_sort_key, classify_allele, detect_annotations,
                        detect_callers, detect_mode, first_float, guess_build, in_par, is_autosome,
                        norm_chrom, open_vcf, parse_header, parse_info, read_header_text, read_ped,
                        snv_bases, to_float, write_json)

IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
GATK_INFO_FIELDS = ["QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum", "VQSLOD", "InbreedingCoeff",
                    "ExcessHet", "BaseQRankSum", "AS_VQSLOD"]
SOMATIC_INFO_FIELDS = ["TLOD", "NLOD", "NALOD", "GERMQ", "POPAF", "ROQ", "STRANDQ", "SEQQ", "CONTQ",
                       "MBQ", "MMQ", "MFRL", "MPOS", "SomaticEVS", "SSC", "SPV"]


class Reservoir:
    """Fixed-size uniform random sample of a numeric stream."""

    def __init__(self, k=200_000, seed=11):
        self.k, self.n, self.v = k, 0, []
        self.rng = random.Random(seed)

    def add(self, x):
        if x is None:
            return
        self.n += 1
        if len(self.v) < self.k:
            self.v.append(x)
        else:
            j = self.rng.randrange(self.n)
            if j < self.k:
                self.v[j] = x

    def summary(self):
        if not self.v:
            return None
        s = sorted(self.v)
        n = len(s)
        q = lambda p: s[min(n - 1, int(p * (n - 1)))]
        return {"n": self.n, "mean": sum(s) / n, "min": s[0], "p05": q(0.05), "p25": q(.25), "median": q(.5),
                "p75": q(.75), "p95": q(.95), "max": s[-1]}


class AnnotationCounter:
    """Summarises VEP CSQ, SnpEff ANN, bcftools BCSQ, ClinVar and gnomAD annotation."""

    def __init__(self, ann_meta):
        self.meta = ann_meta
        self.csq_fields = ann_meta.get("VEP_CSQ", [])
        self.consequence = Counter()
        self.impact = Counter()
        self.biotype = Counter()
        self.genes = Counter()
        self.clinsig = Counter()
        self.gnomad_af = Reservoir(100_000)
        self.n_annotated = 0
        self.gnomad_keys = [k for k in ann_meta.get("gnomAD_fields", []) if re.search(r"AF", k)]

    def add(self, info):
        hit = False
        if "CSQ" in info and self.csq_fields:
            hit = True
            f = {n: i for i, n in enumerate(self.csq_fields)}
            cons, best, bts, genes = set(), 0, set(), set()
            for entry in str(info["CSQ"]).split(","):
                p = entry.split("|")
                g = lambda k: p[f[k]] if k in f and f[k] < len(p) else ""
                for c in g("Consequence").split("&"):
                    if c:
                        cons.add(c)
                best = max(best, IMPACT_RANK.get(g("IMPACT"), 0))
                if g("BIOTYPE"):
                    bts.add(g("BIOTYPE"))
                if g("SYMBOL"):
                    genes.add(g("SYMBOL"))
                for k in ("gnomADg_AF", "gnomAD_AF", "gnomADe_AF", "MAX_AF"):
                    if g(k):
                        self.gnomad_af.add(to_float(g(k).split("&")[0]))
                        break
                if g("CLIN_SIG"):
                    for c in g("CLIN_SIG").split("&"):
                        self.clinsig[c] += 1
            self._commit(cons, best, bts, genes)
        elif "ANN" in info:
            hit = True
            cons, best, bts, genes = set(), 0, set(), set()
            for entry in str(info["ANN"]).split(","):
                p = entry.split("|")
                if len(p) < 8:
                    continue
                for c in p[1].split("&"):
                    cons.add(c)
                best = max(best, IMPACT_RANK.get(p[2], 0))
                if p[3]:
                    genes.add(p[3])
                if p[7]:
                    bts.add(p[7])
            self._commit(cons, best, bts, genes)
        elif "BCSQ" in info:
            hit = True
            cons, genes, bts = set(), set(), set()
            for entry in str(info["BCSQ"]).split(","):
                p = entry.split("|")
                cons.add(p[0].lstrip("@*"))
                if len(p) > 1:
                    genes.add(p[1])
                if len(p) > 3:
                    bts.add(p[3])
            self._commit(cons, 0, bts, genes)
        for k in ("CLNSIG", "CLINSIG"):
            if k in info:
                hit = True
                for c in re.split(r"[,|&/]", str(info[k])):
                    self.clinsig[c] += 1
        for k in self.gnomad_keys:
            if k in info:
                self.gnomad_af.add(first_float(info[k]))
                break
        if hit:
            self.n_annotated += 1

    def _commit(self, cons, best, bts, genes):
        for c in cons:
            self.consequence[c] += 1
        inv = {v: k for k, v in IMPACT_RANK.items()}
        self.impact[inv.get(best, "NA")] += 1
        for b in bts:
            self.biotype[b] += 1
        for g in genes:
            self.genes[g] += 1

    def result(self):
        af = self.gnomad_af
        bins = [0, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.5, 1.0001]
        labels = ["<1e-5", "1e-5–1e-4", "1e-4–1e-3", "0.1–1%", "1–5%", "5–50%", ">50%"]
        afc = Counter()
        for v in af.v:
            for i in range(len(bins) - 1):
                if bins[i] <= v < bins[i + 1]:
                    afc[labels[i]] += 1
                    break
        return {
            "sources": self.meta, "n_annotated_records": self.n_annotated,
            "consequence": self.consequence.most_common(40), "impact": dict(self.impact),
            "biotype": self.biotype.most_common(20), "top_genes": self.genes.most_common(30),
            "clinsig": self.clinsig.most_common(20),
            "gnomad_af_bins": [[l, afc.get(l, 0)] for l in labels], "gnomad_af_summary": af.summary(),
        }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--ref", help="reference FASTA (for REF concordance)")
    ap.add_argument("--ped", help="PED file (for Mendelian errors)")
    ap.add_argument("--bin-size", type=int, default=1_000_000)
    ap.add_argument("--max-gt-sites", type=int, default=100_000,
                    help="max sites kept in the dosage matrix used for kinship / PCA / ROH")
    ap.add_argument("--tumor", help="tumor sample name (somatic VCFs)")
    ap.add_argument("--normal", help="normal sample name (somatic VCFs)")
    ap.add_argument("--max-records", type=int, default=0, help="stop after N records (0 = all)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()

    # ---------------------------------------------------------------- header
    hdr = parse_header(a.vcf)
    htext = read_header_text(a.vcf)
    samples = hdr["samples"]
    nS = len(samples)
    build = guess_build(hdr)
    mode = detect_mode(hdr, htext)
    ann_meta = detect_annotations(hdr)
    callers = detect_callers(hdr, htext)
    contig_len = {c["ID"]: c["length"] for c in hdr["contigs"] if c["length"]}
    contig_order = {c["ID"]: i for i, c in enumerate(hdr["contigs"])}
    fasta = None
    if a.ref and os.path.exists(a.ref):
        try:
            fasta = Fasta(a.ref)
        except Exception as e:  # noqa
            print(f"[profile] cannot open FASTA: {e}", file=sys.stderr)

    tumor = a.tumor or hdr.get("tumor_sample")
    normal = a.normal or hdr.get("normal_sample")
    t_idx = samples.index(tumor) if tumor in samples else None
    n_idx = samples.index(normal) if normal in samples else None

    trios_named = read_ped(a.ped)
    trios = [(samples.index(c), samples.index(f), samples.index(m), c) for c, f, m in trios_named
             if c in samples and f in samples and m in samples]

    # ---------------------------------------------------------------- state
    C = Counter()                       # global scalar counters
    allele_class = Counter()
    site_class = Counter()
    n_alt_hist = Counter()
    indel_len = Counter()
    filters = Counter()
    filter_by_class = defaultdict(Counter)
    subst = Counter()
    qual = Reservoir()
    site_dp = Reservoir()
    gatk_info = {k: Reservoir(100_000) for k in GATK_INFO_FIELDS if k in hdr["info"]}
    som_info = {k: Reservoir(100_000) for k in SOMATIC_INFO_FIELDS if k in hdr["info"]}
    ann = AnnotationCounter(ann_meta)
    known = Counter()
    per_chrom = defaultdict(Counter)
    bins = defaultdict(lambda: [0, 0.0, 0, 0])   # (chrom,bin) -> n, dp_sum, dp_n, n_pass
    site_missing_hist = Counter()
    ref_mismatch_examples = []
    sort_issues = []
    dup_examples = []

    # per-sample arrays
    z = lambda: [0] * nS
    s_called, s_miss, s_homref, s_het, s_homalt, s_hetalt = z(), z(), z(), z(), z(), z()
    s_ts, s_tv, s_snv, s_indel, s_single, s_private = z(), z(), z(), z(), z(), z()
    s_dp_sum, s_dp_n, s_gq_sum, s_gq_n = z(), z(), z(), z()
    s_phased, s_phased_het, s_haploid = z(), z(), z()
    s_x_het, s_x_called, s_y_called, s_y_nonref = z(), z(), z(), z()
    s_ab = [[0] * 20 for _ in range(nS)]
    s_ab_low = z()
    s_ps = [set() for _ in range(nS)]
    s_vaf = [Reservoir(30_000, seed=i) for i in range(nS)] if nS <= 50 else []
    s_dpres = [Reservoir(30_000, seed=100 + i) for i in range(nS)] if nS <= 50 else []
    s_f1r2 = [Reservoir(30_000, seed=200 + i) for i in range(nS)] if mode == "somatic" else []
    s_sb = [Reservoir(30_000, seed=300 + i) for i in range(nS)] if mode == "somatic" else []
    vaf_all = Reservoir()
    chrY_sites = 0
    trio_tested = [0] * len(trios)
    trio_err = [0] * len(trios)

    gt_rng = random.Random(7)
    gt_rows, gt_pos, gt_seen = [], [], 0

    prev_chrom, prev_pos = None, -1
    closed_chroms = set()
    pos_keys = set()
    max_examples = 25

    with open_vcf(a.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                C["malformed_lines"] += 1
                continue
            chrom, pos_s, vid, ref, alt, qual_s, flt, info_s = f[:8]
            try:
                pos = int(pos_s)
            except ValueError:
                C["malformed_lines"] += 1
                continue
            C["records"] += 1
            if a.max_records and C["records"] > a.max_records:
                C["records"] -= 1
                C["truncated"] = 1
                break

            # ---------------- sort order & duplicates
            if chrom != prev_chrom:
                if chrom in closed_chroms:
                    C["noncontiguous_chrom_blocks"] += 1
                    if len(sort_issues) < max_examples:
                        sort_issues.append(f"{chrom} reappears at {chrom}:{pos} after {prev_chrom}")
                if prev_chrom is not None:
                    closed_chroms.add(prev_chrom)
                    if chrom in contig_order and prev_chrom in contig_order and \
                            contig_order[chrom] < contig_order[prev_chrom]:
                        C["chrom_out_of_header_order"] += 1
                prev_chrom, prev_pos = chrom, -1
                pos_keys = set()
            if pos < prev_pos:
                C["unsorted_positions"] += 1
                if len(sort_issues) < max_examples:
                    sort_issues.append(f"{chrom}:{pos} after {chrom}:{prev_pos}")
            if pos != prev_pos:
                pos_keys = set()
            else:
                C["records_sharing_position"] += 1
            key = (ref, alt)
            if key in pos_keys:
                C["duplicate_records_exact"] += 1
                if len(dup_examples) < max_examples:
                    dup_examples.append(f"{chrom}:{pos} {ref}>{alt}")
            pos_keys.add(key)
            prev_pos = max(prev_pos, pos)
            if chrom not in contig_len and hdr["contigs"]:
                C["chrom_not_in_header"] += 1

            alts = alt.split(",")
            real_alts = [x for x in alts if x not in ("*", "<NON_REF>", "<*>", ".")]
            n_alt_hist[len(real_alts)] += 1
            if len(real_alts) == 0:
                C["no_alt_sites"] += 1
            elif len(real_alts) == 1:
                C["biallelic_sites"] += 1
            else:
                C["multiallelic_sites"] += 1
            if "*" in alts:
                C["sites_with_spanning_deletion"] += 1

            # ---------------- REF checks
            if re.search(r"[^ACGTNacgtn]", ref):
                C["ref_non_acgtn"] += 1
            if fasta is not None:
                if fasta.has(chrom):
                    seq = fasta.fetch(chrom, pos, pos + len(ref) - 1)
                    if seq is None:
                        C["ref_out_of_bounds"] += 1
                    else:
                        C["ref_checked"] += 1
                        r = ref.upper()
                        mism = any(x != y and x != "N" and y != "N" for x, y in zip(seq, r))
                        if mism:
                            C["ref_mismatch"] += 1
                            if len(ref_mismatch_examples) < max_examples:
                                ref_mismatch_examples.append(f"{chrom}:{pos} VCF={ref[:20]} FASTA={seq[:20]}")
                        else:
                            C["ref_match"] += 1
                else:
                    C["ref_contig_missing_in_fasta"] += 1

            # ---------------- classes
            classes = []
            is_snv_site = True
            ts_site = None
            for x in alts:
                cl, d = classify_allele(ref, x)
                allele_class[cl] += 1
                classes.append(cl)
                per_chrom[chrom]["allele_" + cl] += 1
                if cl in ("insertion", "deletion", "complex_indel"):
                    indel_len[d] += 1
                if cl == "SNV":
                    b = snv_bases(ref, x)
                    subst[f"{b[0]}>{b[1]}"] += 1
                    if b in TRANSITIONS:
                        C["ts"] += 1
                        ts_site = True if ts_site is None else ts_site
                    else:
                        C["tv"] += 1
                        ts_site = False
                elif cl not in ("spanning_deletion", "ref_block"):
                    is_snv_site = False
            core = sorted({c for c in classes if c not in ("spanning_deletion", "ref_block")})
            if not core:
                sclass = "ref_block" if "ref_block" in classes else "other"
            elif len(core) == 1:
                sclass = core[0]
            elif set(core) <= {"insertion", "deletion", "complex_indel"}:
                sclass = "indel_mixed"
            else:
                sclass = "mixed"
            site_class[sclass] += 1
            per_chrom[chrom]["records"] += 1
            per_chrom[chrom]["site_" + sclass] += 1
            is_indel_site = sclass in ("insertion", "deletion", "complex_indel", "indel_mixed")
            vtype = "SNV" if sclass == "SNV" else ("INDEL" if is_indel_site else "OTHER")

            # ---------------- filters / qual
            flist = ["PASS"] if flt == "PASS" else (["."] if flt == "." else flt.split(";"))
            for x in flist:
                filters[x] += 1
                filter_by_class[vtype][x] += 1
            is_pass = flt in ("PASS", ".")
            if is_pass:
                C["pass_records"] += 1
                per_chrom[chrom]["pass"] += 1
            q = to_float(qual_s)
            qual.add(q)

            info = parse_info(info_s)
            for k, r in gatk_info.items():
                if k in info:
                    r.add(first_float(info[k]))
            for k, r in som_info.items():
                if k in info:
                    r.add(first_float(info[k]))
            ann.add(info)

            # ---------------- known vs novel
            k_rs = vid.startswith("rs")
            k_id = vid != "."
            k_db = "DB" in info
            k_gn = any(k in info and to_float(str(info[k]).split(",")[0]) is not None for k in ann.gnomad_keys)
            k_any = k_rs or k_db or k_gn or k_id
            known["with_ID"] += k_id
            known["with_rsID"] += k_rs
            known["DB_flag"] += k_db
            known["in_gnomAD"] += k_gn
            known["known_any"] += k_any
            known[f"{vtype}_total"] += 1
            known[f"{vtype}_known"] += k_any
            if is_snv_site and len(real_alts) == 1 and ts_site is not None:
                known[("known" if k_any else "novel") + ("_ts" if ts_site else "_tv")] += 1

            # ---------------- depth
            dp = to_float(info.get("DP"))
            fmt = f[8].split(":") if len(f) > 8 else []
            fidx = {k: i for i, k in enumerate(fmt)}
            sample_fields = [s.split(":") for s in f[9:]] if nS else []
            if dp is None and "DP" in fidx:
                tot, got = 0, False
                for sf in sample_fields:
                    if fidx["DP"] < len(sf):
                        v = to_float(sf[fidx["DP"]])
                        if v is not None:
                            tot += v
                            got = True
                dp = tot if got else None
            site_dp.add(dp)

            b = (chrom, pos // a.bin_size)
            bb = bins[b]
            bb[0] += 1
            if dp is not None:
                bb[1] += dp
                bb[2] += 1
            if is_pass:
                bb[3] += 1

            # ---------------- genotypes
            if not nS or "GT" not in fidx:
                continue
            gi = fidx["GT"]
            dpi, gqi, adi, afi = fidx.get("DP"), fidx.get("GQ"), fidx.get("AD"), fidx.get("AF")
            psi, f1i, f2i, sbi = fidx.get("PS"), fidx.get("F1R2"), fidx.get("F2R1"), fidx.get("SB")
            nc = norm_chrom(chrom)
            on_x = nc == "X" and not in_par(pos, build)
            on_y = nc == "Y"
            if on_y:
                chrY_sites += 1
            auto = is_autosome(chrom)
            biallelic = len(alts) == 1 and len(real_alts) == 1
            biallelic_snv = biallelic and sclass == "SNV"
            ac = [0] * (len(alts) + 1)
            carriers = [[] for _ in range(len(alts) + 1)]
            n_missing = 0
            dosages = [] if (biallelic_snv and auto and is_pass) else None
            gts = [None] * nS if trios else None

            for i, sf in enumerate(sample_fields):
                gt = sf[gi] if gi < len(sf) else "."
                phased = "|" in gt
                al = re.split(r"[/|]", gt)
                if "." in al or gt == "":
                    s_miss[i] += 1
                    n_missing += 1
                    if dosages is not None:
                        dosages.append(-1)
                    continue
                try:
                    al = [int(x) for x in al]
                except ValueError:
                    s_miss[i] += 1
                    n_missing += 1
                    if dosages is not None:
                        dosages.append(-1)
                    continue
                s_called[i] += 1
                if gts is not None:
                    gts[i] = al
                if len(al) == 1:
                    s_haploid[i] += 1
                nonref = [x for x in al if x != 0]
                for x in nonref:
                    if x < len(ac):
                        ac[x] += 1
                        if len(carriers[x]) < 2 and (not carriers[x] or carriers[x][-1] != i):
                            carriers[x].append(i)
                if dosages is not None:
                    dosages.append(min(2, len(nonref)) if len(al) == 2 else 2 * len(nonref))
                if not nonref:
                    s_homref[i] += 1
                    gclass = "homref"
                elif len(set(al)) == 1:
                    s_homalt[i] += 1
                    gclass = "homalt"
                    per_chrom[chrom]["homalt"] += 1
                else:
                    s_het[i] += 1
                    gclass = "het"
                    per_chrom[chrom]["het"] += 1
                    if 0 not in al:
                        s_hetalt[i] += 1
                if phased:
                    s_phased[i] += 1
                    if gclass == "het":
                        s_phased_het[i] += 1
                    if psi is not None and psi < len(sf) and sf[psi] not in (".", "") and len(s_ps[i]) < 1_000_000:
                        s_ps[i].add(sf[psi])
                if nonref:
                    if sclass == "SNV":
                        s_snv[i] += 1
                        if biallelic_snv and ts_site is not None:
                            if ts_site:
                                s_ts[i] += 1
                            else:
                                s_tv[i] += 1
                    elif is_indel_site:
                        s_indel[i] += 1
                if on_x:
                    s_x_called[i] += 1
                    if gclass == "het":
                        s_x_het[i] += 1
                if on_y:
                    s_y_called[i] += 1
                    if nonref:
                        s_y_nonref[i] += 1
                if dpi is not None and dpi < len(sf):
                    v = to_float(sf[dpi])
                    if v is not None:
                        s_dp_sum[i] += v
                        s_dp_n[i] += 1
                        if s_dpres:
                            s_dpres[i].add(v)
                if gqi is not None and gqi < len(sf):
                    v = to_float(sf[gqi])
                    if v is not None:
                        s_gq_sum[i] += v
                        s_gq_n[i] += 1
                # VAF / allele balance
                vaf = None
                ad = None
                if adi is not None and adi < len(sf) and sf[adi] not in (".", ""):
                    try:
                        ad = [int(x) if x != "." else 0 for x in sf[adi].split(",")]
                    except ValueError:
                        ad = None
                if nonref:
                    if afi is not None and afi < len(sf):
                        vaf = first_float(sf[afi])
                    if vaf is None and ad and sum(ad) > 0:
                        alt_i = max(nonref)
                        vaf = ad[alt_i] / sum(ad) if alt_i < len(ad) else None
                    if vaf is not None:
                        vaf_all.add(vaf)
                        if s_vaf:
                            s_vaf[i].add(vaf)
                if gclass == "het" and 0 in al and ad and sum(ad) > 0:
                    alt_i = max(al)
                    if alt_i < len(ad):
                        abv = ad[alt_i] / sum(ad)
                        s_ab[i][min(19, int(abv * 20))] += 1
                        if abv < 0.25 or abv > 0.75:
                            s_ab_low[i] += 1
                if s_f1r2 and nonref and f1i is not None and f2i is not None and f1i < len(sf) and f2i < len(sf):
                    try:
                        a1 = [int(x) for x in sf[f1i].split(",")]
                        a2 = [int(x) for x in sf[f2i].split(",")]
                        k = max(nonref)
                        if k < len(a1) and k < len(a2) and a1[k] + a2[k] > 0:
                            s_f1r2[i].add(a1[k] / (a1[k] + a2[k]))
                    except ValueError:
                        pass
                if s_sb and nonref and sbi is not None and sbi < len(sf):
                    try:
                        sb = [int(x) for x in sf[sbi].split(",")]
                        if len(sb) >= 4 and sb[2] + sb[3] > 0:
                            s_sb[i].add(sb[2] / (sb[2] + sb[3]))
                    except ValueError:
                        pass

            site_missing_hist[min(20, int(round(n_missing / nS * 20)))] += 1
            per_chrom[chrom]["missing_gt"] += n_missing
            per_chrom[chrom]["total_gt"] += nS
            for x in range(1, len(alts) + 1):
                if alts[x - 1] in ("*", "<NON_REF>", "<*>"):
                    continue
                if ac[x] == 1:
                    s_single[carriers[x][0]] += 1
                    C["singleton_alleles"] += 1
                if len(carriers[x]) == 1:
                    s_private[carriers[x][0]] += 1
                    C["private_alleles"] += 1
                if ac[x] == 0:
                    C["alleles_with_AC0"] += 1

            if trios and biallelic and auto:
                for ti, (c, fa, mo, _) in enumerate(trios):
                    gc, gf, gm = gts[c], gts[fa], gts[mo]
                    if gc is None or gf is None or gm is None or len(gc) != 2 or len(gf) != 2 or len(gm) != 2:
                        continue
                    trio_tested[ti] += 1
                    ok = (gc[0] in gf and gc[1] in gm) or (gc[1] in gf and gc[0] in gm)
                    if not ok:
                        trio_err[ti] += 1

            if dosages is not None and nS > 1:
                called = sum(1 for d in dosages if d >= 0)
                if called / nS >= 0.9:
                    gt_seen += 1
                    row = (dosages, chrom, pos)
                    if len(gt_rows) < a.max_gt_sites:
                        gt_rows.append(row)
                    else:
                        j = gt_rng.randrange(gt_seen)
                        if j < a.max_gt_sites:
                            gt_rows[j] = row

            if C["records"] % 100_000 == 0:
                el = time.time() - t0
                print(f"[profile] {C['records']:,} records, at {chrom}:{pos:,} ({el:.0f}s, "
                      f"{C['records'] / el:,.0f} rec/s)", file=sys.stderr, flush=True)

    # ---------------------------------------------------------------- outputs
    out = a.out
    fsize = os.path.getsize(a.vcf)
    idx = [p for p in (a.vcf + ".tbi", a.vcf + ".csi") if os.path.exists(p)]

    with open(os.path.join(out, "per_chrom.tsv"), "w") as fh:
        cls_cols = ["SNV", "MNV", "insertion", "deletion", "complex_indel", "indel_mixed", "mixed"]
        fh.write("chrom\tlength\trecords\tpass\tdensity_per_Mb\t" + "\t".join(cls_cols) +
                 "\tother\thet_gt\thomalt_gt\thet_homalt_ratio\tmissing_rate\n")
        for ch in sorted(per_chrom, key=chrom_sort_key):
            d = per_chrom[ch]
            L = contig_len.get(ch)
            dens = d["records"] / (L / 1e6) if L else ""
            cls = [d.get("site_" + c, 0) for c in cls_cols]
            other = d["records"] - sum(cls)
            ratio = d["het"] / d["homalt"] if d["homalt"] else ""
            miss = d["missing_gt"] / d["total_gt"] if d["total_gt"] else ""
            fh.write("\t".join(map(str, [ch, L or "", d["records"], d["pass"], dens] + cls +
                                   [other, d["het"], d["homalt"], ratio, miss])) + "\n")

    with open(os.path.join(out, "density_bins.tsv"), "w") as fh:
        fh.write("chrom\tbin_start\tbin_end\tn_variants\tn_pass\tmean_site_dp\n")
        for (ch, bi) in sorted(bins, key=lambda k: (chrom_sort_key(k[0]), k[1])):
            n, dps, dpn, npass = bins[(ch, bi)]
            fh.write(f"{ch}\t{bi * a.bin_size}\t{(bi + 1) * a.bin_size}\t{n}\t{npass}\t"
                     f"{dps / dpn if dpn else ''}\n")

    sample_rows = []
    for i, s in enumerate(samples):
        nr = s_het[i] + s_homalt[i]
        sample_rows.append({
            "sample": s, "called": s_called[i], "missing": s_miss[i],
            "missing_rate": s_miss[i] / (s_called[i] + s_miss[i]) if (s_called[i] + s_miss[i]) else None,
            "hom_ref": s_homref[i], "het": s_het[i], "hom_alt": s_homalt[i], "het_alt_alt": s_hetalt[i],
            "n_nonref": nr, "n_snv": s_snv[i], "n_indel": s_indel[i],
            "het_homalt_ratio": s_het[i] / s_homalt[i] if s_homalt[i] else None,
            "ts": s_ts[i], "tv": s_tv[i], "ts_tv": s_ts[i] / s_tv[i] if s_tv[i] else None,
            "singletons": s_single[i], "private": s_private[i],
            "mean_dp": s_dp_sum[i] / s_dp_n[i] if s_dp_n[i] else None,
            "mean_gq": s_gq_sum[i] / s_gq_n[i] if s_gq_n[i] else None,
            "phased_gt": s_phased[i], "phased_het": s_phased_het[i],
            "phased_het_frac": s_phased_het[i] / s_het[i] if s_het[i] else None,
            "phase_sets": len(s_ps[i]), "haploid_calls": s_haploid[i],
            "chrX_nonPAR_called": s_x_called[i], "chrX_het": s_x_het[i],
            "chrX_het_rate": s_x_het[i] / s_x_called[i] if s_x_called[i] else None,
            "chrY_called": s_y_called[i], "chrY_nonref": s_y_nonref[i],
            "chrY_call_rate": s_y_called[i] / chrY_sites if chrY_sites else None,
            "het_rate": s_het[i] / s_called[i] if s_called[i] else None,
            "het_ab_outlier_frac": s_ab_low[i] / sum(s_ab[i]) if sum(s_ab[i]) else None,
        })
    with open(os.path.join(out, "per_sample.tsv"), "w") as fh:
        if sample_rows:
            cols = list(sample_rows[0].keys())
            fh.write("\t".join(cols) + "\n")
            for r in sample_rows:
                fh.write("\t".join("" if r[c] is None else (f"{r[c]:.5g}" if isinstance(r[c], float) else str(r[c]))
                                   for c in cols) + "\n")

    if gt_rows:
        try:
            import numpy as np
            gt_rows.sort(key=lambda r: (chrom_sort_key(r[1]), r[2]))
            np.savez_compressed(os.path.join(out, "gt_matrix.npz"),
                                dosage=np.array([r[0] for r in gt_rows], dtype=np.int8),
                                chrom=np.array([r[1] for r in gt_rows]), pos=np.array([r[2] for r in gt_rows]),
                                samples=np.array(samples), n_eligible=gt_seen)
        except ImportError:
            pass

    hist = lambda res, lo, hi, nb: _hist(res.v, lo, hi, nb)
    prof = {
        "input": {"path": os.path.abspath(a.vcf), "size_bytes": fsize,
                  "compressed": open(a.vcf, "rb").read(2) == b"\x1f\x8b", "index": [os.path.basename(p) for p in idx]},
        "runtime_sec": round(time.time() - t0, 1),
        "header": {"fileformat": hdr["fileformat"], "n_header_lines": hdr["n_header_lines"],
                   "n_info": len(hdr["info"]), "n_format": len(hdr["format"]), "n_filter_defs": len(hdr["filter"]),
                   "n_contigs": len(hdr["contigs"]), "sources": hdr["sources"], "commands": hdr["commands"][:30],
                   "reference": hdr["reference"], "info_ids": sorted(hdr["info"]), "format_ids": sorted(hdr["format"]),
                   "filter_defs": {k: v.get("Description", "") for k, v in hdr["filter"].items()},
                   "tumor_sample": tumor, "normal_sample": normal},
        "build": build, "mode": mode, "callers": callers, "annotation_sources": ann_meta,
        "samples": samples, "n_samples": nS,
        "counts": dict(C), "allele_class": dict(allele_class), "site_class": dict(site_class),
        "n_alt_hist": {str(k): v for k, v in sorted(n_alt_hist.items())},
        "indel_length": {str(k): v for k, v in sorted(indel_len.items())},
        "substitutions": dict(subst), "filters": filters.most_common(),
        "filter_by_type": {k: v.most_common() for k, v in filter_by_class.items()},
        "qual": qual.summary(), "qual_hist": hist(qual, 0, max(100, (qual.summary() or {}).get("p95", 100) or 100), 50),
        "site_dp": site_dp.summary(),
        "site_dp_hist": hist(site_dp, 0, max(10, ((site_dp.summary() or {}).get("p95") or 10) * 1.5), 60),
        "gatk_info": {k: {"summary": r.summary(), "hist": _auto_hist(r.v)} for k, r in gatk_info.items() if r.v},
        "somatic_info": {k: {"summary": r.summary(), "hist": _auto_hist(r.v)} for k, r in som_info.items() if r.v},
        "known": dict(known), "annotation": ann.result(),
        "vaf_hist": hist(vaf_all, 0, 1, 50), "vaf_summary": vaf_all.summary(),
        "per_sample_vaf_hist": {samples[i]: _hist(r.v, 0, 1, 50) for i, r in enumerate(s_vaf)} if s_vaf else {},
        "per_sample_dp_summary": {samples[i]: r.summary() for i, r in enumerate(s_dpres)} if s_dpres else {},
        "per_sample_dp_hist": {samples[i]: _auto_hist(r.v, 60) for i, r in enumerate(s_dpres)} if s_dpres else {},
        "per_sample_ab_hist": {samples[i]: s_ab[i] for i in range(nS)} if nS <= 200 else {},
        "orientation_bias_hist": {samples[i]: _hist(r.v, 0, 1, 40) for i, r in enumerate(s_f1r2) if r.v},
        "strand_bias_hist": {samples[i]: _hist(r.v, 0, 1, 40) for i, r in enumerate(s_sb) if r.v},
        "site_missing_hist": [site_missing_hist.get(i, 0) for i in range(21)],
        "chrY_sites": chrY_sites,
        "mendelian": [{"child": t[3], "father": samples[t[1]], "mother": samples[t[2]], "tested": trio_tested[k],
                       "errors": trio_err[k], "error_rate": trio_err[k] / trio_tested[k] if trio_tested[k] else None}
                      for k, t in enumerate(trios)],
        "examples": {"ref_mismatch": ref_mismatch_examples, "sort_issues": sort_issues, "duplicates": dup_examples},
        "bin_size": a.bin_size, "contig_lengths": contig_len,
        "ref_fasta": os.path.abspath(a.ref) if fasta else None,
    }
    write_json(prof, os.path.join(out, "profile.json"))
    print(f"[profile] done: {C['records']:,} records, {nS} samples in {time.time() - t0:.1f}s", file=sys.stderr)


def _hist(vals, lo, hi, nb):
    if not vals or hi <= lo:
        return {"edges": [], "counts": []}
    w = (hi - lo) / nb
    c = [0] * nb
    for v in vals:
        k = int((v - lo) / w)
        c[max(0, min(nb - 1, k))] += 1
    return {"edges": [lo + i * w for i in range(nb + 1)], "counts": c}


def _auto_hist(vals, nb=50):
    if not vals:
        return {"edges": [], "counts": []}
    s = sorted(vals)
    lo, hi = s[int(0.005 * (len(s) - 1))], s[int(0.995 * (len(s) - 1))]
    if hi <= lo:
        hi = lo + 1
    return _hist(vals, lo, hi, nb)


if __name__ == "__main__":
    main()
