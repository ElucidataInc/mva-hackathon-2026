#!/usr/bin/env python3
"""
Rank rare-disease candidates from vcf_hc_vep.sh output, with a recessive (MVA) gene focus.

Input : a per-gene annotation table (one row per variant per gene; VEP + ClinVar columns), built by
        run_ranking.sh from vcf_hc_vep.sh output; optional ClinVar failed-QC list; optional Exomiser
        gene TSV (its combined score sorts genes within each tier)
Output: ranked_variants.tsv  every qualifying rare variant with its evidence class and score
        ranked_genes.tsv     one row per gene: genotype (biallelic / monoallelic), phase, tier, score
        targets.tsv          the follow-up list (tiers 1-4 + top discovery genes)
        targets.bed          variant positions of the targets (+/-50 bp) for IGV / Sanger / RNA checks

Evidence classes follow the corpus: ClinVar (allele-exact, >=1 review star), PVS1-type loss of
function, ClinGen-calibrated REVEL thresholds (Pejaver 2022), SpliceAI >= 0.2 (Walker 2023),
AlphaMissense class, CADD for non-coding. Scores order candidates; they are NOT ACMG classifications.
Only standard library is used (any python3 works).
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from urllib.parse import unquote

MVA_ESTABLISHED = {"BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "SLF2", "SMC5"}
MVA_EMERGING = {"BUB1", "MAD2L1BP", "TRIM37"}

LOF = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
       "start_lost", "transcript_ablation", "stop_lost"}
MISSENSE = {"missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant"}
SPLICE_NEAR = {"splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant",
               "splice_polypyrimidine_tract_variant"}
CLASS_SCORE = {"CLINVAR_PLP": 8, "LOF": 6, "MISSENSE_STRONG": 5, "SPLICE_PREDICTED": 4,
               "MISSENSE_MODERATE": 4, "MISSENSE_SUPPORTING": 3, "MISSENSE_VUS": 2, "MISSENSE_PRED_BENIGN": 1,
               "SPLICE_REGION": 1, "NONCODING_CADD": 1,
               # structural variants (from sv_annotate.py)
               "SV_LOF": 6, "SV_LOF_POSSIBLE": 4, "CNV_DEPTH": 3, "SV_SPLICE": 3, "SV_REGULATORY": 2,
               "SV_INTRONIC": 1, "SV_DOSAGE": 1,
               # AlphaGenome predictions (--alphagenome): not clinically calibrated, so supporting level only;
               # a splice change seen by both score_variant and predict_variant scores 4
               "AG_SPLICE": 3, "AG_EXPRESSION": 3, "AG_REGULATORY": 2}
STRONG = {"CLINVAR_PLP", "LOF", "MISSENSE_STRONG", "MISSENSE_MODERATE", "SPLICE_PREDICTED",
          "SV_LOF", "SV_LOF_POSSIBLE", "AG_SPLICE"}


def fnum(x):
    try:
        vals = [float(v) for v in re.split(r"[&,]", str(x)) if v not in ("", ".", "NA")]
        return max(vals) if vals else None
    except ValueError:
        return None


def get(row, *names):
    for n in names:
        v = row.get(n)
        if v not in (None, "", "."):
            return v
    return ""


def max_af(row):
    vals = [fnum(row.get(k)) for k in ("MAX_AF", "gnomADg_AF", "gnomADe_AF", "AF")]
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None


def classify(row, args):
    """Return (class, details) for one variant-gene row, or (None, reason)."""
    cons = set(get(row, "Consequence").split("&"))
    cs = get(row, "ClinVar_CLNSIG").lower()
    rs = get(row, "ClinVar_CLNREVSTAT").lower()
    clinvar_plp = ("pathogenic" in cs and "conflicting" not in cs and "benign" not in cs
                   and not re.search(r"no_assertion|no_interpretation|no_classification", rs))
    revel = fnum(row.get("REVEL"))
    sai_vals = [fnum(row.get(k)) for k in ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
                                           "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")]
    sai_vals = [v for v in sai_vals if v is not None]
    sai = max(sai_vals) if sai_vals else None
    cadd = fnum(row.get("CADD_PHRED"))
    am = get(row, "am_class")
    notes = []
    if revel is not None:
        notes.append(f"REVEL={revel:g}")
    if am:
        notes.append(f"AM={am}")
    if sai is not None:
        notes.append(f"SpliceAI={sai:g}")
    if cadd is not None:
        notes.append(f"CADD={cadd:g}")
    if clinvar_plp:
        return "CLINVAR_PLP", notes + [f"ClinVar={get(row, 'ClinVar_CLNSIG')} ({rs})"]
    if cons & LOF:
        return "LOF", notes
    if cons & MISSENSE:
        if revel is not None and revel >= 0.932:
            return "MISSENSE_STRONG", notes
        if (revel is not None and revel >= 0.773) or am == "likely_pathogenic":
            return "MISSENSE_MODERATE", notes
        if revel is not None and revel >= 0.644:
            return "MISSENSE_SUPPORTING", notes
        if revel is not None and revel <= 0.290 and am != "likely_pathogenic":
            return None, "missense predicted benign (REVEL <= 0.290)"
        if revel is None:
            # no REVEL: fall back on AlphaMissense + SIFT/PolyPhen (VEP --everything); weaker, so at most supporting
            sift, pp = get(row, "SIFT"), get(row, "PolyPhen")
            if am == "likely_benign":
                return "MISSENSE_PRED_BENIGN", notes + ["AlphaMissense likely benign (no REVEL)"]
            if sift.startswith("deleterious") and "low_confidence" not in sift and pp.startswith("probably_damaging"):
                return "MISSENSE_SUPPORTING", notes + [f"SIFT {sift} + PolyPhen {pp} (no REVEL)"]
        return "MISSENSE_VUS", notes
    if sai is not None and sai >= 0.2:
        return "SPLICE_PREDICTED", notes
    if cons & SPLICE_NEAR:
        return "SPLICE_REGION", notes
    if cadd is not None and cadd >= args.cadd_min:
        return "NONCODING_CADD", notes
    return None, "no qualifying consequence"


def alt_copies(gt):
    al = re.split(r"[/|]", gt)
    if "." in al:
        return 0
    return sum(1 for a in al if a not in ("0",))


def hap(gt):
    """For phased het GT return 0/1 = which haplotype carries ALT, else None."""
    if "|" not in gt:
        return None
    a = gt.split("|")
    if len(a) != 2 or a[0] == a[1]:
        return None
    return 0 if a[0] != "0" else 1


def truthy(x):
    return str(x or "").strip().lower() in ("true", "1", "yes", "y", "t")


def chrom_key(c):
    return c[3:] if c.startswith("chr") else c


def ag_classify(r):
    """AlphaGenome table row -> (class or None, score, note). Flags/thresholds are the AlphaGenome step's own."""
    ev = (r.get("ag_evidence") or "").strip().lower() or "none"
    log2fc = fnum(r.get("sv_expr_own_log2fc"))
    pv_splice = fnum(r.get("pv_splice_max_absdiff"))
    direction = (r.get("pv_sig_direction") or "").strip()
    parts = [f"AlphaGenome {ev}"]
    if truthy(r.get("ag_splicing")):
        q = max([v for v in (fnum(r.get("sv_splice_site_q")), fnum(r.get("sv_splice_usage_q")),
                             fnum(r.get("sv_splice_junction_q"))) if v is not None] or [0])
        both = pv_splice is not None and pv_splice >= 0.2
        parts.append(f"splicing change (q={q:g}" + (f", local splice track change {pv_splice:g}" if pv_splice is not None else "") + ")")
        return "AG_SPLICE", 4 if both else 3, "; ".join(parts + (["confirmed by predict_variant tracks"] if both else []))
    if truthy(r.get("ag_expression")):
        where = r.get("sv_expr_own_where") or ""
        parts.append(f"{r.get('marker_gene')} expression log2FC {log2fc:g} q={r.get('sv_expr_own_q')}" + (f" in {where}" if where else ""))
        if log2fc is not None and log2fc < 0:
            return "AG_EXPRESSION", 3, "; ".join(parts + ["reduced expression (fits loss of function)"])
        return None, 0, "; ".join(parts + ["increased expression: does not fit loss of function, not counted"])
    n_bs = fnum(r.get("pv_n_biosamples_sig")) or 0
    if truthy(r.get("ag_regulatory")) or ev in ("moderate", "strong") or n_bs >= 2:
        parts.append(f"regulatory change {r.get('sv_top_reg_modality') or r.get('pv_sig_modalities') or ''}".strip() +
                     (f" ({direction}, {int(n_bs)} cell types)" if direction or n_bs else ""))
        if direction == "up":                          # gain of activity: not a loss-of-function hit
            return None, 0, "; ".join(parts + ["increase only: not counted for a recessive loss-of-function model"])
        return "AG_REGULATORY", 2, "; ".join(parts)
    return None, 0, "; ".join(parts + ["no AlphaGenome signal"]) if ev != "none" else "AlphaGenome: no signal"


def load_alphagenome(path):
    """(chrom without chr, pos, ref, alt) -> list of AlphaGenome rows (one per marker gene)"""
    ag = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                k = (chrom_key(r["chrom"]), str(int(float(r["pos"]))), r["ref"], r["alt"])
            except (KeyError, ValueError):
                continue
            r["_cls"], r["_score"], r["_note"] = ag_classify(r)
            ag[k].append(r)
    return ag


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="05_tables/high_confidence.annotated.tsv")
    ap.add_argument("--failed-qc", help="05_tables/clinvar_plp_failed_qc.tsv (optional)")
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--max-af", type=float, default=0.01, help="rare = max population AF below this [0.01]")
    ap.add_argument("--cadd-min", type=float, default=20, help="CADD PHRED for non-coding candidates [20]")
    ap.add_argument("--genes", help="extra gene list (one symbol per line) to treat as a target panel")
    ap.add_argument("--top-discovery", type=int, default=20, help="non-panel genes to keep in targets [20]")
    ap.add_argument("--top-mva", type=int, default=10, help="variants (small + SV) in MVA genes listed in top_mva_variants.tsv [10]")
    ap.add_argument("--exomiser-genes", action="append", default=[],
                    help="Exomiser gene TSV (GENE_SYMBOL + EXOMISER_GENE_COMBINED_SCORE); repeatable, max score per gene")
    ap.add_argument("--sv-table", help="sv_candidates.tsv from sv_annotate.py: SVs join the gene genotypes")
    ap.add_argument("--alphagenome", help="AlphaGenome variant table (chrom pos ref alt marker_gene ag_* sv_* pv_* columns)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    panel = set()
    if a.genes:
        panel = {l.strip().split()[0] for l in open(a.genes) if l.strip() and not l.startswith("#")}

    def gene_set(g):
        if g in MVA_ESTABLISHED:
            return "MVA_established"
        if g in MVA_EMERGING:
            return "MVA_emerging"
        if g in panel:
            return "panel"
        return "other"

    exo = {}
    for ef in a.exomiser_genes:
        with open(ef) as fh:
            hdr = fh.readline().lstrip("#").rstrip("\n").split("\t")
            try:
                gi, si = hdr.index("GENE_SYMBOL"), hdr.index("EXOMISER_GENE_COMBINED_SCORE")
            except ValueError:
                sys.exit(f"{ef}: Exomiser TSV needs GENE_SYMBOL and EXOMISER_GENE_COMBINED_SCORE columns")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(gi, si):
                    try:
                        exo[f[gi]] = max(exo.get(f[gi], 0.0), float(f[si] or 0))
                    except ValueError:
                        pass
    if exo:
        print(f"[rank] Exomiser scores loaded for {len(exo):,} genes from {len(a.exomiser_genes)} file(s) (max per gene)")

    csv.field_size_limit(sys.maxsize)
    fh_table = open(a.table)
    rows = csv.DictReader(fh_table, delimiter="\t")
    first_cols = rows.fieldnames or []
    n_rows = 0
    seen = set()
    variants = []          # qualifying variant-gene records
    dropped = defaultdict(int)
    ag = load_alphagenome(a.alphagenome) if a.alphagenome else {}
    ag_seen = {}           # AlphaGenome site -> table genotype/phase (for sites the table assigns to another gene)
    if ag:
        print(f"[rank] AlphaGenome: {sum(len(v) for v in ag.values())} rows at {len(ag)} sites; with a class: "
              f"{sum(1 for v in ag.values() for r in v if r['_cls'])}")
    for r in rows:
        n_rows += 1
        g = get(r, "SYMBOL") or get(r, "Gene")
        k4 = (chrom_key(r["CHROM"]), r["POS"], r["REF"], r["ALT"])
        agrows = ag.get(k4, [])
        if agrows:
            ag_seen[k4] = r
        if not g and agrows:
            g = agrows[0]["marker_gene"]               # intergenic: use the AlphaGenome marker-gene region
        if not g:
            dropped["intergenic (no gene)"] += 1
            continue
        key = (r["CHROM"], r["POS"], r["REF"], r["ALT"], g)
        if key in seen:
            continue
        seen.add(key)
        af = max_af(r)
        if af is not None and af >= a.max_af and "pathogenic" not in get(r, "ClinVar_CLNSIG").lower():
            dropped["common (AF >= %g)" % a.max_af] += 1
            continue
        cls, info = classify(r, a)
        score = CLASS_SCORE.get(cls, 0)
        agr = next((x for x in agrows if x["marker_gene"] == g), None)
        if agr:
            agr["_used"] = True
            if agr["_cls"] and agr["_score"] > score:       # AlphaGenome only raises the class
                cls, score = agr["_cls"], agr["_score"]
                info = (info if isinstance(info, list) else []) + [agr["_note"]]
            elif cls:
                info = info + [agr["_note"]]
        if cls is None:
            dropped[info] += 1
            continue
        flt = get(r, "FILTER") or "PASS"
        failed_qc = flt not in ("PASS", ".")
        if failed_qc and cls != "CLINVAR_PLP":       # only known-pathogenic calls are rescued from the filtered set
            dropped["failed QC (not ClinVar P/LP)"] += 1
            continue
        if failed_qc:
            info = info + [f"FAILED QC ({flt}) - confirm"]
        src = get(r, "SOURCE")
        if src:
            info = info + ([f"in both VCFs"] if src == "both" else [f"only in VCF {src}"])
        if get(r, "GT_NOTE") not in ("", "-"):
            info = info + [get(r, "GT_NOTE")]
        gt = get(r, "GT")
        n = alt_copies(gt)
        if n == 0:
            dropped["no ALT in genotype"] += 1
            continue
        variants.append({"gene": g, "gene_set": gene_set(g), "chrom": r["CHROM"], "pos": r["POS"],
                         "ref": r["REF"], "alt": r["ALT"], "id": get(r, "ID"), "gt": gt, "copies": n,
                         "ps": get(r, "PS"), "class": cls, "score": score,
                         "consequence": get(r, "Consequence"), "impact": get(r, "IMPACT"),
                         "hgvsc": unquote(get(r, "HGVSc")), "hgvsp": unquote(get(r, "HGVSp")),
                         "max_af": "" if af is None else f"{af:.3g}",
                         "clinvar": get(r, "ClinVar_CLNSIG"), "dp": get(r, "DP"), "vaf": get(r, "VAF"),
                         "filter": flt, "ploidy": len(re.split(r"[/|]", gt)), "source": src or "",
                         "gt_note": get(r, "GT_NOTE") if get(r, "GT_NOTE") not in ("", "-") else "",
                         "evidence": "; ".join(info)})

    # AlphaGenome sites not used above: not in the table (e.g. filtered upstream) or assigned to another gene by VEP
    n_ag = 0
    for k4, lst in ag.items():
        for x in lst:
            if x.get("_used"):
                continue
            g = x.get("marker_gene") or ""
            t = ag_seen.get(k4)                        # genotype/phase from the VCF table when available
            row = {"Consequence": x.get("consequence", ""), "IMPACT": x.get("impact", ""),
                   "ClinVar_CLNSIG": x.get("clinvar_exact_clnsig", ""), "am_class": x.get("am_class", ""),
                   "MAX_AF": x.get("MAX_AF", "")}
            cls, info = classify(row, a)
            score = CLASS_SCORE.get(cls, 0)
            info = info if isinstance(info, list) else []
            if x["_cls"] and x["_score"] > score:
                cls, score = x["_cls"], x["_score"]
            if cls is None:
                dropped["AlphaGenome site: no qualifying evidence"] += 1
                continue
            info.append(x["_note"])
            af = fnum(x.get("MAX_AF"))
            if af is not None and af >= a.max_af and cls != "CLINVAR_PLP":
                dropped["common (AF >= %g)" % a.max_af] += 1
                continue
            ok = truthy(x.get("quality_ok")) if x.get("quality_ok") not in (None, "") else (x.get("filter") in ("PASS", ".", ""))
            if not ok and cls != "CLINVAR_PLP":
                dropped["AlphaGenome site failed QC (not PASS / not in high-confidence VCF)"] += 1
                continue
            if not ok:
                info.append("FAILED QC - confirm")
            if t:
                gt, ps = get(t, "GT"), get(t, "PS")
            else:
                z = (x.get("zygosity") or "").lower()
                gt = "1/1" if z.startswith("hom") else "1" if z.startswith("hemi") else "0/1"
                ps = ""
                info.append("not in the VCF table (genotype from AlphaGenome table, unphased)")
            if x.get("region_part") and x["region_part"] != "gene_body":
                info.append(f"{x['region_part']} of {g} (marker region +/-50 kb)")
            variants.append({"gene": g, "gene_set": gene_set(g), "chrom": t["CHROM"] if t else x["chrom"],
                             "pos": str(int(float(x["pos"]))), "ref": x["ref"], "alt": x["alt"], "id": x.get("fid", ""),
                             "gt": gt, "copies": alt_copies(gt) or 1, "ps": ps, "class": cls, "score": score,
                             "consequence": x.get("consequence", ""), "impact": x.get("impact", ""),
                             "hgvsc": unquote(x.get("hgvsc", "") or ""), "hgvsp": unquote(x.get("hgvsp", "") or ""),
                             "max_af": "" if af is None else f"{af:.3g}", "clinvar": x.get("clinvar_exact_clnsig", ""),
                             "dp": get(t, "DP") if t else "", "vaf": get(t, "VAF") if t else "",
                             "filter": x.get("filter") or "PASS", "ploidy": len(re.split(r"[/|]", gt)),
                             "source": get(t, "SOURCE") if t else "", "gt_note": "", "evidence": "; ".join(info)})
            n_ag += 1
    if ag:
        print(f"[rank] AlphaGenome: {n_ag} site(s) added that were not in the table for their marker gene")

    n_sv = 0
    if a.sv_table and os.path.exists(a.sv_table):
        for r in csv.DictReader(open(a.sv_table), delimiter="\t"):
            g = r["gene"]
            gs = gene_set(g)
            if r["filter"] != "PASS" and gs == "other":
                dropped["SV failed caller filters (non-panel gene)"] += 1
                continue
            gt = r["gt"] if r["gt"] not in (".", "./.", "") else "0/1"
            n = alt_copies(gt) or 1
            ev = [r["effect"]] + ([r["notes"]] if r["notes"] not in ("", "-") else [])
            if r.get("gene_info", "-") not in ("", "-"):
                ev.append("gene: " + r["gene_info"])
            try:                                    # sv_annotate.py may lower the score (low breakpoint confidence, benign overlap)
                sv_score = int(float(r.get("score") or CLASS_SCORE.get(r["class"], 1)))
            except ValueError:
                sv_score = CLASS_SCORE.get(r["class"], 1)
            if r["filter"] != "PASS":
                ev.append("SV FAILED caller filters - confirm")
            variants.append({"gene": g, "gene_set": gs, "chrom": r["chrom"], "pos": r["start"], "ref": "-",
                             "alt": f"<{r['svtype']}>", "id": r["ids"], "gt": gt, "copies": n, "ps": "",
                             "class": r["class"], "score": sv_score,
                             "gt_source": r.get("gt_source", "measured"), "gene_flag": r.get("gene_flag", "-"),
                             "depth_ratio": r.get("depth_ratio", "-"),
                             "confidence": r.get("confidence", "-"),
                             "consequence": r["effect"], "impact": "SV", "hgvsc": "", "hgvsp": "",
                             "max_af": r["pop_af"], "clinvar": "", "dp": "", "vaf": "",
                             "evidence": "; ".join(ev) + f" [{r['callers']}; {r['svtype']} {r['svlen']} bp]",
                             "filter": r["filter"], "ploidy": len(re.split(r"[/|]", gt)),
                             "is_sv": True, "svtype": r["svtype"], "sv_end": int(r["end"])})
            n_sv += 1
        print(f"[rank] {n_sv:,} SV-gene rows added from {a.sv_table}")

    failed = defaultdict(list)
    if a.failed_qc and os.path.exists(a.failed_qc):
        for r in csv.DictReader(open(a.failed_qc), delimiter="\t"):
            failed[(r.get("CHROM"), r.get("POS"))].append(r.get("FILTER", ""))

    # ---------------------------------------------------------------- gene genotypes
    by_gene = defaultdict(list)
    for v in variants:
        by_gene[v["gene"]].append(v)
    genes = []
    for g, vs in by_gene.items():
        vs.sort(key=lambda v: -v["score"])
        dels = [v for v in vs if v.get("is_sv") and v["svtype"] in ("DEL", "CNV") and v["copies"] == 1]
        for v in vs:
            if v.get("is_sv") or v["copies"] < 2:
                continue
            for d in dels:
                if v["chrom"] == d["chrom"] and int(d["pos"]) <= int(v["pos"]) <= d["sv_end"]:
                    v["copies"] = 1          # apparent homozygosity: only one copy left over the deletion
                    v["over_del"] = d
                    v["evidence"] = (v["evidence"] + "; " if v["evidence"] else "") + \
                        f"called homozygous inside a heterozygous deletion -> hemizygous, in trans with it"
                    break
        hom = [v for v in vs if v["copies"] >= 2]
        hemi = [v for v in vs if v["ploidy"] == 1 and v["copies"] == 1]          # male chrX/Y haploid calls
        het = [v for v in vs if v["copies"] == 1 and v["ploidy"] >= 2]
        phase = "-"
        if hom:
            geno, alleles = "homozygous", [hom[0], hom[0]]
        elif hemi:
            geno, alleles = "hemizygous", [hemi[0], hemi[0]]
        elif len(het) >= 2:
            # best pair that is not in cis: try pairs in order of combined score
            def pair_phase(p, q):
                hp, hq = hap(p["gt"]), hap(q["gt"])
                if p.get("over_del") is q or q.get("over_del") is p:
                    return "trans"
                if p.get("is_sv") or q.get("is_sv"):
                    return "unknown"
                if p["chrom"] == q["chrom"] and p["pos"] == q["pos"]:
                    return "trans"
                if hp is not None and hq is not None and p["ps"] and p["ps"] == q["ps"]:
                    return "trans" if hp != hq else "cis"
                return "unknown"
            pairs = sorted(((het[i], het[j]) for i in range(len(het)) for j in range(i + 1, len(het))),
                           key=lambda pq: -(pq[0]["score"] + pq[1]["score"]))
            v1, v2 = next((pq for pq in pairs if pair_phase(*pq) != "cis"), pairs[0])
            h1, h2 = hap(v1["gt"]), hap(v2["gt"])
            if v1.get("over_del") is v2 or v2.get("over_del") is v1:
                phase = "trans (SNV hemizygous within the deletion)"
            elif v1.get("is_sv") or v2.get("is_sv"):
                phase = "unphased (SV + small variant)"
            elif v1["chrom"] == v2["chrom"] and v1["pos"] == v2["pos"]:
                phase = "trans (two ALT alleles at one site)"      # split from a 1/2 genotype
            elif h1 is not None and h2 is not None and v1["ps"] and v1["ps"] == v2["ps"]:
                phase = "trans (phased)" if h1 != h2 else "cis (phased)"
            else:
                phase = "unphased"
            geno = "compound_het?" if phase != "cis (phased)" else "two_het_in_cis"
            alleles = [v1, v2]
        else:
            geno, alleles = "monoallelic", [vs[0]]
        mito = vs[0]["chrom"] in ("chrM", "MT", "M", "chrMT")
        if mito:                       # heteroplasmy: recessive (biallelic) logic does not apply
            geno, phase = "mitochondrial", "-"
        biallelic = geno in ("homozygous", "hemizygous", "compound_het?")
        strong = [x["class"] in STRONG and x["score"] >= 4 for x in alleles]
        has_plp = any(x["class"] == "CLINVAR_PLP" for x in alleles)
        gs = vs[0]["gene_set"]
        mva = gs in ("MVA_established", "MVA_emerging", "panel")
        if mva and biallelic and has_plp and all(x["score"] >= 3 for x in alleles):
            tier = 1
        elif gs in ("MVA_established", "panel") and biallelic and all(strong):
            tier = 2
        elif mva and biallelic:
            tier = 3
        elif mva and alleles[0]["score"] >= 4:
            tier = 4
        elif biallelic and all(strong):
            tier = 5
        else:
            tier = 6
        gscore = sum(x["score"] for x in alleles) + (2 if biallelic else 0)
        note = []
        if phase == "unphased":
            note.append("confirm in trans (parents or long reads)")
        elif phase.startswith("unphased (SV"):
            note.append("SV-small variant phase unknown: confirm with parents or long reads")
        if mito:
            note.append(f"mitochondrial: interpret by heteroplasmy level (VAF {vs[0]['vaf']}), not biallelic logic")
        elif tier == 4:
            note.append("second hit? check SV/CNV calls, RNA-seq, non-coding")
        fq = sorted({f"{v['chrom']}:{v['pos']}" for v in vs if not v.get("is_sv") and
                     ((v["chrom"], v["pos"]) in failed or v["filter"] not in ("PASS", "."))})
        if fq:
            note.append("includes ClinVar P/LP that FAILED QC at " + ",".join(fq) + ": confirm (IGV/Sanger)")
        gtn = sorted({p for x in alleles for p in (x.get("gt_note") or "").split("; ") if p.startswith("GT differs")})
        if gtn:
            note.append("genotype differs between the two VCFs (" + "; ".join(gtn) + "): confirm")
        if any((x.get("gt_note") or "").startswith("representation differs") for x in alleles):
            note.append("indel written differently by the two callers (counted once; see ranked_variants.tsv)")
        one = sorted({x["source"] for x in alleles if x.get("source") not in ("", "both", None)})
        if one:
            note.append("allele(s) found in only one VCF (" + ",".join(one) + "): confirm")
        if any(x.get("is_sv") and x["filter"] != "PASS" for x in alleles):
            note.append("includes an SV that failed caller filters: confirm (IGV/split reads/depth)")
        a0 = alleles[0]
        if geno == "monoallelic" and a0.get("is_sv") and a0.get("gt_source") == "placeholder":
            dr = fnum(a0.get("depth_ratio"))
            if dr is None:
                note.append("SV genotype not measured (AnnotSV placeholder GT): if homozygous this gene is biallelic - "
                            "give --bam or check IGV")
            elif (a0["svtype"] in ("DEL", "CNV") and dr >= 0.85) or (a0["svtype"] == "DUP" and dr < 1.15):
                note.append(f"SV not supported by read depth ({dr:.2f}x of flanks): probably a false call - check IGV")
            else:
                note.append(f"SV depth {dr:.2f}x of flanks: partial change (mosaic?) - genotype unresolved, check IGV")
        if geno == "homozygous" and a0.get("gt_source") == "depth":
            note.append(f"homozygous {'deletion' if a0.get('svtype') in ('DEL', 'CNV') else 'SV'} called from read depth"
                        f"{' (' + a0['depth_ratio'] + 'x of flanks)' if a0.get('depth_ratio') not in (None, '-') else ''}"
                        ": confirm (MLPA/qPCR, breakpoint PCR)")
        if geno == "monoallelic" and alleles[0].get("is_sv") and \
           (("loss-sensitive" in alleles[0].get("gene_flag", "") and alleles[0]["svtype"] in ("DEL", "CNV")) or
                ("gain-sensitive" in alleles[0].get("gene_flag", "") and alleles[0]["svtype"] == "DUP")):
            note.append("dosage-sensitive gene (HI/TS/pLI/AD): one SV copy may be enough - review a dominant model")
        if any(x.get("is_sv") and x.get("confidence") == "low" for x in alleles):
            note.append("SV breakpoint in a gap/blacklist/segdup region: confirm (IGV, long reads)")
        if any(x["class"].startswith("AG_") for x in alleles):
            note.append("allele supported only by AlphaGenome prediction: test RNA (RT-PCR / RNA-seq) or a reporter assay")
        if any(x.get("callers") == "depth" or x["class"] == "CNV_DEPTH" for x in alleles):
            note.append("depth-only CNV: confirm by MLPA/qPCR")
        single = int(any(x.get("source") not in ("", "both", None) for x in alleles))
        genes.append({"_single": single, "tier": tier, "gene": g, "gene_set": gs, "genotype": geno, "phase": phase,
                      "gene_score": gscore, "exomiser": f"{exo[g]:.3f}" if g in exo else "-",
                      "n_qualifying": len(vs),
                      "variants": " | ".join(
                          (f"{x['chrom']}:{x['pos']}-{x['sv_end']} {x['alt']} {x['gt']} {x['class']} {x['consequence']}"
                           if x.get("is_sv") else
                           f"{x['chrom']}:{x['pos']} {x['ref']}>{x['alt']} {x['gt']} {x['class']} "
                           f"{x['hgvsp'] or x['hgvsc'] or x['consequence']}")
                          for x in (alleles[:1] if geno in ("homozygous", "hemizygous") else alleles)),
                      "note": "; ".join(note) or "-"})
    # failed-QC ClinVar variants in MVA genes with no passing variant still deserve a row
    # within a tier: alleles seen in both VCFs first, then Exomiser, then gene score
    genes.sort(key=lambda r: (r["tier"], r.pop("_single", 0), -(exo.get(r["gene"], 0.0)), -r["gene_score"], r["gene"]))

    # ---------------------------------------------------------------- write
    vcols = ["gene", "gene_set", "chrom", "pos", "id", "ref", "alt", "gt", "class", "score", "consequence",
             "impact", "hgvsc", "hgvsp", "max_af", "clinvar", "dp", "vaf", "evidence"]
    gtier = {r["gene"]: r["tier"] for r in genes}
    variants.sort(key=lambda v: (gtier[v["gene"]], -v["score"], v["gene"]))
    with open(os.path.join(a.out, "ranked_variants.tsv"), "w") as fh:
        fh.write("tier\t" + "\t".join(vcols) + "\n")
        for v in variants:
            fh.write(f"{gtier[v['gene']]}\t" + "\t".join(str(v[c]) for c in vcols) + "\n")
    gcols = ["tier", "gene", "gene_set", "genotype", "phase", "gene_score", "exomiser", "n_qualifying",
             "variants", "note"]
    with open(os.path.join(a.out, "ranked_genes.tsv"), "w") as fh:
        fh.write("\t".join(gcols) + "\n")
        for r in genes:
            fh.write("\t".join(str(r[c]) for c in gcols) + "\n")
    # top variants (small variants and SVs together) in MVA genes: gene tier first, then the variant's own score
    ginfo = {r["gene"]: r for r in genes}
    mva_v = [v for v in variants if v["gene_set"] in ("MVA_established", "MVA_emerging")]
    mva_v.sort(key=lambda v: (gtier[v["gene"]], -v["score"], v["gene"], int(v["pos"])))
    with open(os.path.join(a.out, "top_mva_variants.tsv"), "w") as fh:
        fh.write("rank\ttier\tgene\tgene_set\ttype\tvariant\tgt\tclass\tscore\tconsequence\thgvsp_or_effect\t"
                 "max_af\tclinvar\tgene_genotype\tgene_phase\tevidence\tgene_note\n")
        for i, v in enumerate(mva_v[:a.top_mva], 1):
            sv = v.get("is_sv")
            vt = "SV" if sv else ("SNV" if len(v["ref"]) == 1 and len(v["alt"]) == 1 else "indel")
            loc = (f"{v['chrom']}:{v['pos']}-{v['sv_end']} {v['alt']}" if sv else f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}")
            gi = ginfo[v["gene"]]
            fh.write("\t".join(str(x) for x in (
                i, gtier[v["gene"]], v["gene"], v["gene_set"], vt, loc, v["gt"], v["class"], v["score"], v["consequence"],
                (v["hgvsp"] or v["hgvsc"] or "-") if not sv else v["consequence"], v["max_af"] or "-", v["clinvar"] or "-",
                gi["genotype"], gi["phase"], v["evidence"] or "-", gi["note"])) + "\n")
    print(f"[rank] MVA genes: {len(mva_v)} qualifying variant(s) ({sum(1 for v in mva_v if v.get('is_sv'))} SV) -> "
          f"top {min(a.top_mva, len(mva_v))} in top_mva_variants.tsv")

    targets = [r for r in genes if r["tier"] <= 4] + [r for r in genes if r["tier"] == 5][:a.top_discovery]
    with open(os.path.join(a.out, "targets.tsv"), "w") as fh:
        fh.write("\t".join(gcols) + "\n")
        for r in targets:
            fh.write("\t".join(str(r[c]) for c in gcols) + "\n")
    tg = {r["gene"] for r in targets}
    with open(os.path.join(a.out, "targets.bed"), "w") as fh:
        for v in variants:
            if v["gene"] in tg:
                p = int(v["pos"])
                e = v.get("sv_end", p)
                fh.write(f"{v['chrom']}\t{max(0, p - 51)}\t{e + 50}\t{v['gene']}_{v['class']}\n")

    tiers = defaultdict(int)
    for r in genes:
        tiers[r["tier"]] += 1
    print(f"[rank] {n_rows:,} table rows -> {len(variants):,} qualifying rare variants in {len(genes):,} genes")
    print("[rank] genes per tier: " + ", ".join(f"T{t}={tiers[t]}" for t in sorted(tiers)))
    print(f"[rank] targets: {len(targets)} genes -> {os.path.join(a.out, 'targets.tsv')}")
    print("[rank] excluded rows: " + "; ".join(f"{k}: {v:,}" for k, v in sorted(dropped.items(), key=lambda kv: -kv[1])))
    missing = [c for c in ("REVEL", "SpliceAI_pred_DS_AG", "CADD_PHRED", "am_class") if c not in first_cols]
    if "REVEL" in missing:
        print("[rank] note: no REVEL - missense evidence uses AlphaMissense (likely_pathogenic -> moderate, likely_benign -> score 1) "
              "and SIFT+PolyPhen agreement (-> supporting)")
    if missing:
        print("[rank] note: no " + ", ".join(missing) + " columns - those evidence classes are unavailable "
              "(add the VEP plugins to use them)")


if __name__ == "__main__":
    main()
