#!/usr/bin/env python3
"""Generate a small synthetic data set to exercise every part of the pipeline:
reference FASTA, a 12-sample germline VCF (trio, duplicate pair, contaminated sample,
sex mismatch, consanguineous sample, phasing, missingness, multi-allelics, deliberate
sort / duplicate / REF errors), a sex file, a PED and a Mutect2-style tumor/normal VCF."""

import os
import random
import sys

out = sys.argv[1] if len(sys.argv) > 1 else "test_data"
os.makedirs(out, exist_ok=True)
rng = random.Random(42)
B = "ACGT"
CHROMS = [("chr1", 6_000_000), ("chr2", 5_000_000), ("chr3", 4_000_000), ("chrX", 4_000_000), ("chrY", 1_000_000)]

# ------------------------------------------------------------------ reference
seqs = {}
with open(f"{out}/ref.fa", "w") as fa:
    for c, L in CHROMS:
        s = "".join(rng.choice(B) for _ in range(L))
        seqs[c] = s
        fa.write(f">{c}\n")
        for i in range(0, L, 60):
            fa.write(s[i:i + 60] + "\n")

S = [f"S{i:02d}" for i in range(1, 13)]
# S01 child of S02 (father) x S03 (mother); S12 = duplicate of S05; S07 contaminated;
# S08 declared F but genetically male; S09 consanguineous (long ROH)
SEX = {"S01": "M", "S02": "M", "S03": "F", "S04": "F", "S05": "M", "S06": "F", "S07": "F", "S08": "F",
       "S09": "M", "S10": "F", "S11": "M", "S12": "M"}
TRUE_SEX = dict(SEX, S08="M")
with open(f"{out}/sex.tsv", "w") as fh:
    for s in S:
        fh.write(f"{s}\t{SEX[s]}\n")
with open(f"{out}/family.ped", "w") as fh:
    for s in S:
        pat, mat = ("S02", "S03") if s == "S01" else ("0", "0")
        fh.write(f"FAM1\t{s}\t{pat}\t{mat}\t{1 if SEX[s] == 'M' else 2}\t1\n")


def gt_from_af(p):
    return (rng.random() < p) + (rng.random() < p)


def fmt_gt(d, phased=False, hap=False):
    if hap:
        return str(min(d, 1))
    sep = "|" if phased else "/"
    return {0: f"0{sep}0", 1: (f"0{sep}1" if rng.random() < .5 or not phased else f"1{sep}0"), 2: f"1{sep}1"}[d]


header = ["##fileformat=VCFv4.2", "##fileDate=20261003",
          '##source=HaplotypeCaller', '##GATKCommandLine=<ID=GenotypeGVCFs,CommandLine="GenotypeGVCFs -R ref.fa -V gendb://db -O cohort.vcf.gz",Version="4.5.0.0">',
          '##FILTER=<ID=PASS,Description="All filters passed">',
          '##FILTER=<ID=LowQual,Description="Low quality">',
          '##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="Truth sensitivity tranche level for SNP model">',
          '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
          '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">',
          '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
          '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
          '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">',
          '##INFO=<ID=QD,Number=1,Type=Float,Description="QualByDepth">',
          '##INFO=<ID=FS,Number=1,Type=Float,Description="FisherStrand">',
          '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">',
          '##INFO=<ID=SOR,Number=1,Type=Float,Description="StrandOddsRatio">',
          '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|BIOTYPE|gnomADg_AF|CLIN_SIG">',
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
          '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
          '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
          '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">']
header += [f"##contig=<ID={c},length={L}>" for c, L in CHROMS]
header.append("##reference=file:///data/ref/GRCh38_test.fa")
header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(S))

CONS = [("missense_variant", "MODERATE"), ("synonymous_variant", "LOW"), ("intron_variant", "MODIFIER"),
        ("intergenic_variant", "MODIFIER"), ("stop_gained", "HIGH"), ("frameshift_variant", "HIGH"),
        ("splice_region_variant", "LOW"), ("upstream_gene_variant", "MODIFIER")]
GENES = ["BRCA1", "TP53", "EGFR", "KRAS", "MYC", "APOE", "CFTR", "LDLR", "PTEN", "ATM"]
lines = []
dup_done = unsorted_done = 0
roh_region = ("chr2", 1_000_000, 4_000_000)

for c, L in CHROMS:
    n = {"chr1": 9000, "chr2": 7000, "chr3": 6000, "chrX": 3000, "chrY": 300}[c]
    # non-uniform density: a dense hotspot on each chromosome
    poss = set()
    while len(poss) < n:
        if rng.random() < 0.15:
            p = int(rng.gauss(L * 0.3, L * 0.02))
        else:
            p = rng.randint(100, L - 100)
        if 100 < p < L - 100:
            poss.add(p)
    for pos in sorted(poss):
        ref = seqs[c][pos - 1]
        r = rng.random()
        if r < 0.82:
            alt = rng.choice([b for b in B if b != ref])
            # bias toward transitions (Ts/Tv ~2)
            if rng.random() < 0.67:
                alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]
            alts = [alt]
        elif r < 0.90:
            k = min(int(rng.expovariate(0.5)) + 1, 40)
            alts = [ref + "".join(rng.choice(B) for _ in range(k))]
        elif r < 0.97:
            k = min(int(rng.expovariate(0.5)) + 1, 40)
            ref = seqs[c][pos - 1: pos + k]
            alts = [ref[0]]
        elif r < 0.985:
            ref = seqs[c][pos - 1:pos + 1]
            alts = ["".join(rng.choice(B) for _ in range(2))]
            if alts[0] == ref:
                continue
        else:
            a1 = rng.choice([b for b in B if b != ref])
            alts = [a1, ref + "T"]
        af = min(0.95, max(0.003, rng.betavariate(0.4, 1.6)))
        dp_site = 0
        samp = []
        ac = 0
        an = 0
        phased_site = rng.random() < 0.6
        for s in S:
            is_male = TRUE_SEX[s] == "M"
            hap = (c in ("chrX", "chrY")) and is_male and not (c == "chrX" and pos < 2_700_000 and False)
            if c == "chrY" and not is_male:
                samp.append("./.:0,0:0:.:.")
                continue
            d = gt_from_af(af)
            if s == "S09" and c == roh_region[0] and roh_region[1] < pos < roh_region[2]:
                d = 2 if d >= 1 else 0
            if c == "chrX" and is_male:
                d = 2 if d >= 1 else 0
            samp.append(d)
        # S12 duplicate of S05 (with a little noise); S01 child of S02/S03
        idx = {s: i for i, s in enumerate(S)}
        if not isinstance(samp[idx["S05"]], str):
            samp[idx["S12"]] = samp[idx["S05"]] if rng.random() > 0.002 else gt_from_af(af)
        f_, m_ = samp[idx["S02"]], samp[idx["S03"]]
        if c not in ("chrX", "chrY") and not isinstance(f_, str) and not isinstance(m_, str):
            tf = (rng.random() < f_ / 2) if f_ != 1 else rng.random() < .5
            tm = (rng.random() < m_ / 2) if m_ != 1 else rng.random() < .5
            ch = int(tf) + int(tm)
            if rng.random() < 0.001:
                ch = (ch + 1) % 3
            samp[idx["S01"]] = ch
        out_s = []
        for s, d in zip(S, samp):
            if isinstance(d, str):
                out_s.append(d)
                an += 0
                continue
            if rng.random() < (0.12 if s == "S10" else 0.01):
                out_s.append("./.:.:.:.:.")
                continue
            hap = c in ("chrX", "chrY") and TRUE_SEX[s] == "M" and not (c == "chrX" and pos < 2_781_479)
            dp = max(3, int(rng.gauss(32 if s != "S11" else 14, 7)))
            dp_site += dp
            if d == 0:
                altc = 0 if rng.random() > 0.03 else 1
            elif d == 1:
                p = 0.5 if s != "S07" else rng.choice([0.2, 0.5, 0.8, 0.3])
                altc = sum(rng.random() < p for _ in range(dp))
            else:
                altc = dp - (0 if rng.random() > 0.05 else 1)
            if s == "S07" and d == 0 and rng.random() < 0.08:
                d = 1
                altc = max(1, int(dp * 0.12))
            gq = min(99, int(rng.gauss(70, 20))) if d != 1 else min(99, int(rng.gauss(55, 20)))
            gq = max(0, gq)
            gt = fmt_gt(d, phased=phased_site and d == 1, hap=hap)
            if len(alts) == 2 and d >= 1 and rng.random() < 0.3:
                gt = "1/2" if d == 1 else "2/2"
                altc_s = f"{dp - altc},{altc // 2},{altc - altc // 2}"
            else:
                altc_s = f"{dp - altc},{altc}" + (",0" if len(alts) == 2 else "")
            ps = str(pos - pos % 50000) if "|" in gt else "."
            ac += d if not hap else min(d, 1)
            an += 1 if hap else 2
            out_s.append(f"{gt}:{altc_s}:{dp}:{gq}:{ps}")
        if ac == 0 and rng.random() < 0.7:
            continue
        qual = round(max(5, rng.gauss(400, 300)), 2)
        flt = "PASS"
        if qual < 50:
            flt = "LowQual"
        elif rng.random() < 0.05:
            flt = "VQSRTrancheSNP99.90to100.00"
        rid = f"rs{rng.randint(1000, 99999999)}" if rng.random() < (0.85 if af > 0.02 else 0.3) else "."
        info = [f"AC={','.join([str(ac)] * len(alts))}", f"AN={an}", f"AF={round(af, 4)}", f"DP={dp_site}",
                f"QD={round(rng.gauss(18, 6), 2)}", f"FS={round(abs(rng.gauss(2, 5)), 3)}",
                f"MQ={round(min(70, rng.gauss(59, 2)), 2)}", f"SOR={round(abs(rng.gauss(0.8, 0.4)), 3)}"]
        if rid != ".":
            info.insert(4, "DB")
        cons, imp = rng.choice(CONS)
        if len(alts[0]) != len(ref) and cons in ("missense_variant", "synonymous_variant", "stop_gained"):
            cons, imp = "frameshift_variant", "HIGH"
        gn = f"{af * rng.uniform(.5, 1.5):.5f}" if rid != "." else ""
        clin = rng.choice(["", "", "", "", "benign", "likely_benign", "uncertain_significance", "pathogenic"])
        info.append("CSQ=" + ",".join(f"{a}|{cons}|{imp}|{rng.choice(GENES)}|ENSG0000{rng.randint(1, 9)}|protein_coding|{gn}|{clin}"
                                      for a in alts))
        ref_out = ref
        if rng.random() < 0.0008 and len(ref) == 1:
            ref_out = rng.choice([b for b in B if b != ref])     # deliberate REF mismatch
        line = "\t".join([c, str(pos), rid, ref_out, ",".join(alts), str(qual), flt, ";".join(info),
                          "GT:AD:DP:GQ:PS"] + out_s)
        lines.append(line)
        if c == "chr1" and not dup_done and pos > 3_000_000:
            lines.append(line)                                   # deliberate duplicate
            dup_done = 1
# one unsorted record
i = len(lines) // 3
lines[i], lines[i + 1] = lines[i + 1], lines[i]

with open(f"{out}/cohort.vcf", "w") as fh:
    fh.write("\n".join(header) + "\n")
    fh.write("\n".join(lines) + "\n")

# ------------------------------------------------------------------ dbSNP-like / gnomAD-like resources
seen = set()
known_lines = []
for l in sorted(lines, key=lambda x: ([c for c, _ in CHROMS].index(x.split("\t")[0]), int(x.split("\t")[1]))):
    f = l.split("\t")
    if (f[0], f[1]) in seen or rng.random() > 0.75:
        continue
    seen.add((f[0], f[1]))
    rid = f[2] if f[2] != "." else f"rs{rng.randint(10 ** 8, 10 ** 9)}"
    known_lines.append("\t".join([f[0], f[1], rid, f[3], f[4], ".", ".", "AF=" + ",".join(f"{rng.uniform(1e-5, 0.6):.6f}" for _ in f[4].split(","))]))
with open(f"{out}/known_sites.vcf", "w") as fh:
    fh.write("##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID={c},length={L}>\n" for c, L in CHROMS) +
             '##INFO=<ID=AF,Number=A,Type=Float,Description="Population allele frequency">\n'
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(known_lines) + "\n")

# ------------------------------------------------------------------ somatic tumor/normal
sh = ["##fileformat=VCFv4.2", '##source=Mutect2', "##tumor_sample=TUMOR", "##normal_sample=NORMAL",
      '##GATKCommandLine=<ID=Mutect2,CommandLine="Mutect2 -R ref.fa -I tumor.bam -I normal.bam -normal NORMAL --germline-resource af-only-gnomad.vcf.gz --panel-of-normals pon.vcf.gz -O somatic.vcf.gz",Version="4.5.0.0">',
      '##GATKCommandLine=<ID=FilterMutectCalls,CommandLine="FilterMutectCalls -V somatic.vcf.gz -O filtered.vcf.gz",Version="4.5.0.0">',
      '##FILTER=<ID=PASS,Description="All filters passed">',
      '##FILTER=<ID=panel_of_normals,Description="Blacklisted site in panel of normals">',
      '##FILTER=<ID=germline,Description="Evidence indicates this site is germline, not somatic">',
      '##FILTER=<ID=strand_bias,Description="Evidence for alt allele comes from one read direction only">',
      '##FILTER=<ID=orientation,Description="orientation bias detected by the orientation bias mixture model">',
      '##FILTER=<ID=weak_evidence,Description="Mutation does not meet likelihood threshold">',
      '##FILTER=<ID=normal_artifact,Description="artifact_in_normal">',
      '##INFO=<ID=TLOD,Number=A,Type=Float,Description="Log 10 likelihood ratio score of variant existing versus not existing">',
      '##INFO=<ID=NLOD,Number=A,Type=Float,Description="Normal log 10 likelihood ratio of diploid het or hom alt genotypes">',
      '##INFO=<ID=NALOD,Number=A,Type=Float,Description="Negative log 10 odds of artifact in normal">',
      '##INFO=<ID=GERMQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not germline">',
      '##INFO=<ID=POPAF,Number=A,Type=Float,Description="negative log 10 population allele frequencies of alt alleles">',
      '##INFO=<ID=STRANDQ,Number=1,Type=Integer,Description="Phred-scaled quality of strand bias artifact">',
      '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
      '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions">',
      '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
      '##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation">',
      '##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation">',
      '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics for strand bias">']
sh += [f"##contig=<ID={c},length={L}>" for c, L in CHROMS]
sh.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL")
sl = []
for c, L in CHROMS[:4]:
    for pos in sorted(rng.sample(range(1000, L - 1000), 700)):
        ref = seqs[c][pos - 1]
        alt = rng.choice([b for b in B if b != ref]) if rng.random() < .9 else ref + "A"
        kind = rng.random()
        tdp, ndp = max(10, int(rng.gauss(80, 20))), max(8, int(rng.gauss(40, 10)))
        if kind < 0.12:
            vaf_t, vaf_n, flt = 0.5, 0.5, "germline"
        elif kind < 0.22:
            vaf_t, vaf_n, flt = rng.uniform(.02, .1), 0.02, "panel_of_normals"
        elif kind < 0.3:
            vaf_t, vaf_n, flt = rng.uniform(.03, .2), 0, "orientation"
        elif kind < 0.36:
            vaf_t, vaf_n, flt = rng.uniform(.03, .2), 0, "strand_bias"
        elif kind < 0.5:
            vaf_t, vaf_n, flt = rng.uniform(.01, .05), 0, "weak_evidence"
        else:
            vaf_t, vaf_n, flt = min(.95, max(.03, rng.betavariate(2, 5))), 0, "PASS"
        ta = sum(rng.random() < vaf_t for _ in range(tdp))
        na = sum(rng.random() < vaf_n for _ in range(ndp))
        ta = max(ta, 2)
        if flt == "orientation":
            f1 = ta
            f2 = 0
        else:
            f1 = sum(rng.random() < .5 for _ in range(ta))
            f2 = ta - f1
        if flt == "strand_bias":
            sbf, sbr = ta, 0
        else:
            sbf = sum(rng.random() < .5 for _ in range(ta))
            sbr = ta - sbf
        tlod = round(ta * 2.9 + rng.gauss(0, 2), 2)
        nlod = round(max(0, (ndp - na) * 0.3 - na * 3 + rng.gauss(0, 1)), 2)
        rr = tdp - ta
        info = f"DP={tdp + ndp};TLOD={tlod};NLOD={nlod};NALOD={round(rng.gauss(1.5, 0.6), 3)};GERMQ={rng.randint(1, 93)};" \
               f"POPAF={round(rng.uniform(1, 7), 2)};STRANDQ={rng.randint(1, 93)}"
        t = f"0/1:{rr},{ta}:{ta / tdp:.3f}:{tdp}:{rr // 2},{f1}:{rr - rr // 2},{f2}:{rr // 2},{rr - rr // 2},{sbf},{sbr}"
        ngt = "0/1" if na > ndp * .3 else "0/0"
        n = f"{ngt}:{ndp - na},{na}:{na / ndp:.3f}:{ndp}:{(ndp - na) // 2},{na // 2}:{(ndp - na) - (ndp - na) // 2},{na - na // 2}:" \
            f"{(ndp - na) // 2},{(ndp - na) - (ndp - na) // 2},{na // 2},{na - na // 2}"
        sl.append("\t".join([c, str(pos), ".", ref, alt, ".", flt, info, "GT:AD:AF:DP:F1R2:F2R1:SB", t, n]))
with open(f"{out}/somatic.vcf", "w") as fh:
    fh.write("\n".join(sh) + "\n" + "\n".join(sl) + "\n")
print(f"wrote test data to {out}: {len(lines)} germline records, {len(sl)} somatic records")
