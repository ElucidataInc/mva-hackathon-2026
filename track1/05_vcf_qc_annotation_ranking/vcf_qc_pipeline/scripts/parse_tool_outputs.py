#!/usr/bin/env python3
"""Parse the outputs of external tools (bcftools, vcftools, GATK, somalier, peddy,
hap.py, RTG vcfeval) found in a pipeline output directory into one JSON file."""

import argparse
import csv
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcf_common import write_json


def _num(x):
    try:
        if re.fullmatch(r"-?\d+", x):
            return int(x)
        return float(x)
    except (ValueError, TypeError):
        return x


def parse_bcftools_stats(path):
    """bcftools stats -> {section: [rows]} (column names taken from the '# XX\t[2]id...' lines)."""
    if not os.path.exists(path):
        return None
    cols, data = {}, {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("# ") and "\t" in line:
                p = line[2:].split("\t")
                sec = p[0]
                cols[sec] = [re.sub(r"^\[\d+\]", "", c) for c in p]
            elif line and not line.startswith("#"):
                p = line.split("\t")
                sec = p[0]
                names = cols.get(sec) or [f"c{i}" for i in range(len(p))]
                data.setdefault(sec, []).append({names[i] if i < len(names) else f"c{i}": _num(v)
                                                 for i, v in enumerate(p)})
    out = {"sections": list(data)}
    for sec in ("SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP", "PSC", "PSI", "HWE", "ID"):
        if sec in data:
            out[sec] = data[sec]
    return out


def read_table(path, delim="\t", limit=None):
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("##")]
    rdr = csv.DictReader(lines, delimiter=delim)
    rows = []
    for i, r in enumerate(rdr):
        rows.append({(k or "").strip().lstrip("#"): _num((v or "").strip()) for k, v in r.items()})
        if limit and i + 1 >= limit:
            break
    return rows


def parse_gatk_metrics(path):
    if not os.path.exists(path):
        return None
    lines = open(path).read().splitlines()
    out, i = [], 0
    while i < len(lines):
        if lines[i].startswith("## METRICS CLASS"):
            hdr = lines[i + 1].split("\t")
            j = i + 2
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("#"):
                out.append(dict(zip(hdr, [_num(x) for x in lines[j].split("\t")])))
                j += 1
            i = j
        i += 1
    return out


def summarise_vcftools(d):
    res = {}
    imiss = read_table(os.path.join(d, "vcftools.imiss"))
    if imiss:
        res["imiss"] = imiss
    het = read_table(os.path.join(d, "vcftools.het"))
    if het:
        res["het"] = het
    idepth = read_table(os.path.join(d, "vcftools.idepth"))
    if idepth:
        res["idepth"] = idepth
    tstv = read_table(os.path.join(d, "vcftools.TsTv.summary"))
    if tstv:
        res["tstv_summary"] = tstv
    lmiss = os.path.join(d, "vcftools.lmiss")
    if os.path.exists(lmiss):
        bins = [0] * 21
        n = 0
        with open(lmiss) as fh:
            next(fh)
            for line in fh:
                try:
                    f = float(line.rstrip().split("\t")[-1])
                except ValueError:
                    continue
                bins[min(20, int(round(f * 20)))] += 1
                n += 1
        res["lmiss_hist"] = bins
        res["lmiss_n"] = n
    hwe = os.path.join(d, "vcftools.hwe")
    if os.path.exists(hwe):
        n = n6 = n4 = 0
        with open(hwe) as fh:
            hdr = next(fh).rstrip().split("\t")
            k = hdr.index("P_HWE") if "P_HWE" in hdr else 5
            for line in fh:
                p = line.rstrip().split("\t")
                try:
                    pv = float(p[k])
                except (ValueError, IndexError):
                    continue
                n += 1
                n6 += pv < 1e-6
                n4 += pv < 1e-4
        res["hwe"] = {"sites_tested": n, "p_lt_1e-4": n4, "p_lt_1e-6": n6}
    ld = os.path.join(d, "vcftools.ldepth.mean")
    if os.path.exists(ld):
        vals = []
        with open(ld) as fh:
            next(fh)
            for line in fh:
                try:
                    vals.append(float(line.split("\t")[2]))
                except (ValueError, IndexError):
                    pass
        if vals:
            vals.sort()
            res["site_mean_depth"] = {"n": len(vals), "median": vals[len(vals) // 2], "mean": sum(vals) / len(vals)}
    return res or None


def parse_norm_log(path):
    if not os.path.exists(path):
        return None
    txt = open(path).read()
    res = {"ref_mismatch_lines": len(re.findall(r"REF_MISMATCH", txt)),
           "warnings": len(re.findall(r"(?i)warn", txt))}
    m = re.search(r"Lines\s+total/split/(?:joined/)?realigned/(?:mismatch_removed/)?(?:dup_removed/)?skipped:\s*([\d/]+)", txt)
    if m:
        res["norm_summary"] = m.group(0)
    res["examples"] = [l for l in txt.splitlines() if "REF_MISMATCH" in l][:15]
    return res


def parse_validate(path):
    if not os.path.exists(path):
        return None
    txt = open(path).read()
    errs = [l for l in txt.splitlines() if re.search(r"(?i)\b(error|exception|invalid|fail)", l)]
    warns = [l for l in txt.splitlines() if "WARN" in l]
    ok = "Shutting down engine" in txt and not re.search(r"A USER ERROR has occurred", txt)
    return {"passed": ok and not re.search(r"(?i)validation.*fail", txt), "n_error_lines": len(errs),
            "n_warn_lines": len(warns), "errors": errs[:20], "warnings": warns[:20]}


def parse_happy(d):
    f = sorted(glob.glob(os.path.join(d, "*.summary.csv")))
    return read_table(f[0], ",") if f else None


def parse_vcfeval(d):
    f = os.path.join(d, "summary.txt")
    if not os.path.exists(f):
        return None
    lines = [l for l in open(f).read().splitlines() if l.strip() and not set(l.strip()) <= {"-"}]
    hdr = lines[0].split()
    return [dict(zip(hdr, [_num(x) for x in l.split()])) for l in lines[1:]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="pipeline output directory")
    ap.add_argument("--out", required=True, help="output JSON")
    a = ap.parse_args()
    o = a.outdir
    res = {}
    res["bcftools_stats"] = parse_bcftools_stats(os.path.join(o, "03_bcftools_stats", "bcftools.stats.txt"))
    res["bcftools_stats_pass"] = parse_bcftools_stats(os.path.join(o, "03_bcftools_stats", "bcftools.stats.PASS.txt"))
    res["norm_check"] = parse_norm_log(os.path.join(o, "02_integrity", "norm_checkref.log"))
    res["validate_variants"] = parse_validate(os.path.join(o, "02_integrity", "gatk_ValidateVariants.log"))
    res["bcftools_index_check"] = (open(os.path.join(o, "02_integrity", "index_check.txt")).read()
                                   if os.path.exists(os.path.join(o, "02_integrity", "index_check.txt")) else None)
    sq = os.path.join(o, "05_sample_qc")
    res["vcftools"] = summarise_vcftools(os.path.join(sq, "vcftools"))
    res["guess_ploidy"] = read_table(os.path.join(sq, "guess_ploidy.tsv"))
    cp = os.path.join(sq, "check_ploidy.txt")
    res["check_ploidy"] = open(cp).read()[:5000] if os.path.exists(cp) else None
    md = os.path.join(sq, "mendelian2.txt")
    res["mendelian2"] = open(md).read()[:8000] if os.path.exists(md) else None
    res["somalier_samples"] = read_table(os.path.join(sq, "somalier", "somalier.samples.tsv"))
    res["somalier_pairs"] = read_table(os.path.join(sq, "somalier", "somalier.pairs.tsv"), limit=500)
    anc = glob.glob(os.path.join(sq, "somalier", "*ancestry*.tsv"))
    res["somalier_ancestry"] = read_table(anc[0]) if anc else None
    pd_ = os.path.join(sq, "peddy")
    for k in ("sex_check", "het_check", "ped_check"):
        f = glob.glob(os.path.join(pd_, f"*.{k}.csv"))
        res[f"peddy_{k}"] = read_table(f[0], ",", limit=2000) if f else None
    res["happy"] = parse_happy(os.path.join(o, "07_concordance", "happy"))
    res["vcfeval"] = parse_vcfeval(os.path.join(o, "07_concordance", "vcfeval"))
    g = os.path.join(o, "08_gatk_metrics")
    res["gatk_cvcm_summary"] = parse_gatk_metrics(os.path.join(g, "cvcm.variant_calling_summary_metrics"))
    res["gatk_cvcm_detail"] = parse_gatk_metrics(os.path.join(g, "cvcm.variant_calling_detail_metrics"))
    an = os.path.join(o, "06_annotation", "annotation_summary.json")
    res["external_annotation"] = json.load(open(an)) if os.path.exists(an) else None
    write_json({k: v for k, v in res.items() if v}, a.out)


if __name__ == "__main__":
    main()
