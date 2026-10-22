#!/usr/bin/env python3
"""Summarise an externally annotated VCF (dbSNP IDs, gnomAD AF, VEP/SnpEff, ClinVar):
known-vs-novel rate overall and per variant type, plus consequence / impact / gene tables."""

import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcf_common import classify_allele, detect_annotations, open_vcf, parse_header, parse_info, to_float, write_json
from vcf_profile import AnnotationCounter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gnomad-tag", default="gnomAD_AF", help="INFO tag holding the gnomAD AF")
    ap.add_argument("--dbsnp-tag", default="", help="optional INFO flag added for dbSNP membership")
    a = ap.parse_args()

    hdr = parse_header(a.vcf)
    meta = detect_annotations(hdr)
    if a.gnomad_tag in hdr["info"]:
        meta.setdefault("gnomAD_fields", [])
        if a.gnomad_tag not in meta["gnomAD_fields"]:
            meta["gnomAD_fields"].insert(0, a.gnomad_tag)
    ac = AnnotationCounter(meta)
    k = Counter()
    with open_vcf(a.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            ref, alt = f[3], f[4].split(",")[0]
            cl, _ = classify_allele(ref, alt)
            vt = "SNV" if cl == "SNV" else ("INDEL" if cl in ("insertion", "deletion", "complex_indel") else "OTHER")
            info = parse_info(f[7])
            rs = f[2].startswith("rs")
            gn = to_float(str(info.get(a.gnomad_tag, ".")).split(",")[0]) is not None
            db = bool(a.dbsnp_tag) and a.dbsnp_tag in info
            k["total"] += 1
            k["rsID"] += rs or db
            k["gnomAD"] += gn
            k["known_any"] += rs or gn or db
            k[f"{vt}_total"] += 1
            k[f"{vt}_known"] += rs or gn or db
            ac.add(info)
    res = {"known": dict(k), "annotation": ac.result()}
    write_json(res, a.out)
    print(f"[annotation_summary] {k['total']:,} records, known={k['known_any']:,}", file=sys.stderr)


if __name__ == "__main__":
    main()
