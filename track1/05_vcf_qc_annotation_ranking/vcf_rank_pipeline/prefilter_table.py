#!/usr/bin/env python3
"""
Stream a per-gene variant table (stdin -> stdout) and keep only rows that can become a ranking
candidate: a gene symbol, an ALT allele in the genotype, rare (population AF below --max-af) or
ClinVar P/LP, and a qualifying evidence class (same rules as rank_variants.py).
Whole-genome tables (millions of rows) are reduced to a few thousand rows without holding them in memory.
"""
import argparse
import csv
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rank_variants import alt_copies, classify, get, max_af  # noqa: E402

csv.field_size_limit(sys.maxsize)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-af", type=float, default=0.01)
    ap.add_argument("--cadd-min", type=float, default=20)
    ap.add_argument("--log", help="write counts here")
    ap.add_argument("--keep-sites", help="TSV with chrom pos ref alt columns (e.g. AlphaGenome table): always kept if ALT is carried")
    a = ap.parse_args()
    rdr = csv.DictReader(sys.stdin, delimiter="\t")
    out = csv.DictWriter(sys.stdout, fieldnames=rdr.fieldnames, delimiter="\t", lineterminator="\n",
                         extrasaction="ignore")
    out.writeheader()
    keep = set()
    if a.keep_sites:
        with open(a.keep_sites) as fh:
            for k in csv.DictReader(fh, delimiter="\t"):
                try:
                    keep.add((k["chrom"][3:] if k["chrom"].startswith("chr") else k["chrom"],
                              str(int(float(k["pos"]))), k["ref"], k["alt"]))
                except (KeyError, ValueError):
                    pass
    n = Counter()
    for r in rdr:
        n["rows in"] += 1
        c = r.get("CHROM", "")
        if keep and (c[3:] if c.startswith("chr") else c, r.get("POS"), r.get("REF"), r.get("ALT")) in keep \
                and alt_copies(get(r, "GT") or "0/1") > 0:
            out.writerow(r)
            n["kept as listed site (AlphaGenome)"] += 1
            continue
        if not (get(r, "SYMBOL") or get(r, "Gene")):
            n["no gene"] += 1
            continue
        if alt_copies(get(r, "GT") or "0/1") == 0:
            n["no ALT allele in genotype"] += 1
            continue
        af = max_af(r)
        if af is not None and af >= a.max_af and "pathogenic" not in get(r, "ClinVar_CLNSIG").lower():
            n["common"] += 1
            continue
        cls, _ = classify(r, a)
        if cls is None:
            n["no qualifying evidence"] += 1
            continue
        out.writerow(r)
        n["rows kept"] += 1
    msg = "[prefilter] " + " | ".join(f"{k}: {v:,}" for k, v in n.items())
    if a.log:
        with open(a.log, "w") as fh:
            fh.write(msg + "\n")
    print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()
