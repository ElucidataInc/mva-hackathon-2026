#!/usr/bin/env python3
"""
Merge the per-gene tables of two phased, VEP-annotated VCFs of the same sample.

Rows are matched on CHROM, POS, REF, ALT and gene (both VCFs must be normalised the same way:
multi-allelic sites split, indels left-aligned - vcf_hc_vep.sh does this).

  --combine union      every variant found in either file (default)
  --combine consensus  only variants found in both files

Added columns: SOURCE (both / <label1> / <label2>) and GT_NOTE (genotype disagreement between files).
Different callers can write the same indel differently (e.g. shifted in a repeat). In union mode a
<label1>-only and a <label2>-only indel/MNV in the same gene within --window bp are treated as one
variant (SOURCE "both", note "representation differs"), so they cannot form a false compound het
(union) and are not lost (consensus).
For a variant in both files the values of file 1 are used; empty fields are filled from file 2.
Phase sets are prefixed with the file label (e.g. "A:123456"): a PS from one file is never compared
with a PS from the other file. If file 1 has no phase for a shared variant but file 2 does, the
genotype and PS of file 2 are used.
"""
import argparse
import csv
import sys


def norm_chrom(c, style_chr):
    b = c[3:] if c.startswith("chr") else c
    b = "M" if b in ("MT", "M") else b
    return ("chr" + b) if style_chr else ("MT" if b == "M" else b)


def read(path, label, style_chr):
    rows, order = {}, []
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        cols = rdr.fieldnames or []
        for r in rdr:
            r["CHROM"] = norm_chrom(r["CHROM"], style_chr)
            ps = r.get("PS", "")
            if ps not in ("", ".", None):
                r["PS"] = f"{label}:{ps}"
            g = r.get("SYMBOL") or r.get("Gene") or "."
            key = (r["CHROM"], r["POS"], r["REF"], r["ALT"], g)
            if key not in rows:
                order.append(key)
            rows[key] = r
    return cols, rows, order


def phased(gt):
    return "|" in (gt or "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t1", required=True)
    ap.add_argument("--t2", required=True)
    ap.add_argument("--label1", default="A")
    ap.add_argument("--label2", default="B")
    ap.add_argument("--combine", choices=["union", "consensus"], default="union")
    ap.add_argument("--style-chr", type=int, default=1)
    ap.add_argument("--window", type=int, default=10, help="bp window to match indel representations [10]")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    c1, r1, o1 = read(a.t1, a.label1, a.style_chr)
    c2, r2, o2 = read(a.t2, a.label2, a.style_chr)
    cols = list(c1) + [c for c in c2 if c not in c1] + ["SOURCE", "GT_NOTE"]
    keys = [k for k in o1] + [k for k in o2 if k not in r1]
    n = {"both": 0, a.label1: 0, a.label2: 0, "gt_discordant": 0}
    out = []
    for k in keys:
        x, y = r1.get(k), r2.get(k)
        if x and y:
            src = "both"
            row = dict(y)
            row.update({c: v for c, v in x.items() if v not in ("", ".", None)})
            note = ""
            g1, g2 = x.get("GT", ""), y.get("GT", "")
            if g1.replace("|", "/") != g2.replace("|", "/") and \
               sorted(g1.replace("|", "/").split("/")) != sorted(g2.replace("|", "/").split("/")):
                note = f"GT differs: {a.label1}={g1} {a.label2}={g2}"
                n["gt_discordant"] += 1
            if not phased(g1) and phased(g2):            # take phase from file 2
                row["GT"], row["PS"] = g2, y.get("PS", "")
            row["GT_NOTE"] = note or "-"
        else:
            src = a.label1 if x else a.label2
            row = dict(x or y)
            row["GT_NOTE"] = "-"
        row["SOURCE"] = src
        n[src] += 1
        out.append(row)

    # collapse <label1>-only / <label2>-only indel or MNV pairs that are probably the same variant
    def complex_(r):
        return len(r["REF"]) != len(r["ALT"]) or len(r["REF"]) > 1
    only1 = [r for r in out if r["SOURCE"] == a.label1]
    only2 = [r for r in out if r["SOURCE"] == a.label2]
    by_gene = {}
    for r in only2:
        by_gene.setdefault((r["CHROM"], r.get("SYMBOL") or r.get("Gene") or "."), []).append(r)
    drop = set()
    for r in only1:
        for q in by_gene.get((r["CHROM"], r.get("SYMBOL") or r.get("Gene") or "."), []):
            if id(q) in drop or not (complex_(r) or complex_(q)):
                continue
            if abs(int(r["POS"]) - int(q["POS"])) <= a.window:
                drop.add(id(q))
                r["SOURCE"] = "both"
                r["GT_NOTE"] = (f"representation differs: {a.label2} has {q['POS']} {q['REF']}>{q['ALT']} {q.get('GT', '')}"
                                + ("" if r["GT_NOTE"] in ("", "-") else "; " + r["GT_NOTE"]))
                n["both"] += 1
                n[a.label1] -= 1
                n[a.label2] -= 1
                n["representation"] = n.get("representation", 0) + 1
                break
    out = [r for r in out if id(r) not in drop]
    if a.combine == "consensus":                     # after collapsing, so re-represented indels count as shared
        out = [r for r in out if r["SOURCE"] == "both"]
        n[a.label1] = n[a.label2] = 0

    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r.get(c, "") if r.get(c, "") is not None else "") for c in cols) + "\n")
    print(f"[merge] {a.combine}: {len(out):,} variant-gene rows | in both: {n['both']:,} | "
          f"{a.label1} only: {n[a.label1]:,} | {a.label2} only: {n[a.label2]:,} | "
          f"genotype differs: {n['gt_discordant']:,} | indel representation differs (merged): {n.get('representation', 0):,}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
