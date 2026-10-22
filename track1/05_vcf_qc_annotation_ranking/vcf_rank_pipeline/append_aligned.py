#!/usr/bin/env python3
"""Append the rows of table B to table A, matching columns by name (columns A lacks in B become ".").
Prints: <rows appended> <comma-separated A columns missing from B, or ->"""
import sys

tab, extra = sys.argv[1], sys.argv[2]
with open(tab) as fh:
    cols = fh.readline().rstrip("\n").split("\t")
with open(extra) as fh:
    ecols = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(ecols)}
    n = 0
    with open(tab, "a") as out:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.write("\t".join(f[idx[c]] if c in idx and idx[c] < len(f) else "." for c in cols) + "\n")
            n += 1
print(n, ",".join(c for c in cols if c not in idx) or "-")
