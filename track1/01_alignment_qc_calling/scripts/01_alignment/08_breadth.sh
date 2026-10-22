#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
O=$QC/coverage/$SAMPLE
run phenotype python - "$REF" "$O.thresholds.bed.gz" <<'PY' > "$O.breadth.tsv"
import gzip, re, sys
from collections import Counter, defaultdict

ref, thr = sys.argv[1:]
main = re.compile(r"chr([0-9]+|X|Y)$")

acgt, chrom = Counter(), None
with open(ref) as f:
    for line in f:
        if line[0] == ">":
            chrom = line[1:].split()[0]
        elif main.match(chrom):
            s = line.rstrip()
            acgt[chrom] += len(s) - s.upper().count("N")

with gzip.open(thr, "rt") as f:
    levels = next(f).split()[4:]
    cov = defaultdict(lambda: [0] * len(levels))
    for line in f:
        c, *rest = line.split()
        if main.match(c):
            cov[c] = [a + int(b) for a, b in zip(cov[c], rest[3:])]

auto = [f"chr{i}" for i in range(1, 23)]
print("contig", "non_N_bp", *(f"pct_ge{x.lower()}" for x in levels), sep="\t")
for name, group in [*((c, [c]) for c in auto + ["chrX", "chrY"]), ("autosomes", auto)]:
    n = sum(acgt[c] for c in group)
    k = [sum(cov[c][i] for c in group) for i in range(len(levels))]
    print(name, n, *(f"{100 * x / n:.2f}" for x in k), sep="\t")
PY
