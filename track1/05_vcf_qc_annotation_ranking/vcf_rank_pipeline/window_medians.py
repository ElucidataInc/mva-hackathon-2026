#!/usr/bin/env python3
"""Per-chromosome median depth of 1-Mb windows (stdin: chrom start end mean_depth).
Windows below 30% of the autosomal median window depth (N gaps, centromeres, acrocentric short arms)
are ignored, so chromosomes with large gaps are not reported as 'low'.
Output: chrom <tab> median_window_depth <tab> windows_used; first line on stderr: autosomal median window depth."""
import re
import sys
from collections import defaultdict
from statistics import median

win = defaultdict(list)
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4 and re.fullmatch(r"(chr)?([0-9]+|X|Y)", f[0]):
        try:
            win[f[0]].append(float(f[3]))
        except ValueError:
            pass
auto = [d for c, v in win.items() if re.fullmatch(r"(chr)?[0-9]+", c) for d in v if d > 0]
if not auto:
    sys.exit("no autosomal windows with depth")
gmed = median(auto)
print(f"{gmed:.2f}", file=sys.stderr)
key = lambda c: (0, int(c.replace("chr", ""))) if c.replace("chr", "").isdigit() else (1, c)
for c in sorted(win, key=key):
    use = [d for d in win[c] if d >= 0.3 * gmed]
    # chrY in a female (or tiny chrY coverage): report the raw median so it is visibly ~0
    print(f"{c}\t{median(use) if use else median(win[c]):.2f}\t{len(use)}")
