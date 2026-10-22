#!/usr/bin/env bash
set -euo pipefail
P=/home/ubuntu/mva-hackathon-2026
run() { local e=$1; shift; "$P/tools/pixi/0.81.0/bin/pixi" run --manifest-path "$P/tools/pixi/pixi.toml" --frozen --no-install -e "$e" "$@"; }

NCBI=$P/references/genome/ncbi_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
GIAB=$P/references/genome/giab_candidate/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz
VCF=$(echo "$P"/data/variants/*.vcf.gz)
D=$P/references/genome/hs38d1_maskedGRC_v2
OUT=$D/GRCh38_no_alt_plus_hs38d1_maskedGRC_exclusions_v2.fna
mkdir -p "$D"

# NCBI contigs in NCBI order, with GIAB masked sequence where GIAB has the contig.
# Stops if lengths differ or if any base differs other than by GIAB masking to N.
run phenotype python - "$NCBI" "$GIAB" "$OUT" > "$D/masked_bases.tsv" <<'PY'
import gzip, sys
import numpy as np

ncbi, giab, out = sys.argv[1:4]

def fasta(f):
    name, seq = None, []
    for line in f:
        if line.startswith(b">"):
            if name:
                yield name, b"".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.rstrip())
    if name:
        yield name, b"".join(seq)

with gzip.open(giab, "rb") as f:
    masked = dict(fasta(f))

total = 0
with open(ncbi, "rb") as f, open(out, "wb") as o:
    for name, seq in fasta(f):
        new = masked.pop(name, seq)
        if len(new) != len(seq):
            sys.exit(f"{name.decode()}: length differs")
        a = np.frombuffer(seq.upper(), np.uint8)
        b = np.frombuffer(new.upper(), np.uint8)
        diff = a != b
        if (b[diff] != ord("N")).any():
            sys.exit(f"{name.decode()}: differs outside masked bases")
        if diff.any():
            print(name.decode(), int(diff.sum()), sep="\t")
            total += int(diff.sum())
        o.write(b">" + name + b"\n")
        for i in range(0, len(new), 60):
            o.write(new[i:i + 60] + b"\n")
if masked:
    sys.exit(f"GIAB contigs not in NCBI: {[k.decode() for k in masked]}")
print("total", total, sep="\t")
PY

run qc samtools faidx "$OUT"
cmp <(cut -f1,2 "$NCBI.fai") <(cut -f1,2 "$OUT.fai")
cmp <(cut -f1,2 "$OUT.fai" | sed 's/^chr//') \
    <(run qc bcftools view -h "$VCF" | sed -n 's/^##contig=<ID=\([^,]*\),length=\([0-9]*\)[,>].*/\1\t\2/p')
run qc samtools dict -o "${OUT%.fna}.dict" "$OUT"
run align bwa index "$OUT" 2> "$D/bwa_index.log"
md5sum "$OUT" > "$OUT.md5"
