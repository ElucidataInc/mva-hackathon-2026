#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
bcftools view -f PASS -Oz -o "$MD/manta.pass.vcf.gz" "$MD/manta_run/results/variants/diploidSV.vcf.gz"
with "$MANTA_BIN" python2 "$MANTA_BIN/convertInversion.py" "$MANTA_BIN/samtools" "$REF" "$MD/manta.pass.vcf.gz" \
    | bcftools view -Oz -o "$MD/manta.pass.inv.vcf.gz" -
for v in "$MD"/manta.pass{,.inv}.vcf.gz; do bcftools index -ft "$v"; done
