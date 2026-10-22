#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen -e qc)"

SRC=$IN
if ! bcftools index -s "$IN" | awk '$1 ~ /^chr([0-9]+|X|Y)$/ {n++} END {exit !n}'; then
    awk -v OFS='\t' '{c = $1; sub(/^chr/, "", c); print c, $1}' "$REF.fai" > "$PH/chr_names.tsv"
    SRC=$PH/$SAMPLE.chr.vcf.gz
    bcftools annotate --rename-chrs "$PH/chr_names.tsv" -Oz -o "$SRC" "$IN"
    bcftools index -t "$SRC"
fi
for c in $CHROMS; do bcftools view -r "$c" -Oz -o "$PH/chrom/$c.vcf.gz" "$SRC"; done
