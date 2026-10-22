#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
A=(--sizemin "$SVMIN" --sizefilt "$SVMIN" --sizemax 3000000 -P 0 -O 0.5 -r 1000 --pick multi)
rm -rf "$TV/tv_manta" "$TV/tv_inv"
with "$TRUVARI_BIN" truvari bench -b "$MD/manta.pass.vcf.gz" -c "$CR/svs.depthvalidated.vcf.gz" -o "$TV/tv_manta" -t "${A[@]}"
with "$TRUVARI_BIN" truvari bench -b "$MD/manta.pass.inv.vcf.gz" -c "$C/cue.inv.vcf.gz" -o "$TV/tv_inv" "${A[@]}"

bcftools concat -a "$TV"/tv_{manta,inv}/tp-comp.vcf.gz | bcftools sort -Oz -o "$OUT/cue.manta_confirmed.sorted.vcf.gz" -
bcftools index -ft "$OUT/cue.manta_confirmed.sorted.vcf.gz"

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL\n' "$TV/tv_manta/fp.vcf.gz" \
    | awk -v OFS='\t' 'NR == FNR {r[$1 ":" $2 + 1] = $7; next} {print $0, $3 - $2, r[$1 ":" $2]}' "$CR/dd.ratio.tsv" - \
    | sort -k6,6n > "$OUT/igv_review.tsv"
