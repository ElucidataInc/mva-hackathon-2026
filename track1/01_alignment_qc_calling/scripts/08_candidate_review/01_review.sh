#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen --no-install -e qc)"
C=$RESULTS/candidate_review
EXV=$RESULTS/exomiser/proband_terms/$SAMPLE.proband_terms.variants.tsv
DV=$RESULTS/deepvariant/$SAMPLE.deepvariant.vcf.gz
PH=$RESULTS/phasing/$SAMPLE.phased.vcf.gz
mkdir -p "$C"

prefix() { bcftools index -s "$1" | awk 'NR == 1 {print ($1 ~ /^chr/) ? "chr" : ""}'; }
PV=$(prefix "$VCF"); PD=$(prefix "$DV"); PP=$(prefix "$PH")
PB=$(awk 'NR == 1 {print ($1 ~ /^chr/) ? "chr" : ""}' "$REF.fai")

awk -F'\t' -v OFS='\t' 'NR == 2 {g = $3} NR > 1 && $3 == g && $11 == 1 {sub(/^chr/, "", $15); print $3, $15, $16, $18, $19}' "$EXV" |
    sort -u | sort -k2,2V -k3,3n > "$C/candidates.tsv"
[[ -s $C/candidates.tsv ]] || { echo "no contributing variants for the top-ranked gene" >&2; exit 1; }

{
    while IFS=$'\t' read -r g c p r a; do
        bcftools query -r "$PV$c:$p" -f "provided\t$g\t$c:$p\t%REF>%ALT\t%FILTER\t[%GT\t%AD\t%DP\t%GQ]\n" "$VCF"
        bcftools query -r "$PD$c:$p" -f "deepvariant\t$g\t$c:$p\t%REF>%ALT\t%FILTER\t[%GT\t%AD\t%DP\t%GQ\t%VAF]\n" "$DV"
        bcftools query -r "$PP$c:$p" -f "phased\t$g\t$c:$p\t%REF>%ALT\t[%GT\t%PS]\n" "$PH"
        for q in 0 20; do
            bcftools mpileup -f "$REF" -r "$PB$c:$p-$p" -q $q -Q 20 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP "$FINAL_BAM" 2>/dev/null |
                bcftools query -f "reads_MAPQ>=$q\t$g\t$c:$p\t%REF>%ALT\t[%DP\t%AD\t%ADF\t%ADR]\n"
        done
    done < "$C/candidates.tsv"
    awk -F'\t' 'NR == 1 {c = $2; s = $3} $2 != c {exit 1} {e = $3} END {print c "\t" s "\t" e}' "$C/candidates.tsv" |
        while read -r c s e; do
            bcftools query -r "$PP$c:$s-$e" -i 'GT="het"' -f '[%PS]\n' "$PH" | sort | uniq -c | awk -v OFS='\t' '{print "het_sites_by_phase_block", $2, $1}'
        done
} > "$C/review.tsv"
