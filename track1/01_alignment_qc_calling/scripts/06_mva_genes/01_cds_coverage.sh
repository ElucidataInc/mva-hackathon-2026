#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
G=$RESULTS/mva_genes
GTF=$(echo "$PROJECT_ROOT"/references/annotation/ensembl_116/*.gtf.gz)
mkdir -p "$G"

zcat "$GTF" | awk -F'\t' -v OFS='\t' -v g="$MVA_GENES" '$3 == "CDS" && /tag "Ensembl_canonical"/ && $9 ~ "gene_name \"(" g ")\";" {
    match($9, /gene_name "[^"]+"/); gene = substr($9, RSTART + 11, RLENGTH - 12)
    match($9, /exon_number "[^"]+"/); print "chr" $1, $4 - 1, $5, gene "_exon" substr($9, RSTART + 13, RLENGTH - 14)}' |
    sort -k1,1 -k2,2n > "$G/mva_cds.bed"

run align mosdepth -t 4 -n -Q 20 -F 3844 -T 10,20,30 -b "$G/mva_cds.bed" "$G/$SAMPLE.cds" "$FINAL_BAM"

paste <(zcat "$G/$SAMPLE.cds.regions.bed.gz") <(zcat "$G/$SAMPLE.cds.thresholds.bed.gz" | grep -v '^#') | awk -v OFS='\t' '
    $2 != $7 || $4 != $9 {print "region order differs: " $4 " " $9 > "/dev/stderr"; exit 1}
    {split($4, a, "_"); g = a[1]; L = $3 - $2; n[g] += L; d[g] += $5 * L; t10[g] += $10; t20[g] += $11; t30[g] += $12
     if ($11 < L) low = low sprintf("%s\t%s:%d-%d\t%d bp\t%.1f%% >=20x\n", $4, $1, $2 + 1, $3, L, 100 * $11 / L)}
    END {print "gene", "cds_bp", "mean_depth", "pct_10x", "pct_20x", "pct_30x"
         for (g in n) printf "%s\t%d\t%.1f\t%.2f\t%.2f\t%.2f\n", g, n[g], d[g] / n[g], 100 * t10[g] / n[g], 100 * t20[g] / n[g], 100 * t30[g] / n[g]
         printf "exons below 20x:\n%s", low ? low : "none\n"}' > "$G/cds_coverage_summary.tsv"
