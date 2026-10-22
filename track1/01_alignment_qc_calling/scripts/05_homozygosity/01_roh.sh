#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen --no-install -e qc)"
R=$RESULTS/homozygosity
GTF=$(echo "$PROJECT_ROOT"/references/annotation/ensembl_116/*.gtf.gz)
AUTO=$(seq -s, -f 'chr%g' 1 22)
mkdir -p "$R"

bcftools view -f PASS -m2 -M2 -v snps -r "$AUTO" -Ou "$RESULTS/deepvariant/$SAMPLE.deepvariant.vcf.gz" |
    bcftools roh -G30 --AF-dflt 0.4 -O r -o "$R/$SAMPLE.roh.txt"

zcat -f "$GTF" | awk -F'\t' -v OFS='\t' -v g="$MVA_GENES" '$3 == "gene" && $9 ~ "gene_name \"(" g ")\";" {
    match($9, /gene_name "[^"]+"/); print "chr" $1, $4 - 1, $5, substr($9, RSTART + 11, RLENGTH - 12)}' > "$R/mva_genes.bed"

{
    awk -v OFS='\t' -v G="$(awk '$1 ~ /^chr([0-9]+)$/ {s += $2} END {print s}' "$REF.fai")" '$1 == "RG" {
            n++; if ($6 > m) m = $6
            if ($6 >= 1e6) {n1++; l1 += $6}
            if ($6 >= 5e6) {n5++; l5 += $6}}
        END {print "segments", n; print "segments_1Mb", n1; print "total_1Mb_Mb", l1 / 1e6
             print "segments_5Mb", n5; print "total_5Mb_Mb", l5 / 1e6; print "longest_Mb", m / 1e6
             printf "F_ROH_1Mb\t%.4f\n", l1 / G}' "$R/$SAMPLE.roh.txt"
    awk -v OFS='\t' 'NR == FNR {if ($1 == "RG") {k++; c[k] = $3; s[k] = $4; e[k] = $5; l[k] = $6}; next}
        {hit = "none"; for (i = 1; i <= k; i++) if (c[i] == $1 && s[i] <= $3 && e[i] > $2) hit = sprintf("%s:%d-%d (%.2f Mb)", c[i], s[i], e[i], l[i] / 1e6)
         print $4, $1 ":" $2 + 1 "-" $3, hit}' "$R/$SAMPLE.roh.txt" "$R/mva_genes.bed"
} > "$R/roh_summary.tsv"
