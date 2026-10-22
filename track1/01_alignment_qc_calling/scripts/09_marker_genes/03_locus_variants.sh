#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
M=$RESULTS/marker_genes
FLANK=50000
GTF=$(echo "$PROJECT_ROOT"/references/annotation/ensembl_116/*.gtf.gz)
DV_TABLE=$RESULTS/annotation/deepvariant/05_tables/high_confidence.annotated.tsv
PV_TABLE=$TRACK1/05_vcf_qc_annotation_ranking/results/hc_out/05_tables/high_confidence.annotated.tsv
mkdir -p "$M"

zcat "$GTF" | awk -F'\t' -v OFS='\t' -v genes="$MVA_GENES" -v f="$FLANK" 'BEGIN {n = split(genes, a, "|"); for (i = 1; i <= n; i++) g[a[i]]}
    $3 == "gene" && match($9, /gene_name "[^"]+"/) && (nm = substr($9, RSTART + 11, RLENGTH - 12)) in g {print nm, "chr" $1, ($4 > f ? $4 - f : 1), $5 + f, $4, $5, $7}' |
    sort -k2,2V -k3,3n > "$M/locus_windows.tsv"

for s in "deepvariant:$DV_TABLE" "provided:$PV_TABLE"; do
    awk -F'\t' -v OFS='\t' -v set="${s%%:*}" '
        FILENAME == ARGV[1] {k++; wg[k] = $1; wc[k] = $2; ws[k] = $3; we[k] = $4; gs[k] = $5; ge[k] = $6; st[k] = $7; next}
        FNR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
        {ch = $(c["CHROM"]); p = $(c["POS"]) + 0
         for (i = 1; i <= k; i++) if (ch == wc[i] && p >= ws[i] && p <= we[i]) {
             r = p < gs[i] ? (st[i] == "+" ? "upstream" : "downstream") : p > ge[i] ? (st[i] == "+" ? "downstream" : "upstream") : "gene_body"
             print wg[i], r, ch, p, $(c["REF"]), $(c["ALT"]), $(c["VARIANT_CLASS"]), $(c["Consequence"]), $(c["SYMBOL"]), $(c["IMPACT"]),
                 $(c["HGVSc"]), $(c["HGVSp"]), $(c["MAX_AF"]), $(c["ClinVar_CLNSIG"]), $(c["am_class"]), set ":" $(c["GT"]) ":" $(c["VAF"]), ($(c["SYMBOL"]) == wg[i] ? 0 : 1)}}' \
        "$M/locus_windows.tsv" "${s#*:}"
done | sort -t$'\t' -k1,1 -k4,4n -k5,5 -k6,6 -k17,17n |
    awk -F'\t' -v OFS='\t' '{k = $1 FS $3 FS $4 FS $5 FS $6; split($16, x, ":")}
        k != p {if (p != "") print r, s; r = $1; for (i = 2; i <= 15; i++) r = r OFS ($i == "" ? "-" : $i); s = $16; delete seen; seen[x[1]]; p = k; next}
        !(x[1] in seen) {s = s "," $16; seen[x[1]]}
        END {if (p != "") print r, s}' |
    { printf 'window_gene\tregion\tchrom\tpos\tref\talt\tvariant_class\tconsequence\tvep_symbol\timpact\thgvsc\thgvsp\tmax_AF\tclinvar\talphamissense\tcall_sets\n'; cat; } > "$M/locus_variants_50kb.tsv"

awk -F'\t' 'NR == 1 || ($9 == $1 && ($10 ~ /^(HIGH|MODERATE)$/ || $8 ~ /splice/))' "$M/locus_variants_50kb.tsv" > "$M/locus_coding_variants.tsv"
