#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
M=$RESULTS/marker_genes
MAX_AF=0.01
DV_TABLE=$RESULTS/annotation/deepvariant/05_tables/high_confidence.annotated.tsv
PV_TABLE=$TRACK1/05_vcf_qc_annotation_ranking/results/hc_out/05_tables/high_confidence.annotated.tsv
mkdir -p "$M"

for s in "deepvariant:$DV_TABLE" "provided:$PV_TABLE"; do
    awk -F'\t' -v OFS='\t' -v genes="$MVA_GENES" -v maxaf="$MAX_AF" -v set="${s%%:*}" '
        NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; n = split(genes, g, "|"); for (i = 1; i <= n; i++) G[g[i]]; next}
        ($(c["SYMBOL"]) in G) && $(c["CANONICAL"]) == "YES" && ($(c["IMPACT"]) ~ /^(HIGH|MODERATE)$/ || $(c["Consequence"]) ~ /splice/) {
            af = $(c["MAX_AF"]); if (af == "" || af == ".") af = "-"; if (af != "-" && af + 0 >= maxaf) next
            print $(c["SYMBOL"]), $(c["CHROM"]), $(c["POS"]), $(c["REF"]), $(c["ALT"]), $(c["Consequence"]), $(c["HGVSc"]), $(c["HGVSp"]),
                $(c["gnomADe_AF"]), $(c["gnomADg_AF"]), af, $(c["ClinVar_CLNSIG"]), $(c["am_class"]), set ":" $(c["GT"]) ":" $(c["VAF"])}' "${s#*:}"
done | sort -t$'\t' -k1,1 -k3,3n -k4,4 -k5,5 |
    awk -F'\t' -v OFS='\t' '{k = $1 FS $2 FS $3 FS $4 FS $5} k != p {if (p != "") print r, s; r = $1; for (i = 2; i <= 13; i++) r = r OFS ($i == "" ? "-" : $i); s = $14; p = k; next} {s = s "," $14} END {if (p != "") print r, s}' |
    { printf 'gene\tchrom\tpos\tref\talt\tconsequence\thgvsc\thgvsp\tgnomADe_AF\tgnomADg_AF\tmax_AF\tclinvar\talphamissense\tcall_sets\n'; cat; } > "$M/rare_coding_variants.tsv"

awk -F'\t' 'NR > 1 {n[$1]++} END {for (g in n) print g, n[g]}' "$M/rare_coding_variants.tsv" | sort
