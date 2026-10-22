#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen --no-install -e qc)"
M=$RESULTS/marker_genes
GTF=$(echo "$PROJECT_ROOT"/references/annotation/ensembl_116/*.gtf.gz)
CLINVAR=$(echo "$PROJECT_ROOT"/references/clinvar/*/*.vcf.gz)
MAP=$RESULTS/deepvariant/compare/chr_map.txt
DV=$RESULTS/deepvariant/$SAMPLE.deepvariant.vcf.gz
EX=$RESULTS/exomiser
mkdir -p "$M"

zcat "$GTF" | awk -F'\t' -v OFS='\t' 'NR == FNR {if (FNR > 1) g[$1]; next}
    $3 == "gene" && match($9, /gene_name "[^"]+"/) && (n = substr($9, RSTART + 11, RLENGTH - 12)) in g {print "chr" $1, $4, $5, n}' "$SCRIPTS/marker_genes.tsv" - |
    sort -k1,1V -k2,2n > "$M/genes.tsv"
[[ $(wc -l < "$M/genes.tsv") -eq $(($(wc -l < "$SCRIPTS/marker_genes.tsv") - 1)) ]] || { echo "not every marker gene found once in the GTF" >&2; exit 1; }
{ cut -f1-3 "$M/genes.tsv"; awk -v OFS='\t' 'NR > 1 {print $1, $2, $2}' "$SCRIPTS/literature_alleles.tsv"; } | sort -k1,1V -k2,2n > "$M/regions.tsv"
sed 's/^chr//' "$M/regions.tsv" > "$M/regions_nochr.tsv"

norm() { bcftools norm -m -any -f "$REF" -Ov | bcftools view -e 'ALT="*"' -Oz -o "$1"; bcftools index -f -t "$1"; }
bcftools view -R "$M/regions_nochr.tsv" -i 'INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic"' -Ov "$CLINVAR" | bcftools annotate --rename-chrs "$MAP" -Ov | norm "$M/clinvar_plp.vcf.gz"
bcftools view -R "$M/regions.tsv" -Ov "$DV" | norm "$M/deepvariant.vcf.gz"
bcftools view -R "$M/regions_nochr.tsv" -Ov "$VCF" | bcftools annotate --rename-chrs "$MAP" -Ov | norm "$M/provided.vcf.gz"

q() { bcftools query -f '%CHROM:%POS:%REF:%ALT\t%FILTER\t[%GT]\n' "$1"; }
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENEINFO\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/MC\n' "$M/clinvar_plp.vcf.gz" |
    awk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] {g[$4]; next} FILENAME == ARGV[2] {d[$1] = $2 ":" $3; next} FILENAME == ARGV[3] {p[$1] = $2 ":" $3; next}
        {n = split($5, a, "|"); hit = ""; for (i = 1; i <= n; i++) {split(a[i], b, ":"); if (b[1] in g) hit = b[1]}
         if (hit == "") next
         k = $1 ":" $2 ":" $3 ":" $4
         print hit, $1, $2, $3, $4, $6, $7, $8, (k in d ? d[k] : "absent"), (k in p ? p[k] : "absent")}' \
        "$M/genes.tsv" <(q "$M/deepvariant.vcf.gz") <(q "$M/provided.vcf.gz") - > "$M/clinvar_plp_status.tsv"

rank() { awk -F'\t' -v g="$1" 'NR > 1 && $3 == g {r = r (r ? "," : "") $1 $5} END {print r ? r : "not ranked"}' "$2"; }
{
    printf 'gene\tphenotype\tclinvar_plp\tdeepvariant_pass\tprovided_pass\texomiser_proband_terms\texomiser_all_terms\n'
    awk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] {n[$1]++; d[$1] += $9 ~ /^PASS:/; p[$1] += $10 ~ /^PASS:/; next}
        FNR > 1 {print $1, $2, n[$1] + 0, d[$1] + 0, p[$1] + 0}' "$M/clinvar_plp_status.tsv" "$SCRIPTS/marker_genes.tsv" |
        while IFS=$'\t' read -r g ph n d p; do
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$g" "$ph" "$n" "$d" "$p" \
                "$(rank "$g" "$EX/proband_terms/$SAMPLE.proband_terms.genes.tsv")" "$(rank "$g" "$EX/all_terms/$SAMPLE.all_terms.genes.tsv")"
        done
    printf '\nClinVar P/LP variants found in either call set\n'
    awk -F'\t' '$9 != "absent" || $10 != "absent"' "$M/clinvar_plp_status.tsv"
    printf '\nLiterature alleles not in ClinVar\n'
    tail -n +2 "$SCRIPTS/literature_alleles.tsv" | while IFS=$'\t' read -r c p r a g s; do
        ref=$(samtools faidx "$REF" "$c:$p-$p" | tail -n +2)
        dv=$(bcftools query -r "$c:$p" -f '%REF>%ALT:%FILTER:[%GT]' "$M/deepvariant.vcf.gz")
        pv=$(bcftools query -r "$c:$p" -f '%REF>%ALT:%FILTER:[%GT]' "$M/provided.vcf.gz")
        reads=$(bcftools mpileup -f "$REF" -r "$c:$p-$p" -q 20 -Q 20 -a FORMAT/AD "$FINAL_BAM" 2>/dev/null | bcftools query -f '%REF>%ALT:[%AD]')
        printf '%s\t%s:%s %s>%s\treference=%s\tdeepvariant=%s\tprovided=%s\treads_MAPQ>=20=%s\t%s\n' "$g" "$c" "$p" "$r" "$a" "$ref" "${dv:-none}" "${pv:-none}" "$reads" "$s"
    done
} > "$M/known_variants_summary.tsv"
