#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen -e qc)"
O=$QC/identity
mkdir -p "$O"

bcftools view -m2 -M2 -v snps -f PASS -i 'FMT/DP >= 20 && FMT/GQ >= 30' "$VCF" |
    bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t[%GT]\n' |
    awk -v OFS='\t' '$1 ~ /^[0-9]+$/ && ++n % 50 == 0 {print "chr" $1, $2, $3, $4}' > "$O/sites.tsv"
cut -f1,2 "$O/sites.tsv" | bgzip > "$O/regions.tsv.gz"
cut -f1-3 "$O/sites.tsv" | bgzip > "$O/alleles.tsv.gz"
tabix -s1 -b2 -e2 "$O/regions.tsv.gz"
tabix -s1 -b2 -e2 "$O/alleles.tsv.gz"

bcftools mpileup -f "$REF" -R "$O/regions.tsv.gz" -a AD,DP -q 20 -Q 20 -Ou "$FINAL_BAM" |
    bcftools call -m -C alleles -T "$O/alleles.tsv.gz" -Oz -o "$O/bam_calls.vcf.gz"

bcftools query -f '%CHROM\t%POS\t[%GT\t%AD]\n' "$O/bam_calls.vcf.gz" |
awk -F'\t' '
    function g(x,  a) { gsub(/\|/, "/", x); split(x, a, "/"); return a[1] <= a[2] ? a[1] "/" a[2] : a[2] "/" a[1] }
    NR == FNR { truth[$1 ":" $2] = g($4); next }
    { k = $1 ":" $2; split($4, ad, ","); d = ad[1] + ad[2] }
    !(k in truth) || $3 ~ /\./ || d < 20 { next }
    { n++; same += g($3) == truth[k] }
    truth[k] == "1/1" { h++; ref += ad[1] / d }
    truth[k] == "0/1" { e++; alt += ad[2] / d }
    END { printf "sites\t%d\nconcordance\t%.4f\nhom_alt_sites\t%d\nref_fraction_at_hom_alt\t%.4f\nhet_sites\t%d\nalt_fraction_at_het\t%.4f\n",
                 n, same / n, h, ref / h, e, alt / e }
' "$O/sites.tsv" - > "$O/$SAMPLE.identity.tsv"
