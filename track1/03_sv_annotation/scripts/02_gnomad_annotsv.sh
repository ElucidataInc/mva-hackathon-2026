#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
(cd "$OUT" && svafotate annotate -v "$OUT/combined.mane.vcf.gz" -o "$OUT/combined.mane.gnomad.vcf" \
    -b "$SVREF/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz" -s gnomAD -f 0.5 -a best --cpu "$THREADS")
bcftools view "$OUT/combined.mane.gnomad.vcf" \
    | awk -v OFS='\t' '/^##/ {print; next} /^#CHROM/ {print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; print $0, "FORMAT", "SAMPLE"; next} {print $0, "GT", "0/1"}' \
    | bgzip > "$OUT/combined.annot_in.vcf.gz"
rm -rf "$AS"
AnnotSV -SVinputFile "$OUT/combined.annot_in.vcf.gz" -SVinputInfo 1 -genomeBuild GRCh38 \
    -annotationsDir "$SVREF" -annotationMode full -outputFile combined.annotsv.tsv -outputDir "$AS"
(cd "$AS" && python "$STEP/extract_all_annotations.py")
