#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
run qc bcftools concat -Oz -o "$PH/$SAMPLE.phased.vcf.gz" $(for c in $CHROMS; do echo "$PH/chrom/$c.phased.vcf.gz"; done)
run qc bcftools index -t "$PH/$SAMPLE.phased.vcf.gz"
run align whatshap stats --tsv "$PH/$SAMPLE.phase_stats.tsv" "$PH/$SAMPLE.phased.vcf.gz" > "$PH/$SAMPLE.phase_stats.txt"
