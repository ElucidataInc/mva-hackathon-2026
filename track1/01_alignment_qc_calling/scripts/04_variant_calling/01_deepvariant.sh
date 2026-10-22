#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
DV=$RESULTS/deepvariant
mkdir -p "$DV/tmp"

sudo docker run --rm --user "$(id -u):$(id -g)" -v "$PROJECT_ROOT:$PROJECT_ROOT" google/deepvariant:1.10.0 \
    /opt/deepvariant/bin/run_deepvariant --model_type WGS --ref "$REF" --reads "$FINAL_BAM" \
    --output_vcf "$DV/$SAMPLE.deepvariant.vcf.gz" --output_gvcf "$DV/$SAMPLE.deepvariant.g.vcf.gz" \
    --intermediate_results_dir "$DV/tmp" --logging_dir "$RUN_DIR" --num_shards "$THREADS"
