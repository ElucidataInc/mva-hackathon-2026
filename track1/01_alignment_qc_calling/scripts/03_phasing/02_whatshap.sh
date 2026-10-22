#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen -e align)"
export REF PH FINAL_BAM RUN_DIR

xargs -P $((THREADS / 2)) -I{} sh -c 'whatshap phase --reference "$REF" --indels -o "$PH/chrom/{}.phased.vcf.gz" "$PH/chrom/{}.vcf.gz" "$FINAL_BAM" 2> "$RUN_DIR/whatshap.{}.log"' <<< "$CHROMS"
