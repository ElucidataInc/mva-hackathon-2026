#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
O=$QC/alignment
run align samtools flagstat -@ $THREADS "$FINAL_BAM" > "$O/$SAMPLE.flagstat.txt"
run align samtools idxstats "$FINAL_BAM" > "$O/$SAMPLE.idxstats.txt"
run align samtools stats -@ $THREADS -r "$REF" "$FINAL_BAM" > "$O/$SAMPLE.stats.txt"
