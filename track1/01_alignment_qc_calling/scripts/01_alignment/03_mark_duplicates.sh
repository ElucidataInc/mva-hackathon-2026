#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
mkdir -p "$QC/alignment"
G=$A/$SAMPLE.qgroup.bam

{ printf '@HD\tVN:1.6\tSO:unsorted\tGO:query\n'
  run align samtools view -H --no-PG "$A/$SAMPLE.${LANES%%[[:space:]]*}.bam" | grep -v -e '^@HD' -e '^@RG'
  for L in $LANES; do run align samtools view -H --no-PG "$A/$SAMPLE.$L.bam" | grep '^@RG'; done
} > "$A/qgroup.header.sam"
run align samtools cat -h "$A/qgroup.header.sam" -o "$G" $(for L in $LANES; do echo "$A/$SAMPLE.$L.bam"; done)

run gatk "$GATK" --java-options "-Xmx48g" MarkDuplicates -I "$G" -O "$A/$SAMPLE.md.unsorted.bam" \
    -M "$QC/alignment/$SAMPLE.md.metrics.txt" --ASSUME_SORT_ORDER queryname \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --TMP_DIR "$A/tmp"
