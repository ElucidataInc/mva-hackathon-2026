#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
run align samtools sort -@ $((THREADS / 2)) -m 3G -T "$A/tmp/md" -o "$A/md.part.bam" "$A/$SAMPLE.md.unsorted.bam"
run align samtools quickcheck "$A/md.part.bam"
mv "$A/md.part.bam" "$FINAL_BAM"
run align samtools index -@ $THREADS "$FINAL_BAM"
(cd "$A" && md5sum "$SAMPLE.md.bam" "$SAMPLE.md.bam.bai" > "$SAMPLE.md.bam.md5")
