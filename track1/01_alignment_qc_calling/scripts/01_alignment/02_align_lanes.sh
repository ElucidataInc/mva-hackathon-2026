#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"

for L in $LANES; do
    out=$A/$SAMPLE.$L.bam
    [[ -s $out ]] && continue
    rg=$(awk -F'\t' -v l=$L '$1 == l {print $5}' "$RG_TSV")
    run align bash -c "set -o pipefail
        bwa mem -t $((THREADS - 4)) -K 100000000 -Y -R '$rg' '$REF' '$(fq $L R1)' '$(fq $L R2)' 2> '$RUN_DIR/bwa.$L.log' |
        samtools view -@ 4 -b -o '$A/$L.part.bam' -"
    run align samtools quickcheck "$A/$L.part.bam"
    mv "$A/$L.part.bam" "$out"
done
