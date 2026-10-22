#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"

ls "$REF"{,.fai,.amb,.ann,.bwt,.pac,.sa} "${REF%.fna}.dict" "$GATK" "$FQ_DIR"/*.fastq.gz > /dev/null
df -h "$A" | tail -1

for L in $LANES; do
    head -n 400000 < <(zcat "$(fq $L R1)" 2>/dev/null) | awk -v L=$L -v s=$SAMPLE '
        NR == 1 {split(substr($1, 2), h, ":")}
        NR == 2 {len = length($0)}
        NR % 4 == 1 {n = split($2, a, ":"); c[a[n]]++; t++}
        END {for (k in c) if (c[k] > m) {m = c[k]; b = k}
             printf "%s\t%s\t%.4f\t%d\t@RG\\tID:%s.%s.%s\\tPU:%s.%s.%s\\tSM:%s\\tLB:%s_LIB1\\tPL:ILLUMINA\n", L, b, m / t, len, s, h[3], h[4], h[3], h[4], b, s, s}'
done > "$RG_TSV"
cp "$RG_TSV" "$RUN_DIR/"

[[ $(cut -f2 "$RG_TSV" | sort -u | wc -l) -eq 1 ]] || { echo "lanes differ in dominant index" >&2; exit 1; }
awk '$3 < 0.9 {exit 1}' "$RG_TSV" || { echo "dominant index below 90% of reads" >&2; exit 1; }

set +o pipefail
{ run qc fastqc --version; run qc multiqc --version; run qc samtools --version | head -1
  run align bwa 2>&1 | grep Version; run align mosdepth --version
  run gatk "$GATK" --version 2>&1 | grep -E 'GATK\)|Picard'; java -version 2>&1 | head -1; } > "$RUN_DIR/tool_versions.txt"
