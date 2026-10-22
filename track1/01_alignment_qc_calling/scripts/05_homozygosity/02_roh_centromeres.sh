#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
R=$RESULTS/homozygosity

{
    awk -v OFS='\t' -v G="$(awk '$1 ~ /^chr[0-9]+$/ {s += $2} END {print s}' "$REF.fai")" '
        FNR == NR {if (!($2 in s) || $3 < s[$2]) s[$2] = $3; if ($4 > e[$2]) e[$2] = $4; next}
        $1 == "RG" && $6 >= 1e6 {
            o = 0; if ($3 in s) {a = $4 > s[$3] ? $4 : s[$3]; b = $5 < e[$3] ? $5 : e[$3]; if (b > a) o = b - a}
            if (!o) l += $6
            printf "%s\t%d\t%d\t%.2f\t%.0f\t%s:%d-%d\t%.2f\n", $3, $4, $5, $6 / 1e6, $7 / ($6 / 1e6), $3, s[$3], e[$3], o / 1e6}
        END {printf "runs_1Mb_no_centromere_Mb\t%.2f\nF_ROH_1Mb_no_centromere\t%.4f\n", l / 1e6, l / G > "/dev/stderr"}' \
        <(zcat "$PROJECT_ROOT/references/ucsc_hg38/centromeres.txt.gz") "$R/$SAMPLE.roh.txt" | sort -k7,7g -k1,1V
} > "$R/roh_centromere.tsv" 2>&1
