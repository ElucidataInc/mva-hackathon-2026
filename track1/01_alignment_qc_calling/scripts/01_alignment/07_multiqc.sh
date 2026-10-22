#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
C=$QC/custom
mkdir -p "$C"

awk -F'\t' -v OFS='\t' -v s="$SAMPLE" '
    NR == 1 { n = NF; for (i = 3; i <= NF; i++) x[i] = substr($i, 7); next }
    $1 ~ /^(autosomes|chrX|chrY)$/ { for (i = 3; i <= n; i++) { h = h OFS $1 " >=" x[i] " (%)"; v = v OFS $i } }
    END { print "# plot_type: \x27generalstats\x27"; print "Sample" h; print s v }
' "$QC/coverage/$SAMPLE.breadth.tsv" > "$C/breadth_mqc.tsv"

awk -v OFS='\t' -v s="$SAMPLE" '
    { h = h OFS $1; v = v OFS $2 }
    END { print "# plot_type: \x27generalstats\x27"; print "Sample" h; print s v }
' "$QC/identity/$SAMPLE.identity.tsv" > "$C/identity_mqc.tsv"

run qc multiqc -f --cl-config "extra_fn_clean_exts: [.qgroup]" -n mva_qc_report -o "$QC/multiqc" "$QC"/{fastqc,alignment,coverage,custom}
