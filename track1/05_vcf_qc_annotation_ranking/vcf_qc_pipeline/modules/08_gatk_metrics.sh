#!/usr/bin/env bash
# Module 08 – GATK/Picard CollectVariantCallingMetrics (needs dbSNP + sequence dictionary).

mod_gatk() {
    local d="${OUTDIR}/08_gatk_metrics"; mkdir -p "$d"
    have gatk || { mnote "gatk not found"; return 3; }
    [[ -n "$DBSNP" && -n "$REF" ]] || { mnote "needs --dbsnp and --ref"; return 3; }
    local dict="${REF%.*}.dict"; [[ "$REF" == *.gz ]] && dict="${REF%.*.*}.dict"
    [[ -f "$dict" ]] || { mnote "sequence dictionary $dict missing"; return 3; }
    gatk CollectVariantCallingMetrics -I "$VCF_IN" --DBSNP "$DBSNP" -SD "$dict" -O "$d/cvcm" \
        ${TARGETS_INTERVAL_LIST:+-TI "$TARGETS_INTERVAL_LIST"} --THREAD_COUNT "$THREADS" >"$d/cvcm.log" 2>&1
}
