#!/usr/bin/env bash
# Module 04 – single-pass Python profiler (no external tools needed): classes, sort/dup/REF checks,
# per-chromosome & per-sample metrics, density/depth bins, somatic metrics, genotype matrix.

mod_profile() {
    local d="${OUTDIR}/04_profile"; mkdir -p "$d"
    local in="$VCF_ORIG"
    # the profiler reads text VCF(.gz); BCF must be decoded first. Use the ORIGINAL file otherwise,
    # so sort-order problems in the input are reported faithfully.
    [[ "$VCF_ORIG" == *.bcf ]] && in="$VCF_IN"
    "$PYTHON" "${SCRIPTS}/vcf_profile.py" --vcf "$in" --out "$d" \
        ${REF:+--ref "$REF"} ${PED:+--ped "$PED"} ${TUMOR:+--tumor "$TUMOR"} ${NORMAL:+--normal "$NORMAL"} \
        --bin-size "$BIN_SIZE" --max-gt-sites "$MAX_GT_SITES" ${MAX_RECORDS:+--max-records "$MAX_RECORDS"}
}
