#!/usr/bin/env bash
# Shared helpers for the VCF QC workflow. Sourced by run_vcf_qc.sh.

log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
warn() { printf '[%s] WARNING: %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
have() { command -v "$1" >/dev/null 2>&1; }

MODULE_NOTE=""

# run_module <id> <function>
# Return codes from module functions: 0 = OK, 3 = SKIPPED, anything else = FAILED.
run_module() {
    local id="$1" fn="$2" t0 rc status
    if [[ -n "${MODULES:-}" && ",${MODULES}," != *",${id},"* && ",${MODULES}," != *",all,"* ]]; then
        return 0
    fi
    MODULE_NOTE=""
    CURRENT_MODULE="$id"
    rm -f "${LOGDIR}/${id}.note"
    log "=== module ${id} (${fn}) ==="
    t0=$(date +%s)
    ( set -o pipefail; "$fn" ) >>"${LOGDIR}/${id}.log" 2>&1
    rc=$?
    # notes are written by the subshell to a side file
    [[ -f "${LOGDIR}/${id}.note" ]] && MODULE_NOTE=$(paste -sd';' "${LOGDIR}/${id}.note")
    case $rc in
        0) status=OK ;;
        3) status=SKIPPED ;;
        *) status=FAILED ;;
    esac
    printf '%s\t%s\t%s\t%s\t%s\n' "$id" "$fn" "$status" "$(( $(date +%s) - t0 ))" "${MODULE_NOTE:-}" >>"${LOGDIR}/module_status.tsv"
    log "    -> ${status} ($(( $(date +%s) - t0 ))s) ${MODULE_NOTE:-}"
    [[ $status == FAILED ]] && warn "module ${id} failed – see ${LOGDIR}/${id}.log"
    return 0
}

# Inside a module (subshell): record a note visible in the report's module log
mnote() { echo "$*" >>"${LOGDIR}/${CURRENT_MODULE}.note"; echo "NOTE: $*"; }

# Pick the X / Y contig names used by this VCF
contig_name() {   # contig_name X|Y
    local want="$1"
    local hdr
    hdr=$(bcf_header | grep '^##contig' | sed -E 's/.*ID=([^,>]+).*/\1/')
    for c in "chr${want}" "${want}"; do
        grep -qx "$c" <<<"$hdr" && { echo "$c"; return; }
    done
    echo "chr${want}"
}

bcf_header() {
    if have bcftools; then bcftools view -h "$VCF_IN"
    else (gzip -cdf "$VCF_IN" 2>/dev/null || cat "$VCF_IN") | awk '/^#/{print; next}{exit}'
    fi
}

n_samples() {
    bcf_header | awk -F'\t' '/^#CHROM/{print (NF>9?NF-9:0); exit}'
}

# genome tag for bcftools +guess-ploidy / peddy from the profile or header
genome_tag() {
    local b
    b=$(bcf_header | grep -m1 -E '^##contig=<ID=(chr)?1,' | grep -oE 'length=[0-9]+' | cut -d= -f2)
    local chr=""; bcf_header | grep -q '^##contig=<ID=chr1,' && chr=1
    case "$b" in
        248956422) [[ -n $chr ]] && echo hg38 || echo b38 ;;
        249250621) [[ -n $chr ]] && echo hg19 || echo b37 ;;
        *) echo "" ;;
    esac
}
