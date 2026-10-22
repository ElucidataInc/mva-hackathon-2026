#!/usr/bin/env bash
# =============================================================================
#  VCF QC workflow – modular shell pipeline producing an HTML QC report.
#
#  Modules (run in this order; select with --modules a,b,c):
#    prepare      00  tool inventory, bgzip/index/sort copy of the input
#    inspect      01  size, checksum, header, samples, records per contig
#    integrity    02  REF concordance (bcftools norm -c w), GATK ValidateVariants
#    stats        03  bcftools stats -s - (+ PASS-only) and plot-vcfstats
#    profile      04  Python single-pass profiler (always available)
#    sampleqc     05  vcftools, +guess-ploidy, roh, +mendelian2, somalier, peddy, kinship/PCA/sex/contamination
#    annotation   06  dbSNP / gnomAD (bcftools annotate) + VEP or SnpEff
#    concordance  07  hap.py / RTG vcfeval against a truth set
#    gatk         08  GATK CollectVariantCallingMetrics
#    report       09  HTML report
# =============================================================================
set -uo pipefail

PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="${PIPE_DIR}/scripts"
# shellcheck source=modules/lib.sh
source "${PIPE_DIR}/modules/lib.sh"
for m in "${PIPE_DIR}"/modules/[0-9][0-9]_*.sh; do source "$m"; done

usage() {
cat <<EOF
Usage: $(basename "$0") -i input.vcf[.gz|.bcf] -o outdir [options]

Required
  -i, --vcf FILE             input VCF / VCF.gz / BCF
  -o, --outdir DIR           output directory

Reference & cohort information (optional – each unlocks more checks)
  -r, --ref FASTA            reference FASTA (+ .fai; .dict for GATK)    REF concordance, indel context, GATK
      --ped FILE             PED file                                     Mendelian errors, peddy, sex
      --sex FILE             sample<TAB>sex (M/F/1/2)                     sex check
      --tumor NAME           tumor sample name                            somatic section
      --normal NAME          normal sample name
      --targets BED          restrict bcftools stats / hap.py to targets (exome/panel)
      --targets-interval-list FILE   Picard interval list for CollectVariantCallingMetrics

Annotation / known-vs-novel
      --dbsnp VCF.gz         dbSNP (indexed)                              rsIDs, GATK metrics, ValidateVariants IDs
      --gnomad VCF.gz        gnomAD sites (indexed)                       known-vs-novel, AF bins
      --gnomad-af-tag TAG    INFO tag with AF in the gnomAD file (default: AF)
      --vep-cache DIR        Ensembl VEP cache dir                        consequences / impact
      --vep-assembly NAME    GRCh38 (default) or GRCh37
      --snpeff-db NAME       SnpEff database (used if VEP not configured)

Relatedness / ancestry
      --somalier-sites VCF   somalier sites file (e.g. sites.hg38.vcf.gz)
      --somalier-labels TSV  1000G ancestry labels (ancestry-labels-1kg.tsv)
      --somalier-ref-dir DIR directory of 1000G *.somalier files

Truth-set concordance
      --truth VCF.gz         truth VCF (e.g. GIAB)
      --truth-bed BED        high-confidence regions
      --rtg-sdf DIR          RTG SDF of the reference (for vcfeval)
      --query-sample NAME    sample in the input to compare

Run control
  -t, --threads N            threads (default 4)
      --bin-size BP          bin size for density / depth maps (default 1000000)
      --max-gt-sites N       sites kept for kinship / PCA / ROH fallback (default 100000)
      --max-records N        profile only the first N records (quick look)
      --modules LIST         comma list of modules to run (default: all)
      --title TEXT           report title
      --report-name FILE     report filename inside outdir (default vcf_qc_report.html)
  -h, --help
EOF
}

# ------------------------------------------------------------------ defaults
VCF_ORIG="" OUTDIR="" REF="" PED="" SEX_FILE="" TUMOR="" NORMAL="" TARGETS_BED="" TARGETS_INTERVAL_LIST=""
DBSNP="" GNOMAD="" GNOMAD_AF_TAG="AF" VEP_CACHE="" VEP_ASSEMBLY="GRCh38" SNPEFF_DB=""
SOMALIER_SITES="" SOMALIER_LABELS="" SOMALIER_REF_DIR=""
TRUTH="" TRUTH_BED="" RTG_SDF="" QUERY_SAMPLE=""
THREADS=4 BIN_SIZE=1000000 MAX_GT_SITES=100000 MAX_RECORDS="" MODULES="all"
TITLE="VCF Quality Control Report" REPORT_NAME="vcf_qc_report.html"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--vcf) VCF_ORIG="$2"; shift ;;
        -o|--outdir) OUTDIR="$2"; shift ;;
        -r|--ref) REF="$2"; shift ;;
        --ped) PED="$2"; shift ;;
        --sex) SEX_FILE="$2"; shift ;;
        --tumor) TUMOR="$2"; shift ;;
        --normal) NORMAL="$2"; shift ;;
        --targets) TARGETS_BED="$2"; shift ;;
        --targets-interval-list) TARGETS_INTERVAL_LIST="$2"; shift ;;
        --dbsnp) DBSNP="$2"; shift ;;
        --gnomad) GNOMAD="$2"; shift ;;
        --gnomad-af-tag) GNOMAD_AF_TAG="$2"; shift ;;
        --vep-cache) VEP_CACHE="$2"; shift ;;
        --vep-assembly) VEP_ASSEMBLY="$2"; shift ;;
        --snpeff-db) SNPEFF_DB="$2"; shift ;;
        --somalier-sites) SOMALIER_SITES="$2"; shift ;;
        --somalier-labels) SOMALIER_LABELS="$2"; shift ;;
        --somalier-ref-dir) SOMALIER_REF_DIR="$2"; shift ;;
        --truth) TRUTH="$2"; shift ;;
        --truth-bed) TRUTH_BED="$2"; shift ;;
        --rtg-sdf) RTG_SDF="$2"; shift ;;
        --query-sample) QUERY_SAMPLE="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        --bin-size) BIN_SIZE="$2"; shift ;;
        --max-gt-sites) MAX_GT_SITES="$2"; shift ;;
        --max-records) MAX_RECORDS="$2"; shift ;;
        --modules) MODULES="$2"; shift ;;
        --title) TITLE="$2"; shift ;;
        --report-name) REPORT_NAME="$2"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

[[ -n "$VCF_ORIG" && -n "$OUTDIR" ]] || { usage; exit 1; }
[[ -f "$VCF_ORIG" ]] || { echo "Input not found: $VCF_ORIG" >&2; exit 1; }
PYTHON="${PYTHON:-python3}"
have "$PYTHON" || { echo "python interpreter not found: $PYTHON (set PYTHON=/path/to/python)" >&2; exit 1; }
if ! py_err=$("$PYTHON" -c "import numpy, matplotlib" 2>&1); then
    echo "ERROR: $PYTHON ($(command -v "$PYTHON")) cannot import numpy + matplotlib:" >&2
    echo "  ${py_err##*$'\n'}" >&2
    echo "Fix: $PYTHON -m pip install numpy matplotlib   (or activate the env that has them," >&2
    echo "     or point the pipeline at it: PYTHON=/path/to/env/bin/python ./run_vcf_qc.sh ...)" >&2
    exit 1
fi
for f in REF PED SEX_FILE DBSNP GNOMAD TRUTH TRUTH_BED TARGETS_BED SOMALIER_SITES SOMALIER_LABELS; do
    v="${!f}"; [[ -z "$v" || -e "$v" ]] || { echo "$f file not found: $v" >&2; exit 1; }
done

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"
VCF_ORIG="$(cd "$(dirname "$VCF_ORIG")" && pwd)/$(basename "$VCF_ORIG")"
LOGDIR="${OUTDIR}/logs"; mkdir -p "$LOGDIR"
printf 'module\tfunction\tstatus\tseconds\tnote\n' >"${LOGDIR}/module_status.tsv"
VCF_IN="$VCF_ORIG"
export SCRIPTS OUTDIR LOGDIR PYTHON

{
    echo "date:        $(date)"
    echo "host:        $(hostname)"
    echo "pipeline:    ${PIPE_DIR}"
    echo "python:      $(command -v "$PYTHON") ($("$PYTHON" --version 2>&1))"
    echo "command:     $0 $*"
    for v in VCF_ORIG REF PED SEX_FILE TUMOR NORMAL TARGETS_BED DBSNP GNOMAD GNOMAD_AF_TAG VEP_CACHE VEP_ASSEMBLY \
             SNPEFF_DB SOMALIER_SITES TRUTH TRUTH_BED RTG_SDF QUERY_SAMPLE THREADS BIN_SIZE MAX_GT_SITES MODULES; do
        printf '%-15s %s\n' "$v:" "${!v}"
    done
} >"${LOGDIR}/params.txt"

log "VCF QC pipeline: $VCF_ORIG -> $OUTDIR"
T0=$(date +%s)

run_module prepare mod_prepare
[[ -e "${OUTDIR}/00_prepare/input.use.vcf.gz" ]] && VCF_IN="${OUTDIR}/00_prepare/input.use.vcf.gz"
log "tool-module input: $VCF_IN"

run_module inspect     mod_inspect
run_module integrity   mod_integrity
run_module stats       mod_stats
run_module profile     mod_profile
run_module sampleqc    mod_sampleqc
run_module annotation  mod_annotation
run_module concordance mod_concordance
run_module gatk        mod_gatk
run_module report      mod_report

log "finished in $(( $(date +%s) - T0 ))s"
column -t -s $'\t' "${LOGDIR}/module_status.tsv" >&2 2>/dev/null || cat "${LOGDIR}/module_status.tsv" >&2
[[ -f "${OUTDIR}/${REPORT_NAME}" ]] && log "report: ${OUTDIR}/${REPORT_NAME}"
