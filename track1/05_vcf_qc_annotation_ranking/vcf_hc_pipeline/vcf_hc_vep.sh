#!/usr/bin/env bash
# =============================================================================
#  Single-sample human VCF -> QC -> normalisation -> filtering -> VEP annotation
#  -> high-confidence variant set (VCF + TSV) with a QC summary.
#
#  Steps (each writes to its own numbered folder in --outdir):
#    01_prepare     check tools/inputs, select sample, match contig names to the FASTA,
#                   keep primary chromosomes, sort, index, raw bcftools stats
#    02_normalize   REF check vs FASTA, split multi-allelics, left-align + trim,
#                   remove exact duplicates, '*' / symbolic / non-variant records
#    03_filter      soft-filter (labels kept in FILTER for auditing):
#                   caller FILTER, QUAL, GQ, DP (min + max), allele balance,
#                   GATK hard filters (when no VQSR/CNN), exclude/include BEDs;
#                   then keep PASS -> high-confidence set
#    04_vep         Ensembl VEP (offline cache or GTF), optional plugins + ClinVar
#    05_tables      one-row-per-variant/gene TSV (MANE/canonical pick) and an
#                   optional "prioritised" subset (rare + damaging / ClinVar P/LP)
#    06_qc          step-by-step counts, Ts/Tv, Het/Hom, known rate, filter breakdown,
#                   bcftools stats before/after, summary.md
# =============================================================================
set -Eeuo pipefail

VERSION="1.2"
usage() {
cat <<EOF
Usage: $(basename "$0") -i sample.vcf.gz -r GRCh38.fa -o outdir --vep-cache ~/.vep [options]

Required
  -i, --vcf FILE            input VCF/VCF.gz/BCF (single sample, or use --sample)
  -r, --ref FASTA           reference FASTA used for calling (with .fai)
  -o, --outdir DIR          output directory
  One of:
      --vep-cache DIR       VEP offline cache directory (recommended)
      --vep-gtf FILE        bgzipped+tabix-indexed GTF (no-cache mode; fewer annotations)
      --skip-vep            stop after filtering

Input options
      --sample NAME         sample to extract from a multi-sample VCF
      --assembly NAME       GRCh38 (default) or GRCh37
      --keep-all-contigs    keep alt/decoy/unplaced contigs (default: 1-22, X, Y, M only)
      --exome               exome/panel data: disables the automatic max-depth filter

Filter thresholds (defaults in brackets)
      --min-qual N          minimum QUAL [20]
      --min-gq N            minimum FORMAT/GQ [20]
      --min-dp N            minimum depth [10]
      --max-dp N            maximum depth; 0 = auto (3 x median depth, WGS only) [0]
      --max-dp-factor F     multiplier for the automatic max depth [3]
      --het-ab-min F        minimum ALT fraction for heterozygous calls [0.2]
      --het-ab-max F        maximum ALT fraction for heterozygous calls [0.8]
      --homalt-ab-min F     minimum ALT fraction for hom-alt calls [0.9]
      --gatk-hard-filters auto|yes|no
                            GATK best-practice hard filters on QD/FS/SOR/MQ/rank sums;
                            auto = only if the caller did not filter (no VQSR/CNN) [auto]
      --keep-caller-filters treat non-PASS FILTER values from the caller as PASS (not advised)
      --exclude-bed BED     regions to exclude (e.g. ENCODE blacklist, low-complexity, segdups)
      --include-bed BED     regions to require (e.g. GIAB high-confidence, exome targets)
      --bed-padding N       pad --include-bed intervals by N bp each side (exome: 50-100) [0]

VEP options
      --vep-fork N          VEP --fork [= --threads]
      --vep-plugins-dir DIR VEP plugins directory
      --clinvar VCF.gz      ClinVar VCF: normalised like the sample, then matched by exact allele
                            (adds INFO/CLNSIG, CLNREVSTAT, CLNDN, ALLELEID)
      --cadd-snv FILE  --cadd-indel FILE        CADD plugin files
      --spliceai-snv FILE --spliceai-indel FILE SpliceAI plugin files
      --revel FILE          REVEL plugin file
      --alphamissense FILE  AlphaMissense plugin file
      --annotate-all        annotate all normalised variants (with FILTER labels), not just PASS
      --max-af F            population AF cut-off for the prioritised subset [0.01]

Other
  -t, --threads N           threads [4]
      --resume              reuse an existing VEP output in --outdir if it is newer than its input
      --cleanup             delete large intermediate files at the end
      --force-build         continue even if VCF and FASTA contig lengths disagree
      --qc-report           also run the HTML QC report (../vcf_qc_pipeline/run_vcf_qc.sh) on the final VCF
  -h, --help
EOF
}

# ------------------------------------------------------------------ defaults
VCF="" REF="" OUT="" VEP_CACHE="" VEP_GTF="" SKIP_VEP=0 SAMPLE="" ASSEMBLY="GRCh38" KEEP_ALL_CONTIGS=0 EXOME=0
MIN_QUAL=20 MIN_GQ=20 MIN_DP=10 MAX_DP=0 MAX_DP_FACTOR=3 HET_AB_MIN=0.2 HET_AB_MAX=0.8 HOMALT_AB_MIN=0.9
GATK_HF=auto KEEP_CALLER_FILTERS=0 EXCLUDE_BED="" INCLUDE_BED=""
VEP_FORK="" VEP_PLUGINS_DIR="" CLINVAR="" CADD_SNV="" CADD_INDEL="" SPLICEAI_SNV="" SPLICEAI_INDEL="" REVEL="" ALPHAMISSENSE=""
ANNOTATE_ALL=0 MAX_AF=0.01 THREADS=4 QC_REPORT=0 BED_PAD=0 RESUME=0 CLEANUP=0 FORCE_BUILD=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--vcf) VCF="$2"; shift ;;
        -r|--ref) REF="$2"; shift ;;
        -o|--outdir) OUT="$2"; shift ;;
        --vep-cache) VEP_CACHE="$2"; shift ;;
        --vep-gtf) VEP_GTF="$2"; shift ;;
        --skip-vep) SKIP_VEP=1 ;;
        --sample) SAMPLE="$2"; shift ;;
        --assembly) ASSEMBLY="$2"; shift ;;
        --keep-all-contigs) KEEP_ALL_CONTIGS=1 ;;
        --exome) EXOME=1 ;;
        --min-qual) MIN_QUAL="$2"; shift ;;
        --min-gq) MIN_GQ="$2"; shift ;;
        --min-dp) MIN_DP="$2"; shift ;;
        --max-dp) MAX_DP="$2"; shift ;;
        --max-dp-factor) MAX_DP_FACTOR="$2"; shift ;;
        --het-ab-min) HET_AB_MIN="$2"; shift ;;
        --het-ab-max) HET_AB_MAX="$2"; shift ;;
        --homalt-ab-min) HOMALT_AB_MIN="$2"; shift ;;
        --gatk-hard-filters) GATK_HF="$2"; shift ;;
        --keep-caller-filters) KEEP_CALLER_FILTERS=1 ;;
        --exclude-bed) EXCLUDE_BED="$2"; shift ;;
        --include-bed) INCLUDE_BED="$2"; shift ;;
        --vep-fork) VEP_FORK="$2"; shift ;;
        --vep-plugins-dir) VEP_PLUGINS_DIR="$2"; shift ;;
        --clinvar) CLINVAR="$2"; shift ;;
        --cadd-snv) CADD_SNV="$2"; shift ;;
        --cadd-indel) CADD_INDEL="$2"; shift ;;
        --spliceai-snv) SPLICEAI_SNV="$2"; shift ;;
        --spliceai-indel) SPLICEAI_INDEL="$2"; shift ;;
        --revel) REVEL="$2"; shift ;;
        --alphamissense) ALPHAMISSENSE="$2"; shift ;;
        --annotate-all) ANNOTATE_ALL=1 ;;
        --max-af) MAX_AF="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        --qc-report) QC_REPORT=1 ;;
        --bed-padding) BED_PAD="$2"; shift ;;
        --resume) RESUME=1 ;;
        --cleanup) CLEANUP=1 ;;
        --force-build) FORCE_BUILD=1 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

# ------------------------------------------------------------------ helpers
log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" | tee -a "${OUT}/logs/pipeline.log" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date '+%H:%M:%S')" "$*" | tee -a "${OUT}/logs/pipeline.log" >&2; exit 1; }
warn() { printf '[%s] WARNING: %s\n' "$(date '+%H:%M:%S')" "$*" | tee -a "${OUT}/logs/pipeline.log" >&2
         echo "- $*" >>"${OUT}/logs/warnings.txt"; }
have() { command -v "$1" >/dev/null 2>&1; }
n_rec() { bcftools index -n "$1" 2>/dev/null || bcftools view -H "$1" | wc -l | tr -d ' '; }
hdr_has() { bcftools view -h "$1" | grep -q "^##$2=<ID=$3,"; }   # hdr_has file INFO|FORMAT ID
count_step() { printf '%s\t%s\t%s\n' "$1" "$(n_rec "$2")" "$3" >>"${OUT}/06_qc/step_counts.tsv"; }

[[ -n "$VCF" && -n "$REF" && -n "$OUT" ]] || { usage; exit 1; }
[[ -n "$VEP_CACHE$VEP_GTF" || $SKIP_VEP -eq 1 ]] || { echo "Give --vep-cache, --vep-gtf or --skip-vep" >&2; exit 1; }
for f in VCF REF; do [[ -f "${!f}" ]] || { echo "$f not found: ${!f}" >&2; exit 1; }; done
for f in EXCLUDE_BED INCLUDE_BED CLINVAR VEP_GTF CADD_SNV CADD_INDEL SPLICEAI_SNV SPLICEAI_INDEL REVEL ALPHAMISSENSE; do
    v="${!f}"; [[ -z "$v" || -f "$v" ]] || { echo "$f not found: $v" >&2; exit 1; }
done
[[ -z "$VEP_CACHE" || -d "$VEP_CACHE" ]] || { echo "VEP cache dir not found: $VEP_CACHE" >&2; exit 1; }

mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
mkdir -p "$OUT"/{logs,01_prepare,02_normalize,03_filter,04_vep,05_tables,06_qc}
: >"${OUT}/logs/warnings.txt"
# never fail silently: report the line and command of any unexpected error
trap 'rc=$?; printf "[%s] ERROR: unexpected failure (exit %s) at line %s: %s\n" "$(date +%H:%M:%S)" "$rc" "$LINENO" "$BASH_COMMAND" | tee -a "${OUT}/logs/pipeline.log" >&2; exit $rc' ERR
printf 'step\trecords\tdescription\n' >"${OUT}/06_qc/step_counts.tsv"
VEP_FORK="${VEP_FORK:-$THREADS}"

# ------------------------------------------------------------------ tools
have bcftools || die "bcftools not found"
BV=$(bcftools --version | head -1 | awk '{print $2}')
awk -v v="$BV" 'BEGIN{split(v,a,"."); exit !(a[1]>1 || (a[1]==1 && a[2]>=13))}' || die "bcftools >= 1.13 required (found $BV)"
have bgzip && have tabix || die "bgzip/tabix (htslib) not found"
if [[ $SKIP_VEP -eq 0 ]]; then have vep || die "vep not found (conda install -c bioconda ensembl-vep) or use --skip-vep"; fi

{
    echo "vcf_hc_vep.sh v${VERSION}  $(date)"
    echo "command: $0 $*"
    for v in VCF REF OUT SAMPLE ASSEMBLY EXOME MIN_QUAL MIN_GQ MIN_DP MAX_DP MAX_DP_FACTOR HET_AB_MIN HET_AB_MAX \
             HOMALT_AB_MIN GATK_HF KEEP_CALLER_FILTERS EXCLUDE_BED INCLUDE_BED BED_PAD VEP_CACHE VEP_GTF CLINVAR MAX_AF THREADS; do
        printf '  %-20s %s\n' "$v" "${!v}"
    done
    echo "bcftools: $BV"
    have vep && echo "vep: $(vep --help 2>/dev/null | grep -m1 -i 'ensembl-vep' | sed 's/^ *//')"
    if [[ -n "$VEP_CACHE" ]]; then
        for i in "$VEP_CACHE"/homo_sapiens*/*_"$ASSEMBLY"/info.txt; do
            [[ -f "$i" ]] && { echo "vep cache: $(dirname "$i")"; grep -E '^(source_|version|assembly|sift|polyphen|COSMIC|ClinVar|gnomAD|dbSNP)' "$i" | sed 's/^/    /'; }
        done
    fi
} >"${OUT}/logs/run_info.txt"
log "output: $OUT  (bcftools $BV)"

# =============================================================================
# 01 PREPARE
# =============================================================================
P="${OUT}/01_prepare"
log "01 prepare: inputs, sample, contig names, primary chromosomes"
if [[ ! -f "${REF}.fai" ]]; then
    have samtools || die "${REF}.fai missing and samtools not found: install samtools and run 'samtools faidx $REF'"
    log "indexing FASTA (samtools faidx)"; samtools faidx "$REF" || die "samtools faidx failed for $REF"
fi

samples=$(bcftools query -l "$VCF")
ns=$(printf '%s\n' "$samples" | grep -c . || true)
[[ "$ns" -ge 1 ]] || die "VCF has no samples (sites-only) – a genotyped single-sample VCF is required"
if [[ -z "$SAMPLE" ]]; then
    [[ "$ns" -eq 1 ]] || die "VCF has $ns samples; choose one with --sample ($(printf '%s ' $samples | cut -c1-200))"
    SAMPLE="$samples"
else
    printf '%s\n' "$samples" | grep -qx "$SAMPLE" || die "sample $SAMPLE not in VCF"
fi
log "  sample: $SAMPLE"

if bcftools view -h "$VCF" | grep -q '<NON_REF>\|<\*>'; then
    warn "input looks like a gVCF: reference blocks are dropped, but QUAL/filters are only meaningful after joint genotyping (GATK GenotypeGVCFs / glnexus)"
fi

# contig naming: map VCF contig names onto the FASTA's naming (1 <-> chr1, MT <-> chrM)
cut -f1 "${REF}.fai" >"$P/fasta_contigs.txt"
bcftools view -h "$VCF" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' >"$P/vcf_header_contigs.txt"
if [[ ! -s "$P/vcf_header_contigs.txt" ]]; then
    bcftools query -f '%CHROM\n' "$VCF" | uniq | sort -u >"$P/vcf_header_contigs.txt"
fi
awk 'NR==FNR{fa[$1]=1; next}
     { c=$1; t="";
       if (c in fa) t=c;
       else if (("chr" c) in fa) t="chr" c;
       else if (c ~ /^chr/ && (substr(c,4) in fa)) t=substr(c,4);
       if (t=="" && (c=="MT" || c=="chrMT") && ("chrM" in fa)) t="chrM";
       if (t=="" && (c=="chrM" || c=="M") && ("MT" in fa)) t="MT";
       if (t!="" && t!=c) print c "\t" t }' "$P/fasta_contigs.txt" "$P/vcf_header_contigs.txt" >"$P/contig_rename.tsv"

# primary contigs present in the FASTA
if [[ $KEEP_ALL_CONTIGS -eq 1 ]]; then
    cp "$P/fasta_contigs.txt" "$P/keep_contigs.txt"
else
    grep -E '^(chr)?([0-9]{1,2}|X|Y|M|MT)$' "$P/fasta_contigs.txt" >"$P/keep_contigs.txt" || true
fi
[[ -s "$P/keep_contigs.txt" ]] || die "no primary chromosomes found in FASTA index"
awk 'NR==FNR{k[$1]=1;next} ($1 in k)' "$P/keep_contigs.txt" "${REF}.fai" >"$P/keep_contigs.fai"

# genome build check: contig lengths in the VCF header must equal those in the FASTA
bcftools view -h "$VCF" | sed -n 's/^##contig=<ID=\([^,>]*\).*length=\([0-9]*\).*/\1\t\2/p' >"$P/vcf_contig_lengths.tsv"
awk -F'\t' -v map="$P/contig_rename.tsv" '
    BEGIN { while ((getline l < map) > 0) { split(l, a, "\t"); m[a[1]] = a[2] } }
    NR==FNR { fa[$1] = $2; next }
    { c = ($1 in m) ? m[$1] : $1; if ((c in fa) && fa[c] != $2) print $1 "\t" $2 "\t" c "\t" fa[c] }' \
    "${REF}.fai" "$P/vcf_contig_lengths.tsv" >"$P/contig_length_mismatch.tsv"
if [[ -s "$P/contig_length_mismatch.tsv" ]]; then
    msg="contig lengths differ between VCF header and FASTA (e.g. $(head -1 "$P/contig_length_mismatch.tsv" | awk '{print $1"="$2" vs "$3"="$4}')): VCF and FASTA are probably different genome builds"
    if [[ $FORCE_BUILD -eq 1 ]]; then warn "$msg (continuing: --force-build)"; else die "$msg. Use the FASTA the VCF was called on, or --force-build"; fi
fi
# mitochondrial contig name in FASTA style (depth / allele-balance filters do not apply to it)
MT=$(grep -m1 -xE '(chr)?(M|MT)' "$P/keep_contigs.txt" || true)

strip=()
if bcftools view -h "$VCF" | grep -q '^##INFO=<ID=CSQ,'; then
    warn "input already has INFO/CSQ – removed so it cannot clash with the new VEP annotation"
    strip=(-x INFO/CSQ)
fi
{
    # -I: do not recompute INFO on subsetting (inconsistent AC/AN in the input would abort);
    # AC/AN/AF are recomputed for the single sample with +fill-tags instead
    bcftools view --threads "$THREADS" -s "$SAMPLE" -I -Ou "$VCF" |
    bcftools +fill-tags -Ou -- -t AC,AN,AF |
    if [[ ${#strip[@]} -gt 0 ]]; then bcftools annotate "${strip[@]}" -Ou; else cat; fi |
    if [[ -s "$P/contig_rename.tsv" ]]; then bcftools annotate --rename-chrs "$P/contig_rename.tsv" -Ou; else cat; fi |
    bcftools view -t "$(paste -sd, "$P/keep_contigs.txt")" -Ou |
    bcftools reheader --fai "$P/keep_contigs.fai" |
    bcftools sort -T "$P/sort_tmp" -Oz -o "$P/prepared.vcf.gz"
} 2>"$P/prepare.log" || { cat "$P/prepare.log" >&2; die "prepare step failed"; }
bcftools index -t -f "$P/prepared.vcf.gz"
[[ -s "$P/contig_rename.tsv" ]] && log "  renamed $(wc -l <"$P/contig_rename.tsv" | tr -d ' ') contigs to match FASTA (e.g. $(head -1 "$P/contig_rename.tsv" | tr '\t' '>'))"
n_in=$( (bcftools index -n "$VCF" 2>/dev/null || bcftools view -H "$VCF" 2>/dev/null | wc -l) | tr -d ' ' || true)
printf 'input\t%s\traw records in input VCF\n' "$n_in" >>"${OUT}/06_qc/step_counts.tsv"
count_step "01_prepared" "$P/prepared.vcf.gz" "single sample, primary chromosomes, sorted"
[[ "$(n_rec "$P/prepared.vcf.gz")" -gt 0 ]] || die "no records left after contig selection – check contig names vs FASTA ($P/contig_rename.tsv)"
bcftools stats --threads "$THREADS" -s - "$P/prepared.vcf.gz" >"${OUT}/06_qc/stats.01_prepared.txt"

# =============================================================================
# 02 NORMALIZE
# =============================================================================
N="${OUT}/02_normalize"
log "02 normalize: REF check, split multi-allelics, left-align, de-duplicate"
# REF concordance report (does not change data)
bcftools norm --threads "$THREADS" -f "$REF" -c w -Ou -o /dev/null "$P/prepared.vcf.gz" 2>"$N/ref_check.log" || true
n_refmm=$(grep -c 'REF_MISMATCH' "$N/ref_check.log" || true)
n_prep=$(n_rec "$P/prepared.vcf.gz")
log "  REF mismatches vs FASTA: $n_refmm of $n_prep"
if [[ "$n_prep" -gt 0 ]] && awk -v a="$n_refmm" -v b="$n_prep" 'BEGIN{exit !(a/b > 0.01)}'; then
    warn "$n_refmm REF mismatches (>1%): the VCF was probably called on a different reference/build than $REF"
fi

{
    # -m -any : split multi-allelic records (one ALT per line)
    # -f/-c x : left-align & trim indels against the reference; drop records whose REF disagrees
    # --force : tolerate FORMAT fields with the wrong number of values (e.g. AD on no-call genotypes)
    bcftools norm --threads "$THREADS" -f "$REF" -c x -m -any --force -Ou "$P/prepared.vcf.gz" |
    # exact duplicates (same CHROM/POS/REF/ALT) created by splitting or present in the input
    bcftools norm -d exact -Ou |
    # spanning deletions '*', symbolic alleles (<NON_REF>, <DEL>, ...), and sites the sample does not carry
    bcftools view -e 'ALT="*" || ALT~"<" || GT!="alt"' -Ou |
    bcftools sort -T "$N/sort_tmp" -Oz -o "$N/normalized.vcf.gz"
} 2>"$N/normalize.log" || { grep -m5 -iE 'error|wrong' "$N/normalize.log" >&2; die "normalisation failed (see $N/normalize.log)"; }
grep -q 'wrong number of fields' "$N/normalize.log" && warn "some records have malformed FORMAT fields (wrong number of values); kept with bcftools norm --force"
bcftools index -t -f "$N/normalized.vcf.gz"
grep -E 'Lines +total' "$N/normalize.log" | head -2 >"$N/norm_summary.txt" || true
count_step "02_normalized" "$N/normalized.vcf.gz" "split, left-aligned, REF-checked, de-duplicated, ALT-carrying only"
bcftools stats --threads "$THREADS" -s - -F "$REF" "$N/normalized.vcf.gz" >"${OUT}/06_qc/stats.02_normalized.txt"

# =============================================================================
# 03 FILTER
# =============================================================================
F="${OUT}/03_filter"
log "03 filter: soft-filter with labels, then keep PASS"
cur="$N/normalized.vcf.gz"

# FORMAT/VAF from AD (consistent allele fraction for every caller)
if hdr_has "$cur" FORMAT AD; then
    bcftools +fill-tags "$cur" -Oz -o "$F/with_vaf.vcf.gz" -- -t FORMAT/VAF 2>"$F/fill_tags.log"
    bcftools index -t -f "$F/with_vaf.vcf.gz"; cur="$F/with_vaf.vcf.gz"; HAS_VAF=1
else
    warn "no FORMAT/AD: allele-balance filters skipped"; HAS_VAF=0
fi

# caller FILTER values present?
bcftools query -f '%FILTER\n' "$cur" | sort | uniq -c | sort -rn >"$F/input_filter_values.txt"
caller_filtered=$(awk '$2!="PASS" && $2!="."{n+=$1} END{print n+0}' "$F/input_filter_values.txt")
hdr=$(bcftools view -h "$cur")

# chrM: depth is ~100x nuclear and heteroplasmy gives any allele fraction, so depth and
# allele-balance filters are not applied there (median depth is also computed without it)
NOT_MT=""; [[ -n "$MT" ]] && NOT_MT="CHROM!=\"${MT}\" && "

# depth source
DPX=""
if hdr_has "$cur" FORMAT DP; then DPX="FMT/DP"; elif hdr_has "$cur" INFO DP; then DPX="INFO/DP"; fi
if [[ -n "$DPX" && "$MAX_DP" == "0" && $EXOME -eq 0 ]]; then
    med=$(bcftools query ${MT:+-e "CHROM=\"$MT\""} -f "$( [[ $DPX == FMT/DP ]] && echo '[%DP]' || echo '%INFO/DP')\n" "$cur" |
          awk '$1!="."' | sort -n | awk '{a[NR]=$1} END{if(NR) print a[int((NR+1)/2)]; else print 0}')
    MAX_DP=$(awk -v m="$med" -v f="$MAX_DP_FACTOR" 'BEGIN{printf "%d", m*f}')
    log "  median depth $med -> max depth filter DP > $MAX_DP"
    echo "median_depth	$med" >"$F/depth.tsv"
fi

# each filter: label<TAB>expression (records matching the expression FAIL)
: >"$F/filters.tsv"
addf() { printf '%s\t%s\t%s\n' "$1" "$2" "$3" >>"$F/filters.tsv"; }
[[ "$MIN_QUAL" != "0" ]] && addf LowQual "QUAL<${MIN_QUAL}" "QUAL < $MIN_QUAL"
hdr_has "$cur" FORMAT GQ && addf LowGQ "FMT/GQ<${MIN_GQ}" "GQ < $MIN_GQ" || warn "no FORMAT/GQ: GQ filter skipped"
if [[ -n "$DPX" ]]; then
    addf LowDP "${DPX}<${MIN_DP}" "depth < $MIN_DP"
    [[ "$MAX_DP" != "0" ]] && addf HighDP "${NOT_MT}${DPX}>${MAX_DP}" "depth > $MAX_DP (collapsed repeats / CNVs; not applied to ${MT:-chrM})"
else
    warn "no depth field: DP filters skipped"
fi
if [[ $HAS_VAF -eq 1 ]]; then
    addf HetAB "${NOT_MT}GT=\"het\" && (FMT/VAF<${HET_AB_MIN} || FMT/VAF>${HET_AB_MAX})" "het with ALT fraction outside ${HET_AB_MIN}-${HET_AB_MAX} (not applied to ${MT:-chrM})"
    addf HomAltAB "${NOT_MT}GT=\"AA\" && FMT/VAF<${HOMALT_AB_MIN}" "hom-alt with ALT fraction < ${HOMALT_AB_MIN} (not applied to ${MT:-chrM})"
fi

# GATK hard filters (https://gatk.broadinstitute.org/hc/en-us/articles/360035890471)
use_hf=0
if [[ "$GATK_HF" == "yes" ]]; then use_hf=1
elif [[ "$GATK_HF" == "auto" ]] && grep -q '^##INFO=<ID=QD,' <<<"$hdr" && \
     ! grep -qE '^##INFO=<ID=(VQSLOD|AS_VQSLOD|CNN_1D|CNN_2D),' <<<"$hdr" && [[ "$caller_filtered" -eq 0 ]]; then
    use_hf=1
fi
if [[ $use_hf -eq 1 ]]; then
    snp=() ind=()
    grep -q '^##INFO=<ID=QD,' <<<"$hdr"             && { snp+=("QD<2"); ind+=("QD<2"); }
    grep -q '^##INFO=<ID=FS,' <<<"$hdr"             && { snp+=("FS>60"); ind+=("FS>200"); }
    grep -q '^##INFO=<ID=SOR,' <<<"$hdr"            && { snp+=("SOR>3"); ind+=("SOR>10"); }
    grep -q '^##INFO=<ID=MQ,' <<<"$hdr"             && snp+=("MQ<40")
    grep -q '^##INFO=<ID=MQRankSum,' <<<"$hdr"      && snp+=("MQRankSum<-12.5")
    grep -q '^##INFO=<ID=ReadPosRankSum,' <<<"$hdr" && { snp+=("ReadPosRankSum<-8"); ind+=("ReadPosRankSum<-20"); }
    j() { local IFS="|"; echo "$*" | sed 's/|/ || /g'; }
    [[ ${#snp[@]} -gt 0 ]] && addf GATK_SNP_HF "TYPE!=\"indel\" && ($(j "${snp[@]}"))" "GATK SNP hard filter: $(j "${snp[@]}")"
    [[ ${#ind[@]} -gt 0 ]] && addf GATK_INDEL_HF "TYPE=\"indel\" && ($(j "${ind[@]}"))" "GATK indel hard filter: $(j "${ind[@]}")"
    log "  GATK hard filters applied"
else
    log "  GATK hard filters not applied (mode=$GATK_HF; caller-filtered records: $caller_filtered)"
fi

# region filters via INFO flags
prep_bed() {   # bgzip + tabix a BED (renaming contigs to FASTA style, optional padding)
    local in="$1" out="$2" pad="${3:-0}"
    ( if [[ "$in" == *.gz ]]; then gzip -cd "$in"; else cat "$in"; fi ) | grep -vE '^(#|track|browser)' |
        awk -v map="$P/contig_rename.tsv" 'BEGIN{ OFS="\t"; while ((getline l < map) > 0) { split(l, a, "\t"); m[a[1]]=a[2] } }
             { if ($1 in m) $1=m[$1]; s=$2-pad; if (s<0) s=0; print $1,s,$3+pad }' pad="$pad" |
        sort -k1,1 -k2,2n | bgzip -c >"$out"
    tabix -f -p bed "$out"
}
if [[ -n "$EXCLUDE_BED" ]]; then
    prep_bed "$EXCLUDE_BED" "$F/exclude.bed.gz"
    bcftools annotate -a "$F/exclude.bed.gz" -c CHROM,FROM,TO --mark-sites +EXCLUDED_REGION -Oz -o "$F/regions1.vcf.gz" "$cur"
    bcftools index -t -f "$F/regions1.vcf.gz"; cur="$F/regions1.vcf.gz"
    addf ExcludedRegion "INFO/EXCLUDED_REGION=1" "in --exclude-bed regions"
fi
if [[ -n "$INCLUDE_BED" ]]; then
    prep_bed "$INCLUDE_BED" "$F/include.bed.gz" "$BED_PAD"
    bcftools annotate -a "$F/include.bed.gz" -c CHROM,FROM,TO --mark-sites +INCLUDED_REGION -Oz -o "$F/regions2.vcf.gz" "$cur"
    bcftools index -t -f "$F/regions2.vcf.gz"; cur="$F/regions2.vcf.gz"
    addf OutsideIncludeBed "INFO/INCLUDED_REGION!=1" "outside --include-bed regions (padding ${BED_PAD} bp)"
fi

# apply all soft filters (labels appended to FILTER, so every reason is kept)
src="$cur"
if [[ $KEEP_CALLER_FILTERS -eq 1 ]]; then
    bcftools annotate -x FILTER -Oz -o "$F/caller_filters_cleared.vcf.gz" "$src"; src="$F/caller_filters_cleared.vcf.gz"
    warn "--keep-caller-filters: caller FILTER values were discarded"
fi
FL=() FX=()
while IFS=$'\t' read -r label expr desc; do FL+=("$label"); FX+=("$expr"); done <"$F/filters.tsv"
chain() {   # chain i : stream through bcftools filter for labels i..end (single pass, no recompression)
    local i=$1
    if [[ $i -ge ${#FL[@]} ]]; then cat; return; fi
    bcftools filter -m + -s "${FL[$i]}" -e "${FX[$i]}" -Ou - | chain $((i + 1))
}
{ bcftools view -Ou "$src" | chain 0 | bcftools view -Oz -o "$F/soft_filtered.vcf.gz"; } 2>"$F/filter.log" ||
    { cat "$F/filter.log" >&2; die "soft filtering failed (see $F/filter.log; expressions in $F/filters.tsv)"; }
bcftools index -t -f "$F/soft_filtered.vcf.gz"

# hard selection: PASS only
bcftools view -f PASS -Oz -o "$F/high_confidence.vcf.gz" "$F/soft_filtered.vcf.gz"
bcftools index -t -f "$F/high_confidence.vcf.gz"
count_step "03_high_confidence" "$F/high_confidence.vcf.gz" "FILTER=PASS after all soft filters"

# filter breakdown: each label, alone and as the only reason
bcftools query -f '%FILTER\n' "$F/soft_filtered.vcf.gz" >"$F/filter_column.txt"
{
    printf 'filter\trecords_failing\trecords_failing_only_this\n'
    tr ';' '\n' <"$F/filter_column.txt" | grep -vx PASS | sort | uniq -c | sort -rn | while read -r n lab; do
        only=$(grep -cx "$lab" "$F/filter_column.txt" || true)
        printf '%s\t%s\t%s\n' "$lab" "$n" "$only"
    done
} >"${OUT}/06_qc/filter_breakdown.tsv"
rm -f "$F/filter_column.txt"
bcftools stats --threads "$THREADS" -s - -F "$REF" "$F/high_confidence.vcf.gz" >"${OUT}/06_qc/stats.03_high_confidence.txt"
log "  high-confidence: $(n_rec "$F/high_confidence.vcf.gz") of $(n_rec "$F/soft_filtered.vcf.gz") normalised variants"

# =============================================================================
# 04 VEP
# =============================================================================
V="${OUT}/04_vep"
FINAL_VCF="$F/high_confidence.vcf.gz"

# ClinVar: normalise exactly like the sample (contig names, split, left-align) and match on
# CHROM/POS/REF/ALT, so a pathogenic entry for a different allele at the same position never transfers
if [[ -n "$CLINVAR" ]]; then
    log "04 ClinVar: normalising and annotating by exact allele"
    bcftools view -h "$CLINVAR" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' >"$V/clinvar_contigs.txt"
    [[ -s "$V/clinvar_contigs.txt" ]] || bcftools query -f '%CHROM\n' "$CLINVAR" | uniq | sort -u >"$V/clinvar_contigs.txt"
    awk 'NR==FNR{fa[$1]=1; next}
         { c=$1; t="";
           if (c in fa) t=c; else if (("chr" c) in fa) t="chr" c; else if (c ~ /^chr/ && (substr(c,4) in fa)) t=substr(c,4);
           if (t=="" && (c=="MT"||c=="chrMT") && ("chrM" in fa)) t="chrM"; if (t=="" && (c=="chrM"||c=="M") && ("MT" in fa)) t="MT";
           if (t!="" && t!=c) print c "\t" t }' "$P/fasta_contigs.txt" "$V/clinvar_contigs.txt" >"$V/clinvar_rename.tsv"
    {
        bcftools view -Ou "$CLINVAR" |
        if [[ -s "$V/clinvar_rename.tsv" ]]; then bcftools annotate --rename-chrs "$V/clinvar_rename.tsv" -Ou; else cat; fi |
        bcftools view -t "$(paste -sd, "$P/keep_contigs.txt")" -Ou |
        bcftools norm -f "$REF" -c x -m -any -Ou |
        bcftools norm -d exact -Ou |
        bcftools sort -T "$V/sort_tmp" -Oz -o "$V/clinvar.norm.vcf.gz"
    } 2>"$V/clinvar_norm.log" || { tail "$V/clinvar_norm.log" >&2; die "ClinVar normalisation failed"; }
    bcftools index -t -f "$V/clinvar.norm.vcf.gz"
    cv_tags=""
    for t in CLNSIG CLNREVSTAT CLNDN ALLELEID CLNSIGCONF; do
        hdr_has "$V/clinvar.norm.vcf.gz" INFO "$t" && cv_tags="${cv_tags:+$cv_tags,}INFO/$t"
    done
    [[ -n "$cv_tags" ]] || die "ClinVar file has none of CLNSIG/CLNREVSTAT/CLNDN/ALLELEID"
    bcftools annotate -a "$V/clinvar.norm.vcf.gz" -c "$cv_tags" -Oz -o "$V/high_confidence.clinvar.vcf.gz" \
        "$F/high_confidence.vcf.gz" 2>"$V/clinvar_annotate.log" || { cat "$V/clinvar_annotate.log" >&2; die "ClinVar annotation failed"; }
    bcftools index -t -f "$V/high_confidence.clinvar.vcf.gz"
    FINAL_VCF="$V/high_confidence.clinvar.vcf.gz"
    log "  ClinVar entries matched: $(bcftools view -H -i 'INFO/CLNSIG!="."' "$FINAL_VCF" | wc -l | tr -d ' ')"
    # ClinVar pathogenic/likely pathogenic variants that FAILED QC are not silently lost:
    # they are listed for orthogonal confirmation (e.g. Sanger / IGV review)
    bcftools annotate -a "$V/clinvar.norm.vcf.gz" -c "$cv_tags" -Ou "$F/soft_filtered.vcf.gz" |
        # exclusion in its own step: inside a combined -i expression, bcftools' !~ does not exclude
        bcftools view -e 'INFO/CLNSIG~"onflicting"' -Ou |
        bcftools view -i 'FILTER!="PASS" && INFO/CLNSIG~"athogenic"' -Ou |
        bcftools query -H -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN[\t%GT$(hdr_has "$F/soft_filtered.vcf.gz" FORMAT GQ && echo '\t%GQ')$(hdr_has "$F/soft_filtered.vcf.gz" FORMAT DP && echo '\t%DP')$(hdr_has "$F/soft_filtered.vcf.gz" FORMAT VAF && echo '\t%VAF')]\n" |
        sed '1s/^# *//; 1s/\[[0-9]*\]//g' >"${OUT}/05_tables/clinvar_plp_failed_qc.tsv"
    n_plp_fail=$(($(wc -l <"${OUT}/05_tables/clinvar_plp_failed_qc.tsv") - 1))
    [[ $n_plp_fail -gt 0 ]] && warn "$n_plp_fail ClinVar pathogenic/likely-pathogenic variant(s) FAILED QC filters – review 05_tables/clinvar_plp_failed_qc.tsv (confirm orthogonally before discarding)"
fi
if [[ $SKIP_VEP -eq 0 ]]; then
    vin="$FINAL_VCF"
    if [[ $ANNOTATE_ALL -eq 1 ]]; then
        vin="$F/soft_filtered.vcf.gz"
        if [[ -n "$CLINVAR" ]]; then
            bcftools annotate -a "$V/clinvar.norm.vcf.gz" -c "$cv_tags" -Oz -o "$V/soft_filtered.clinvar.vcf.gz" "$vin"
            bcftools index -t -f "$V/soft_filtered.clinvar.vcf.gz"; vin="$V/soft_filtered.clinvar.vcf.gz"
        fi
    fi
    log "04 VEP: annotating $(n_rec "$vin") variants ($( [[ -n "$VEP_CACHE" ]] && echo cache || echo GTF ) mode, fork $VEP_FORK)"
    vep_args=(--input_file "$vin" --output_file "$V/annotated.vcf.gz" --vcf --compress_output bgzip --force_overwrite
              --species homo_sapiens --assembly "$ASSEMBLY" --fasta "$REF" --fork "$VEP_FORK" --buffer_size 5000
              --stats_file "$V/vep_summary.html" --warning_file "$V/vep_warnings.txt"
              --flag_pick_allele_gene)
    if [[ -n "$VEP_CACHE" ]]; then
        vv=$(vep --help 2>/dev/null | sed -n 's/.*ensembl-vep *: *\([0-9]*\).*/\1/p' | head -1)
        if [[ -n "$vv" ]] && ! ls -d "$VEP_CACHE"/homo_sapiens*/"${vv}"_"${ASSEMBLY}" >/dev/null 2>&1; then
            warn "no VEP cache for version ${vv}_${ASSEMBLY} in $VEP_CACHE (found: $(ls "$VEP_CACHE"/homo_sapiens* 2>/dev/null | paste -sd' ' -)); use a cache matching the VEP version"
        fi
        # --everything: symbol, biotype, canonical, MANE, HGVS, SIFT, PolyPhen, gnomAD/1000G AFs,
        # existing variants (rsIDs), ClinVar significance, domains, regulatory, etc.
        vep_args+=(--offline --cache --dir_cache "$VEP_CACHE" --everything)
    else
        vep_args+=(--gtf "$VEP_GTF" --symbol --biotype --canonical --numbers --hgvs --variant_class)
    fi
    [[ -n "$VEP_PLUGINS_DIR" ]] && vep_args+=(--dir_plugins "$VEP_PLUGINS_DIR")
    [[ -n "$CADD_SNV" ]] && vep_args+=(--plugin "CADD,snv=${CADD_SNV}${CADD_INDEL:+,indels=${CADD_INDEL}}")
    [[ -n "$SPLICEAI_SNV" ]] && vep_args+=(--plugin "SpliceAI,snv=${SPLICEAI_SNV}${SPLICEAI_INDEL:+,indel=${SPLICEAI_INDEL}}")
    [[ -n "$REVEL" ]] && vep_args+=(--plugin "REVEL,file=${REVEL}")
    [[ -n "$ALPHAMISSENSE" ]] && vep_args+=(--plugin "AlphaMissense,file=${ALPHAMISSENSE}")
    printf '%q ' vep "${vep_args[@]}" >"$V/vep_command.txt"; echo >>"$V/vep_command.txt"
    # fingerprint = VEP input records + VEP command; --resume reuses a finished VEP run with the same fingerprint
    md5() { if have md5sum; then md5sum | cut -d' ' -f1; else command md5 -q; fi; }
    fp=$( { bcftools view -H "$vin"; sed "s#$V/[^ ]*##g" "$V/vep_command.txt"; } | md5 )
    if [[ $RESUME -eq 1 && -s "$V/annotated.vcf.gz" && -f "$V/annotated.fingerprint" && "$(cat "$V/annotated.fingerprint")" == "$fp" ]] &&
       bcftools view -h "$V/annotated.vcf.gz" >/dev/null 2>&1; then
        log "  --resume: VEP input and options unchanged – reusing $V/annotated.vcf.gz"
    else
        rm -f "$V/annotated.fingerprint"
        vep "${vep_args[@]}" >"$V/vep.log" 2>&1 || { tail -20 "$V/vep.log" >&2; die "VEP failed (see $V/vep.log)"; }
        echo "$fp" >"$V/annotated.fingerprint"
    fi
    bcftools index -t -f "$V/annotated.vcf.gz"
    count_step "04_vep_annotated" "$V/annotated.vcf.gz" "VEP-annotated"
    FINAL_VCF="$V/annotated.vcf.gz"
    if [[ $ANNOTATE_ALL -eq 1 ]]; then
        bcftools view -f PASS -Oz -o "$V/annotated.high_confidence.vcf.gz" "$V/annotated.vcf.gz"
        bcftools index -t -f "$V/annotated.high_confidence.vcf.gz"
        FINAL_VCF="$V/annotated.high_confidence.vcf.gz"
    fi
else
    log "04 VEP skipped"
fi

# =============================================================================
# 05 TABLES
# =============================================================================
T="${OUT}/05_tables"
log "05 tables"
# sample columns: GT GQ DP VAF AD (whichever exist)
fmt_q="" fmt_h=""
for x in GT GQ DP VAF AD; do
    hdr_has "$FINAL_VCF" FORMAT "$x" && { fmt_q="${fmt_q}\t%${x}"; fmt_h="${fmt_h}\t${x}"; }
done
info_q="" info_h=""
for x in CLNSIG CLNREVSTAT CLNDN ALLELEID; do
    hdr_has "$FINAL_VCF" INFO "$x" && { info_q="${info_q}\t%INFO/${x}"; info_h="${info_h}\tClinVar_${x}"; }
done
base="%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER${info_q}${fmt_q:+[${fmt_q}]}"
hcols="CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER${info_h}${fmt_h}"
if bcftools view -h "$FINAL_VCF" | grep -q '^##INFO=<ID=CSQ,'; then
    bcftools +split-vep -l "$FINAL_VCF" | cut -f2 >"$T/csq_fields_available.txt"
    want="Consequence IMPACT SYMBOL Gene Feature BIOTYPE MANE_SELECT CANONICAL EXON INTRON HGVSc HGVSp
          Protein_position Amino_acids Existing_variation VARIANT_CLASS SIFT PolyPhen gnomADe_AF gnomADg_AF
          AF MAX_AF MAX_AF_POPS CLIN_SIG SOMATIC PHENO PUBMED
          CADD_PHRED CADD_RAW SpliceAI_pred_SYMBOL SpliceAI_pred_DS_AG SpliceAI_pred_DS_AL SpliceAI_pred_DS_DG
          SpliceAI_pred_DS_DL REVEL am_class am_pathogenicity PICK"
    fields=""
    for w in $want; do grep -qx "$w" "$T/csq_fields_available.txt" && fields="$fields $w"; done
    csq_f=$(for w in $fields; do printf '\\t%%%s' "$w"; done)
    csq_h=$(for w in $fields; do printf '\t%s' "$w"; done)
    sel=()
    grep -qx PICK "$T/csq_fields_available.txt" && sel=(-i 'PICK="1"')
    { printf "${hcols}${csq_h}\n"
      bcftools +split-vep "$FINAL_VCF" -d ${sel[@]+"${sel[@]}"} -f "${base}${csq_f}\n" 2>"$T/split_vep.log" || true
    } >"$T/high_confidence.annotated.tsv"
    # fall back to all transcripts if PICK filtering produced nothing
    if [[ $(wc -l <"$T/high_confidence.annotated.tsv") -le 1 ]]; then
        warn "PICK-based table empty ($(head -1 "$T/split_vep.log")); using the most severe consequence per variant"
        { printf "${hcols}${csq_h}\n"; bcftools +split-vep "$FINAL_VCF" -s worst -f "${base}${csq_f}\n" 2>>"$T/split_vep.log"; } \
            >"$T/high_confidence.annotated.tsv" || die "bcftools +split-vep failed (see $T/split_vep.log)"
    fi
    # prioritised subset (a starting point for review, not a diagnosis):
    #   rare_damaging : max population AF < MAX_AF (or absent) AND HIGH/MODERATE impact
    #   clinvar_plp   : ClinVar (allele-matched) pathogenic/likely pathogenic, not conflicting,
    #                   with assertion criteria (>= 1 review star)
    #   vep_clin_sig  : VEP cache CLIN_SIG pathogenic (co-located, allele-checked by VEP)
    awk -F'\t' -v OFS='\t' -v maxaf="$MAX_AF" '
        function col(n) { return (n in c) ? $(c[n]) : "" }
        NR==1 { for (i=1;i<=NF;i++) c[$i]=i; print $0, "priority_reason"; next }
        { af=-1; n=split("MAX_AF gnomADe_AF gnomADg_AF AF", ks, " ")
          for (k=1;k<=n;k++) { v0=col(ks[k]); if (v0!="" && v0!=".") { split(v0,v,"&"); for (j in v) if (v[j]!="" && v[j]+0>af) af=v[j]+0 } }
          rare = (af<0 || af<maxaf); imp = col("IMPACT")
          cs = tolower(col("ClinVar_CLNSIG")); rs = tolower(col("ClinVar_CLNREVSTAT"))
          plp = (cs ~ /pathogenic/ && cs !~ /conflicting/ && cs !~ /benign/ && rs !~ /no_assertion|no_interpretation|no_classification/ && cs != ".")
          vs = tolower(col("CLIN_SIG"))
          vplp = (vs ~ /(^|&)(likely_)?pathogenic/ && vs !~ /conflicting|benign/)
          r = ""
          if (rare && (imp=="HIGH" || imp=="MODERATE")) r = r "rare_" tolower(imp) ";"
          if (plp) r = r "clinvar_plp;"
          if (vplp) r = r "vep_clin_sig;"
          if (r != "") { sub(/;$/, "", r); print $0, r } }' \
        "$T/high_confidence.annotated.tsv" >"$T/prioritised.tsv"
    log "  $(($(wc -l <"$T/high_confidence.annotated.tsv")-1)) annotated rows; $(($(wc -l <"$T/prioritised.tsv")-1)) prioritised rows"
else
    { printf "${hcols}\n"; bcftools query -f "${base}\n" "$FINAL_VCF"; } >"$T/high_confidence.tsv"
fi
cp "$FINAL_VCF" "${OUT}/final.high_confidence.vcf.gz"
cp "${FINAL_VCF}.tbi" "${OUT}/final.high_confidence.vcf.gz.tbi"

# =============================================================================
# 06 QC SUMMARY
# =============================================================================
Q="${OUT}/06_qc"
log "06 QC summary"
sn() { awk -F'\t' -v k="$2" '$1=="SN" && $3==k {print $4}' "$1"; }
tstv() { awk -F'\t' '$1=="TSTV"{print $5; exit}' "$1"; }
hethom() { awk -F'\t' '$1=="PSC"{ if ($5>0) printf "%.3f", $6/$5; else print "NA"; exit}' "$1"; }
{
    printf 'stage\trecords\tSNPs\tindels\tMNPs\tmultiallelic_sites\tTs/Tv\tHet/HomAlt\n'
    for s in 01_prepared 02_normalized 03_high_confidence; do
        f="$Q/stats.$s.txt"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$(sn "$f" 'number of records:')" "$(sn "$f" 'number of SNPs:')" \
            "$(sn "$f" 'number of indels:')" "$(sn "$f" 'number of MNPs:')" "$(sn "$f" 'number of multiallelic sites:')" \
            "$(tstv "$f")" "$(hethom "$f")"
    done
} >"$Q/metrics.tsv"

# known = rsID in the VCF ID column (caller/dbSNP) or in VEP Existing_variation (cache mode only)
known="NA (needs VEP cache or rsIDs in the ID column)"
TAB=$(ls "$T"/high_confidence*.tsv 2>/dev/null | head -1)
if [[ -n "$TAB" ]]; then
    known=$(awk -F'\t' -v cache="$([[ -n "$VEP_CACHE" && $SKIP_VEP -eq 0 ]] && echo 1 || echo 0)" '
        NR==1 { for (i=1;i<=NF;i++) if ($i=="Existing_variation") k=i; next }
        { key=$1":"$2":"$4":"$5; if (key in s) next; s[key]=1
          t = (length($4)==1 && length($5)==1) ? "SNV" : "indel/other"; n[t]++
          if ($3 ~ /rs[0-9]+/ || (k && $k ~ /rs[0-9]+/)) { r[t]++; any=1 } }
        END { if (!any && !cache) { print "NA (needs VEP cache or rsIDs in the ID column)"; exit }
              out=""; for (t in n) out = out sprintf("%s %.2f%% (%d/%d); ", t, 100*r[t]/n[t], r[t], n[t]); sub(/; $/,"",out); print out }' "$TAB")
fi
have plot-vcfstats && plot-vcfstats -P -p "$Q/plots_high_confidence/" "$Q/stats.03_high_confidence.txt" >/dev/null 2>&1 || true

hc_tstv=$(tstv "$Q/stats.03_high_confidence.txt")
hc_hh=$(hethom "$Q/stats.03_high_confidence.txt")
{
    echo "# High-confidence variant pipeline – summary"
    echo
    echo "- Sample: **$SAMPLE**   Input: \`$VCF\`   Reference: \`$REF\` ($ASSEMBLY)"
    echo "- Date: $(date)"
    echo "- Final VCF: \`final.high_confidence.vcf.gz\`  Table: \`05_tables/$(ls "$T" | grep -m1 -E '^high_confidence.*tsv$')\`"
    echo
    echo "## Variants at each step"
    echo; echo '| step | records | description |'; echo '|---|---:|---|'
    tail -n +2 "$Q/step_counts.tsv" | awk -F'\t' '{printf "| %s | %s | %s |\n",$1,$2,$3}'
    echo
    echo "## Metrics"
    echo; echo '| stage | records | SNPs | indels | MNPs | multiallelic | Ts/Tv | Het/HomAlt |'; echo '|---|---:|---:|---:|---:|---:|---:|---:|'
    tail -n +2 "$Q/metrics.tsv" | awk -F'\t' '{printf "| %s | %s | %s | %s | %s | %s | %s | %s |\n",$1,$2,$3,$4,$5,$6,$7,$8}'
    echo
    echo "- REF mismatches vs FASTA (removed): $n_refmm"
    [[ -n "$MT" ]] && echo "- ${MT}: depth and allele-balance filters not applied (heteroplasmy, very high depth); for mitochondrial analysis use a dedicated caller (e.g. GATK Mutect2 --mitochondria-mode)"
    if [[ -n "$CLINVAR" ]]; then echo "- ClinVar (allele-exact) matches in final set: $(bcftools view -H -i 'INFO/CLNSIG!="."' "${OUT}/final.high_confidence.vcf.gz" | wc -l | tr -d ' ')"; fi
    [[ -f "$T/prioritised.tsv" ]] && echo "- Prioritised rows (\`05_tables/prioritised.tsv\`): $(($(wc -l <"$T/prioritised.tsv")-1)) – $(tail -n +2 "$T/prioritised.tsv" | awk -F'\t' '{print $NF}' | tr ';' '\n' | sort | uniq -c | awk '{printf "%s%s %s", (NR>1?"; ":""), $2, $1}')"
    echo "- Known variants (rsID in ID column or VEP Existing_variation), high-confidence set: $known"
    echo
    echo "## Why variants were filtered"
    echo; echo '| filter | failing | failing only this filter |'; echo '|---|---:|---:|'
    tail -n +2 "$Q/filter_breakdown.tsv" | awk -F'\t' '{printf "| %s | %s | %s |\n",$1,$2,$3}'
    echo
    echo "Filter definitions:"; echo
    awk -F'\t' '{printf "- **%s**: %s (`%s`)\n",$1,$3,$2}' "$F/filters.tsv"
    echo
    echo "## Sanity checks (expected ranges for a germline human sample)"
    echo
    echo "| metric | value | expected |"; echo "|---|---|---|"
    echo "| Ts/Tv (high-confidence SNVs) | ${hc_tstv:-NA} | WGS ≈ 2.0–2.1; exome ≈ 2.8–3.3 |"
    echo "| Het/HomAlt | ${hc_hh:-NA} | ≈ 1.4–2.0 (lower with consanguinity, higher with contamination) |"
    echo "| Known (rsID) rate | $known | SNVs > 95–99 % in WGS; indels lower (~85–95 %) |"
    echo "| SNVs | $(sn "$Q/stats.03_high_confidence.txt" 'number of SNPs:') | WGS ≈ 3.5–4.5 M; exome ≈ 20–25 k (coding targets) |"
    echo
    if [[ -s "${OUT}/logs/warnings.txt" ]]; then echo "## Warnings"; echo; cat "${OUT}/logs/warnings.txt"; fi
} >"${OUT}/summary.md"

if [[ $QC_REPORT -eq 1 ]]; then
    qc="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../vcf_qc_pipeline/run_vcf_qc.sh"
    if [[ -x "$qc" ]]; then
        log "running HTML QC report"
        "$qc" -i "${OUT}/final.high_confidence.vcf.gz" -o "${OUT}/qc_report" -r "$REF" -t "$THREADS" \
            --title "High-confidence variants – $SAMPLE" >"${OUT}/logs/qc_report.log" 2>&1 || warn "HTML QC report failed"
    else
        warn "--qc-report: $qc not found"
    fi
fi

if [[ $CLEANUP -eq 1 ]]; then
    rm -f "$F"/with_vaf.vcf.gz* "$F"/regions[12].vcf.gz* "$F"/caller_filters_cleared.vcf.gz* "$V"/soft_filtered.clinvar.vcf.gz*
    rm -rf "$P"/sort_tmp* "$N"/sort_tmp* "$V"/sort_tmp*
    log "intermediate files removed (--cleanup)"
fi

log "done. Final: ${OUT}/final.high_confidence.vcf.gz ($(n_rec "${OUT}/final.high_confidence.vcf.gz") variants); summary: ${OUT}/summary.md"
