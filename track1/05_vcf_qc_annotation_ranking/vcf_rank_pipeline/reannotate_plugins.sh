#!/usr/bin/env bash
# =============================================================================
#  Re-annotate an already phased, VEP-annotated VCF with extra VEP plugins
#  (REVEL, SpliceAI, CADD, AlphaMissense) for run_ranking.sh.
#
#  Keeps genotypes, phase sets (PS) and ClinVar INFO tags; only the CSQ field is rebuilt,
#  with the same VEP options as vcf_hc_vep.sh (--everything --flag_pick_allele_gene).
#  By default only variants that the ranking can use are re-annotated (population AF < --max-af,
#  no AF, or ClinVar pathogenic in VEP's CLIN_SIG): about 20% of a genome, so VEP runs ~5x faster.
#
#  Usage:
#    reannotate_plugins.sh -i sample.phased.vcf.gz -o sample.phased.plugins.vcf.gz -r GRCh38.fa \
#        --vep-cache ~/.vep --plugins-dir <dir with REVEL.pm SpliceAI.pm CADD.pm AlphaMissense.pm> \
#        --revel new_tabbed_revel_grch38.tsv.gz \
#        --spliceai-snv spliceai_scores.masked.snv.hg38.vcf.gz --spliceai-indel spliceai_scores.masked.indel.hg38.vcf.gz \
#        --cadd-snv whole_genome_SNVs.tsv.gz --cadd-indel gnomad.genomes.r4.0.indel.tsv.gz \
#        --alphamissense AlphaMissense_hg38.tsv.gz -t 8
#  Any plugin file can be left out; the run reports which score fields the output contains.
# =============================================================================
set -Eeuo pipefail

usage() { sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
IN="" OUT="" REF="" CACHE="" GTF="" PLUGDIR="" REVEL="" SAI_SNV="" SAI_INDEL="" CADD_SNV="" CADD_INDEL="" AM=""
THREADS=4 MAX_AF=0.01 ALL=0 ASSEMBLY=GRCh38
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) IN="$2"; shift ;;
        -o) OUT="$2"; shift ;;
        -r|--ref) REF="$2"; shift ;;
        --vep-cache) CACHE="$2"; shift ;;
        --vep-gtf) GTF="$2"; shift ;;
        --plugins-dir) PLUGDIR="$2"; shift ;;
        --revel) REVEL="$2"; shift ;;
        --spliceai-snv) SAI_SNV="$2"; shift ;;
        --spliceai-indel) SAI_INDEL="$2"; shift ;;
        --cadd-snv) CADD_SNV="$2"; shift ;;
        --cadd-indel) CADD_INDEL="$2"; shift ;;
        --alphamissense) AM="$2"; shift ;;
        --max-af) MAX_AF="$2"; shift ;;
        --all) ALL=1 ;;
        --assembly) ASSEMBLY="$2"; shift ;;
        -t) THREADS="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "unknown option: $1" >&2; usage 1 ;;
    esac
    shift
done
log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }
[[ -n "$IN" && -n "$OUT" && -n "$REF" ]] || usage 1
[[ -n "$CACHE" || -n "$GTF" ]] || die "give --vep-cache (or --vep-gtf)"
for t in bcftools vep filter_vep; do command -v "$t" >/dev/null || die "$t not on PATH (activate the VEP env)"; done
for f in "$IN" "$REF" ${REVEL:+"$REVEL"} ${SAI_SNV:+"$SAI_SNV"} ${SAI_INDEL:+"$SAI_INDEL"} ${CADD_SNV:+"$CADD_SNV"} \
         ${CADD_INDEL:+"$CADD_INDEL"} ${AM:+"$AM"}; do
    [[ -s "$f" ]] || die "file not found: $f"
done
for f in ${REVEL:+"$REVEL"} ${SAI_SNV:+"$SAI_SNV"} ${SAI_INDEL:+"$SAI_INDEL"} ${CADD_SNV:+"$CADD_SNV"} ${CADD_INDEL:+"$CADD_INDEL"} ${AM:+"$AM"}; do
    [[ -f "$f.tbi" || -f "$f.csi" ]] || die "plugin file has no tabix index (.tbi/.csi): $f"
done
if [[ -n "$PLUGDIR" ]]; then
    for p in ${REVEL:+REVEL} ${SAI_SNV:+SpliceAI} ${CADD_SNV:+CADD} ${AM:+AlphaMissense}; do
        [[ -f "$PLUGDIR/$p.pm" ]] || die "$p.pm not found in --plugins-dir $PLUGDIR"
    done
fi
[[ -n "$REVEL$SAI_SNV$CADD_SNV$AM" ]] || log "WARNING: no plugin file given: only the standard VEP annotation will be rebuilt"

tmp="${OUT%.vcf.gz}.tmp"; mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT

# 1. variants the ranking can use (rare or ClinVar pathogenic) - selected on the existing CSQ
if [[ $ALL -eq 0 ]] && bcftools +split-vep -l "$IN" 2>/dev/null | cut -f2 | grep -cx MAX_AF >/dev/null; then
    log "1/3 selecting variants with MAX_AF < $MAX_AF, no MAX_AF, or ClinVar pathogenic (CLIN_SIG)"
    bcftools view "$IN" |
        filter_vep --format vcf --filter "MAX_AF < $MAX_AF or not MAX_AF or CLIN_SIG match pathogenic" --force_overwrite 2>"$tmp/filter_vep.log" |
        bcftools annotate -x INFO/CSQ -Oz -o "$tmp/in.vcf.gz"
else
    log "1/3 using all variants$( [[ $ALL -eq 1 ]] && echo " (--all)" || echo " (input has no MAX_AF to select on)")"
    bcftools annotate -x INFO/CSQ -Oz -o "$tmp/in.vcf.gz" "$IN"
fi
n_in=$(bcftools view -H "$IN" | wc -l | tr -d ' ')
n_sel=$(bcftools view -H "$tmp/in.vcf.gz" | wc -l | tr -d ' ')
log "    $n_sel of $n_in records will be re-annotated"
[[ "$n_sel" -gt 0 ]] || die "no records selected (see $tmp/filter_vep.log)"

# 2. VEP with the same options as vcf_hc_vep.sh + plugins
args=(--input_file "$tmp/in.vcf.gz" --output_file "$OUT" --vcf --compress_output bgzip --force_overwrite
      --species homo_sapiens --assembly "$ASSEMBLY" --fasta "$REF" --fork "$THREADS" --buffer_size 5000
      --stats_file "${OUT%.vcf.gz}.vep_summary.html" --warning_file "${OUT%.vcf.gz}.vep_warnings.txt"
      --flag_pick_allele_gene)
if [[ -n "$CACHE" ]]; then args+=(--offline --cache --dir_cache "$CACHE" --everything)
else args+=(--gtf "$GTF" --symbol --biotype --canonical --numbers --hgvs --variant_class); fi
[[ -n "$PLUGDIR" ]] && args+=(--dir_plugins "$PLUGDIR")
[[ -n "$REVEL" ]] && args+=(--plugin "REVEL,file=$REVEL")
[[ -n "$SAI_SNV" ]] && args+=(--plugin "SpliceAI,snv=$SAI_SNV${SAI_INDEL:+,indel=$SAI_INDEL}")
[[ -n "$CADD_SNV" ]] && args+=(--plugin "CADD,snv=$CADD_SNV${CADD_INDEL:+,indels=$CADD_INDEL}")
[[ -n "$AM" ]] && args+=(--plugin "AlphaMissense,file=$AM")
printf '%q ' vep "${args[@]}" >"${OUT%.vcf.gz}.vep_command.txt"; echo >>"${OUT%.vcf.gz}.vep_command.txt"
log "2/3 VEP ($THREADS forks): $(sed 's/.*--flag_pick_allele_gene//' "${OUT%.vcf.gz}.vep_command.txt" | cut -c1-200)"
vep "${args[@]}" >"${OUT%.vcf.gz}.vep.log" 2>&1 || { tail -20 "${OUT%.vcf.gz}.vep.log" >&2; die "VEP failed (see ${OUT%.vcf.gz}.vep.log)"; }
bcftools index -t -f "$OUT"
# VEP only warns when a plugin cannot start (wrong file, wrong assembly, missing Perl module): show it
if grep -q 'Failed to instantiate plugin' "${OUT%.vcf.gz}.vep.log"; then
    grep -A1 'Failed to instantiate plugin' "${OUT%.vcf.gz}.vep.log" | grep -v '^--$' | sort -u | sed 's/^/    PLUGIN ERROR: /' >&2
fi

# 3. check: score fields in CSQ and how many records got a value
log "3/3 score fields in the new CSQ:"
fields=$(bcftools +split-vep -l "$OUT" | cut -f2)
for x in REVEL SpliceAI_pred_DS_AG CADD_PHRED am_class; do
    if grep -qx "$x" <<<"$fields"; then
        n=$(bcftools +split-vep "$OUT" -f "%$x\n" -d 2>/dev/null | awk '$1!="." && $1!=""' | wc -l | tr -d ' ')
        log "    [x] $x  ($n transcript rows with a value)"
    else
        log "    [ ] $x  (not added)"
    fi
done
n_ps=$(bcftools view -h "$OUT" | grep -c '^##FORMAT=<ID=PS,' || true)
log "    phase sets (FORMAT/PS) kept: $( [[ "$n_ps" -gt 0 ]] && echo yes || echo NO)"
log "done: $OUT  ->  use it as --vcf1/--vcf2 in run_ranking.sh"
