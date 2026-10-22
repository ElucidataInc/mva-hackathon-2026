#!/usr/bin/env bash
# =============================================================================
#  Candidate ranking module: two phased VEP-annotated VCFs + annotated SVs
#
#  Steps (folders in --outdir):
#    01_table      per-gene table of each VCF, merged (--combine union | consensus);
#                  ClinVar P/LP calls that failed QC (from -i hc_out) added back, flagged
#    02_coverage   [--bam] per-exon depth (CNV screen) + per-chromosome aneuploidy screen
#    03_sv         SV candidates per gene from the AnnotSV TSV (or plain SV VCFs);
#                  depth support; depth-only exon CNVs
#    04_rank       evidence classes + gene genotypes (SNVs, indels, SVs together) -> tiers
#    05_targets    every small variant (both VCFs) and SV in the target genes
#    summary.md
# =============================================================================
set -Eeuo pipefail
VERSION="3.8"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
cat <<EOF
Usage: $(basename "$0") --vcf1 A.vcf.gz [--vcf2 B.vcf.gz] [--sv annotsv.tsv] -o outdir [options]

Small variants (phased, VEP-annotated, same sample)
      --vcf1 FILE         first VCF (required)
      --vcf2 FILE         second VCF (optional)
      --label1 / --label2 names used in the output [A / B]
      --combine MODE      union = variants in either VCF (flagged by source) | consensus = in both [union]
  -i, --hc-out DIR        vcf_hc_vep.sh output folder of VCF 1 (optional): failed-QC ClinVar calls and
                          filtered variants (if run with --annotate-all)
Structural variants
      --sv FILE           AnnotSV TSV (e.g. merged Manta + cue) or SV VCF; repeatable
      --sv-sample NAME    sample column in the AnnotSV TSV [first sample column]
      --sv-gt MODE        AnnotSV genotype: auto | placeholder | measured [auto: placeholder if every SV
                          has the same GT-only genotype, e.g. an inserted 0/1]
      --gene-gtf FILE     gene models, only for plain SV VCFs (MANE GTF)
      --sv-pop FILE       gnomAD-SV sites VCF, only for plain SV VCFs
Coverage (optional)
      --bam FILE          sample BAM/CRAM (indexed): exon CNV screen, per-chromosome aneuploidy screen, and
                          depth inside every candidate DEL/DUP vs its flanks -> SV genotype (needs samtools
                          >=1.10; mosdepth used for depth if runnable, else samtools)
      --bam-ref FASTA     FASTA of the BAM (needed for CRAM)
      --exon-bed FILE     exons to screen [bundled MVA exons]
      --mosdepth "CMD"    [mosdepth]  e.g. "pixi run -e align --manifest-path M mosdepth"
Ranking
      --genes FILE        gene panel: one gene symbol per line, ranked like the MVA genes
      --exomiser-genes F  Exomiser gene TSV; repeatable (one per VCF/mode): max score per gene
      --alphagenome FILE  AlphaGenome variant table (chrom pos ref alt marker_gene ag_* / sv_* / pv_* columns):
                          splicing / expression / regulatory predictions raise non-coding variant classes
      --max-af F          rare = population AF below F (small variants and SVs) [0.01]
      --cadd-min N        CADD PHRED for non-coding candidates [20]
      --top-discovery N   non-panel biallelic genes kept in targets [20]
  -o, --outdir DIR        output folder [./rank_out]
  -t, --threads N         [4]
EOF
}

VCF1="" VCF2="" L1="A" L2="B" COMBINE="union" HC="" OUT="./rank_out" GENES="" MAX_AF=0.01 CADD_MIN=20 TOPDISC=20
BAM="" BAM_REF="" EXON_BED="" THREADS=4 GTF="" SVPOP="" SVSAMPLE="" SVGT="auto" MOSDEPTH="mosdepth"
SVS=() EXOS=() AG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf1|--vcf) VCF1="$2"; shift ;;
        --vcf2) VCF2="$2"; shift ;;
        --label1) L1="$2"; shift ;;
        --label2) L2="$2"; shift ;;
        --combine) COMBINE="$2"; shift ;;
        -i|--hc-out) HC="$2"; shift ;;
        --sv|--sv-vcf) SVS+=("$2"); shift ;;
        --sv-sample) SVSAMPLE="$2"; shift ;;
        --sv-gt) SVGT="$2"; shift ;;
        --gene-gtf) GTF="$2"; shift ;;
        --sv-pop) SVPOP="$2"; shift ;;
        --bam) BAM="$2"; shift ;;
        --bam-ref) BAM_REF="$2"; shift ;;
        --exon-bed) EXON_BED="$2"; shift ;;
        --mosdepth) MOSDEPTH="$2"; shift ;;
        --genes) GENES="$2"; shift ;;
        --exomiser-genes) EXOS+=("$2"); shift ;;
        --alphagenome) AG="$2"; shift ;;
        --max-af) MAX_AF="$2"; shift ;;
        --cadd-min) CADD_MIN="$2"; shift ;;
        --top-discovery) TOPDISC="$2"; shift ;;
        -o|--outdir) OUT="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
    shift
done
[[ -n "$VCF1" ]] || { usage; exit 1; }
[[ "$COMBINE" == union || "$COMBINE" == consensus ]] || { echo "--combine must be union or consensus" >&2; exit 1; }
mkdir -p "$OUT"/{logs,01_table,02_coverage,03_sv,04_rank,05_targets}
OUT="$(cd "$OUT" && pwd)"
: >"$OUT/logs/warnings.txt"

log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$OUT/logs/ranking.log" >&2; }
warn() { log "WARNING: $*"; echo "- $*" >>"$OUT/logs/warnings.txt"; }
die()  { log "ERROR: $*"; exit 1; }
trap 'rc=$?; log "ERROR: unexpected failure (exit $rc) at line $LINENO: $BASH_COMMAND"; exit $rc' ERR
have() { command -v "$1" >/dev/null 2>&1; }
runs() { bash -c "$1 --version" >/dev/null 2>&1 || bash -c "$1 --help" >/dev/null 2>&1; }

have bcftools || die "bcftools not found (source ~/vcf_env.sh)"
have python3 || die "python3 not found"
if [[ -n "$AG" ]]; then
    [[ -s "$AG" ]] || die "AlphaGenome table not found: $AG"
    AG="$(cd "$(dirname "$AG")" && pwd)/$(basename "$AG")"
    for c in chrom pos ref alt marker_gene; do
        head -1 "$AG" | tr '\t' '\n' | grep -cx "$c" >/dev/null || die "AlphaGenome table has no '$c' column: $AG"
    done
fi
[[ -z "$HC" ]] || HC="$(cd "$HC" && pwd)"
ALLVCF="${HC:+$HC/04_vep/annotated.vcf.gz}"
FAILED="${HC:+$HC/05_tables/clinvar_plp_failed_qc.tsv}"
EXON_BED="${EXON_BED:-$HERE/mva_exons_GRCh38.bed}"
for f in VCF1 VCF2 GENES BAM BAM_REF GTF SVPOP EXON_BED; do v="${!f}"; [[ -z "$v" || -f "$v" ]] || die "$f not found: $v"; done
for x in ${SVS[@]+"${SVS[@]}"} ${EXOS[@]+"${EXOS[@]}"}; do [[ -f "${x#*=}" ]] || die "file not found: ${x#*=}"; done
for v in "$VCF1" ${VCF2:+"$VCF2"}; do
    bcftools view -h "$v" | grep -c '^##INFO=<ID=CSQ,' >/dev/null || die "$v has no VEP annotation (CSQ)"
    [[ -f "$v.tbi" || -f "$v.csi" ]] || bcftools index -t "$v" 2>/dev/null || die "cannot index $v (sort it first)"
    nm=$( { bcftools view -H "$v" 2>/dev/null || true; } | awk -F'\t' 'NR>200000{exit} $5 ~ /,/{n++} END{print n+0}')
    [[ "$nm" -eq 0 ]] || warn "$(basename "$v") has multi-allelic records ($nm in the first 200k): normalise (bcftools norm -m -any) so both VCFs match"
done
style_chr=0; bcftools view -h "$VCF1" | grep -c '^##contig=<ID=chr' >/dev/null && style_chr=1

phase_stat() {   # phase_stat <vcf>: "phased_het/het"
    { bcftools query -f '[%GT]\n' "$1" 2>/dev/null || true; } |
        awk 'NR>200000{exit} { n=split($1,a,/[\/|]/); if (n==2 && a[1]!=a[2] && a[1]!="." && a[2]!=".") { h++; if ($1 ~ /\|/) p++ } }
             END{printf "%d/%d", p+0, h+0}'
}
log "ranking module v$VERSION -> $OUT"
log "  VCF $L1: $VCF1 (phased hets $(phase_stat "$VCF1"))"
[[ -n "$VCF2" ]] && log "  VCF $L2: $VCF2 (phased hets $(phase_stat "$VCF2")); combine = $COMBINE"

# ----------------------------------------------------------------------------
build_table() {   # build_table <vcf> <out.tsv>: one row per variant per gene (VEP PICK), pre-filtered while streaming
    local vcf="$1" out="$2" hdr fq="" fh="" iq="" ih="" fields="" x cq ch sel=()
    hdr=$(bcftools view -h "$vcf")
    for x in GT GQ DP VAF AD PS; do
        grep -c "^##FORMAT=<ID=$x," <<<"$hdr" >/dev/null && { fq="${fq}\t%${x}"; fh="${fh}\t${x}"; }
    done
    for x in CLNSIG CLNREVSTAT CLNDN ALLELEID; do
        grep -c "^##INFO=<ID=$x," <<<"$hdr" >/dev/null && { iq="${iq}\t%INFO/${x}"; ih="${ih}\tClinVar_${x}"; }
    done
    bcftools +split-vep -l "$vcf" | cut -f2 >"$out.csq_fields"
    for x in Consequence IMPACT SYMBOL Gene Feature BIOTYPE MANE_SELECT CANONICAL EXON INTRON HGVSc HGVSp \
             Protein_position Amino_acids Existing_variation VARIANT_CLASS SIFT PolyPhen gnomADe_AF gnomADg_AF \
             AF MAX_AF MAX_AF_POPS CLIN_SIG CADD_PHRED CADD_RAW SpliceAI_pred_DS_AG SpliceAI_pred_DS_AL \
             SpliceAI_pred_DS_DG SpliceAI_pred_DS_DL REVEL am_class am_pathogenicity PICK; do
        grep -qx "$x" "$out.csq_fields" && fields="$fields $x"
    done
    # CSQ AF (1000 Genomes) is already part of MAX_AF; asking for %AF also triggers bcftools' "ambiguous key" note
    [[ " $fields " == *" MAX_AF "* ]] && fields="${fields/ AF / }"
    cq=$(for x in $fields; do printf '\\t%%%s' "$x"; done)
    ch=$(for x in $fields; do printf '\t%s' "$x"; done)
    grep -qx PICK "$out.csq_fields" && sel=(-i 'PICK="1"')
    { printf "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER${ih}${fh}${ch}\n"
      bcftools +split-vep "$vcf" -d ${sel[@]+"${sel[@]}"} \
          -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER${iq}${fq:+[${fq}]}${cq}\n"
    } | python3 "$HERE/prefilter_table.py" --max-af "$MAX_AF" --cadd-min "$CADD_MIN" --log "$out.prefilter.log" \
        ${AG:+--keep-sites "$AG"} 2>/dev/null >"$out"
}

add_failed_clinvar() {   # ClinVar P/LP calls that FAILED a filter (needs -i hc_out from --annotate-all)
    local tab="$1" sub="$OUT/01_table/failed_clinvar.vcf.gz" n
    [[ -n "$ALLVCF" && -f "$ALLVCF" ]] || return 0
    bcftools view -h "$ALLVCF" | grep -c '^##INFO=<ID=CLNSIG,' >/dev/null || return 0
    { bcftools view -e 'INFO/CLNSIG~"onflicting"' -Ou "$ALLVCF" |
      bcftools view -i 'FILTER!="PASS" && INFO/CLNSIG~"athogenic"' -Oz -o "$sub"; } 2>/dev/null || return 0
    n=$(bcftools view -H "$sub" | wc -l | tr -d ' ')
    [[ "$n" -gt 0 ]] || return 0
    build_table "$sub" "$OUT/01_table/failed_clinvar.tsv"
    # the hc_out VCF may come from another VEP run than VCF 1: align its columns to VCF 1's table by name
    local added missing
    local res; res=$(python3 "$HERE/append_aligned.py" "$tab" "$OUT/01_table/failed_clinvar.tsv")
    added="${res%% *}" missing="${res#* }"
    log "  ClinVar P/LP calls that failed QC: $n found, $added gene row(s) added to VCF $L1 table (flagged)$( [[ $added -lt $n ]] && echo "; the rest are outside genes or not rare" || true)"
    if [[ "$missing" != "-" ]]; then
        # PS (phase) is expected to be missing: the failed calls come from the unphased hc_out VCF
        if [[ -z "$(tr ',' '\n' <<<"$missing" | grep -vxE 'PS|GQ|DP|VAF|AD' || true)" ]]; then
            log "  (failed-QC rows have no $missing: they come from the unphased hc_out VCF, so they are unphased)"
        else
            warn "failed-QC ClinVar rows come from a different VEP run than VCF $L1: columns missing there were left empty ($missing)"
        fi
    fi
}

bam_contigs() { samtools view -H "$BAM" | awk -F'\t' '$1=="@SQ"{for(i=2;i<=NF;i++) if ($i ~ /^SN:/) print substr($i,4)}'; }
to_style() {
    awk -F'\t' -v OFS='\t' -v sc="$1" '{c=$1; if (sc && c !~ /^chr/) {c="chr" c; if (c=="chrMT") c="chrM"}
        if (!sc && c ~ /^chr/) {sub(/^chr/,"",c); if (c=="M") c="MT"} $1=c; print}'
}

# =============================================================================
log "01 table"
build_table "$VCF1" "$OUT/01_table/table_$L1.tsv"
add_failed_clinvar "$OUT/01_table/table_$L1.tsv"
log "  VCF $L1: $(sed 's/^\[prefilter\] //' "$OUT/01_table/table_$L1.tsv.prefilter.log")"
if [[ -n "$VCF2" ]]; then
    build_table "$VCF2" "$OUT/01_table/table_$L2.tsv"
    log "  VCF $L2: $(sed 's/^\[prefilter\] //' "$OUT/01_table/table_$L2.tsv.prefilter.log")"
    python3 "$HERE/merge_tables.py" --t1 "$OUT/01_table/table_$L1.tsv" --t2 "$OUT/01_table/table_$L2.tsv" \
        --label1 "$L1" --label2 "$L2" --combine "$COMBINE" --style-chr "$style_chr" \
        --out "$OUT/01_table/annotated.tsv" 2>"$OUT/01_table/merge.log"
    sed 's/^/  /' "$OUT/01_table/merge.log" >&2
else
    cp "$OUT/01_table/table_$L1.tsv" "$OUT/01_table/annotated.tsv"
fi

# =============================================================================
C="$OUT/02_coverage"
BAM_OK=0 COV_TOOL=""
if [[ -n "$BAM" ]]; then
    if ! have samtools; then
        warn "--bam given but samtools is not on PATH: coverage checks and SV depth genotypes skipped (activate the env with samtools)"
    elif ! samtools --version >/dev/null 2>&1; then
        warn "--bam given but '$(command -v samtools)' is samtools 0.x (too old; the VEP env often ships one): put samtools >=1.10 first on PATH; coverage checks skipped"
    elif [[ ! -f "$BAM.bai" && ! -f "${BAM%.*}.bai" && ! -f "$BAM.csi" && ! -f "$BAM.crai" && ! -f "${BAM%.*}.crai" ]]; then
        warn "--bam given but no index found next to it (.bai/.csi/.crai): run 'samtools index $BAM'; coverage checks skipped"
    elif [[ "$BAM" == *.cram && -z "$BAM_REF" ]]; then
        warn "--bam is a CRAM: also give --bam-ref (the FASTA it was written with); coverage checks skipped"
    elif ! samtools idxstats "$BAM" >"$C/idxstats.txt" 2>"$C/idxstats.log"; then
        warn "--bam given but samtools idxstats failed (index .bai/.crai missing? CRAM without --bam-ref?): see $C/idxstats.log; coverage checks skipped"
    else
        BAM_OK=1
        if runs "$MOSDEPTH"; then COV_TOOL=mosdepth
        else
            COV_TOOL=samtools
            warn "--bam given but mosdepth is not runnable ('$MOSDEPTH'): using samtools bedcov instead (same results, slower). Install mosdepth or pass --mosdepth \"pixi run ... mosdepth\" to speed this up"
        fi
    fi
fi
if [[ "$BAM_OK" == 1 ]]; then
    log "02 coverage: exon CNV screen + per-chromosome screen ($COV_TOOL); BAM: $BAM"
    bam_contigs >"$C/bam_contigs.txt"
    bchr=0; grep -q '^chr' "$C/bam_contigs.txt" && bchr=1
    grep -v '^#' "$EXON_BED" | cut -f1-4 | to_style $bchr |
        awk -F'\t' -v list="$C/bam_contigs.txt" 'BEGIN{while((getline l<list)>0) ok[l]=1} ($1 in ok)' >"$C/exons.bam.bed"
    [[ -s "$C/exons.bam.bed" ]] || warn "no exon of $EXON_BED lies on a contig of the BAM; exon CNV screen skipped"
    gm="" ref_label=""
    if [[ "$COV_TOOL" == mosdepth ]]; then
        fa=""; [[ "$BAM" == *.cram ]] && fa="--fasta '$BAM_REF'"
        bash -c "$MOSDEPTH -t $THREADS -n --fast-mode -Q 1 $fa --by 1000000 '$C/genome' '$BAM'" >"$C/mosdepth_genome.log" 2>&1 || warn "mosdepth (genome) failed: see $C/mosdepth_genome.log"
        [[ -f "$C/genome.regions.bed.gz" ]] && gzip -cd "$C/genome.regions.bed.gz" >"$C/windows_1Mb.tsv"
        if [[ -s "$C/exons.bam.bed" ]]; then
            rm -f "$C"/exons.regions.bed.gz*
            bash -c "$MOSDEPTH -t $THREADS -n --fast-mode -Q 1 $fa --by '$C/exons.bam.bed' '$C/exons' '$BAM'" >"$C/mosdepth_exons.log" 2>&1 ||
                { warn "mosdepth (exons) failed: see $C/mosdepth_exons.log"; rm -f "$C"/exons.regions.bed.gz*; }
            [[ -f "$C/exons.regions.bed.gz" ]] && gzip -cd "$C/exons.regions.bed.gz" | cut -f1-5 >"$C/exon_means.tsv"
        fi

    else
        if [[ -s "$C/exons.bam.bed" ]]; then
            # bedcov: summed per-base depth (MAPQ >= 1) -> mean depth per exon
            samtools bedcov -Q 1 ${BAM_REF:+--reference "$BAM_REF"} "$C/exons.bam.bed" "$BAM" 2>"$C/bedcov_exons.log" |
                awk -F'\t' -v OFS='\t' '{l=$3-$2; print $1,$2,$3,$4,(l>0? $5/l : 0)}' >"$C/exon_means.tsv" ||
                warn "samtools bedcov (exons) failed: see $C/bedcov_exons.log"
        fi
        # 1-Mb windows over the main chromosomes (one pass over the BAM, like mosdepth --by 1000000)
        awk -F'\t' -v OFS='\t' '$1 ~ /^(chr)?([0-9]+|X|Y)$/ {for (s=0; s<$2; s+=1000000) print $1, s, (s+1000000<$2? s+1000000 : $2)}' \
            "$C/idxstats.txt" >"$C/windows_1Mb.bed"
        samtools bedcov -Q 1 ${BAM_REF:+--reference "$BAM_REF"} "$C/windows_1Mb.bed" "$BAM" 2>"$C/bedcov_windows.log" |
            awk -F'\t' -v OFS='\t' '{l=$3-$2; print $1,$2,$3,(l>0? $4/l : 0)}' >"$C/windows_1Mb.tsv" ||
            warn "samtools bedcov (1-Mb windows) failed: see $C/bedcov_windows.log"
    fi
    # exon reference = median depth of the screened exons: exons share GC content / mappability, so this cancels
    # most of the exon-vs-genome bias (the genome mean also counts N gaps and is too low)
    gm=$(awk -F'\t' '$5>0{print $5}' "$C/exon_means.tsv" 2>/dev/null | sort -n | awk '{a[NR]=$1} END{if(NR) print a[int((NR+1)/2)]}')
    ref_label="median depth of the screened exons"
    # chromosome screen: median of 1-Mb windows per chromosome, gap windows ignored
    sig="median_1Mb_window_depth" gwin=""
    if [[ -s "$C/windows_1Mb.tsv" ]]; then
        gwin=$(python3 "$HERE/window_medians.py" <"$C/windows_1Mb.tsv" 2>&1 >"$C/chrom_means.txt" | tail -1) || { warn "chromosome window medians failed"; rm -f "$C/chrom_means.txt"; }
    fi
    if [[ -s "$C/exon_means.tsv" ]]; then
        { printf 'chrom\tstart\tend\texon\tmean_depth\tratio_to_reference\tflag\n'
          awk -F'\t' -v OFS='\t' -v gm="${gm:-0}" '{
              r=(gm>0)? $5/gm : -1; f=""
              if ($5 < 10) f="low coverage (<10x)"
              else if (r>=0 && r<0.65) f="possible deletion"
              else if (r>1.35) f="possible duplication"
              printf "%s\t%s\t%s\t%s\t%.1f\t%s\t%s\n", $1,$2,$3,$4,$5, (r>=0? sprintf("%.2f",r):"NA"), f}' "$C/exon_means.tsv"
        } >"$C/exon_depth.tsv"
        log "  exons: $(($(wc -l <"$C/exon_depth.tsv") - 1)) screened (reference: $ref_label ${gm:-?}x; genome median 1-Mb window ${gwin:-?}x); $(awk -F'\t' 'NR>1 && $7 ~ /deletion|duplication/' "$C/exon_depth.tsv" | wc -l | tr -d ' ') possible CNV, $(awk -F'\t' 'NR>1 && $7 ~ /low/' "$C/exon_depth.tsv" | wc -l | tr -d ' ') low coverage"
    fi
    if [[ -s "$C/chrom_means.txt" ]]; then
        med=$(awk '$1 ~ /^(chr)?[0-9]+$/ {print $2}' "$C/chrom_means.txt" | sort -n | awk '{a[NR]=$1} END{print (NR? a[int((NR+1)/2)] : 0)}')
        if awk -v m="$med" 'BEGIN{exit !(m+0 > 0)}'; then
            { printf 'chrom\t%s\tratio_to_autosomal_median\twindows_used\tflag\n' "$sig"
              awk -F'\t' -v m="$med" '{r=$2/m; f="";
                   if ($1 ~ /^(chr)?[0-9]+$/ && (r>1.08 || r<0.92)) f="possible mosaic gain/loss"
                   if ($1 ~ /X|Y/) f="sex chromosome: compare with expected sex"
                   printf "%s\t%.2f\t%.3f\t%s\t%s\n", $1, $2, r, $3, f}' "$C/chrom_means.txt"
            } >"$C/karyotype_screen.tsv"
            log "  chromosomes: $(awk -F'\t' 'NR>1 && $5 ~ /mosaic/' "$C/karyotype_screen.tsv" | wc -l | tr -d ' ') autosome(s) outside 0.92-1.08x of the median ($sig)"
        else
            warn "autosomal median is 0; aneuploidy screen skipped"
        fi
    fi
elif [[ -z "$BAM" ]]; then
    log "02 coverage: skipped (no --bam)"
else
    log "02 coverage: skipped (BAM not usable; see warnings)"
fi

# =============================================================================
S="$OUT/03_sv"
dep=(); [[ -f "$C/exon_depth.tsv" ]] && dep=(--exon-depth "$C/exon_depth.tsv")
svbam=(); [[ "$BAM_OK" == 1 ]] && svbam=(--bam "$BAM" ${BAM_REF:+--bam-ref "$BAM_REF"})
if [[ ${#SVS[@]} -gt 0 ]]; then
    log "03 SV: ${#SVS[@]} file(s)"
    svargs=()
    for s in "${SVS[@]}"; do svargs+=(--sv "$s"); done
    python3 "$HERE/sv_annotate.py" "${svargs[@]}" --out "$S/sv_candidates.tsv" --max-af "$MAX_AF" --style-chr "$style_chr" --sv-gt "$SVGT" \
        ${SVSAMPLE:+--sample "$SVSAMPLE"} ${GTF:+--gtf "$GTF"} ${SVPOP:+--pop "$SVPOP"} ${dep[@]+"${dep[@]}"} ${svbam[@]+"${svbam[@]}"} \
        2>&1 | tee "$S/sv_annotate.log" | sed 's/^/  /' >&2
elif [[ ${#dep[@]} -gt 0 ]]; then
    python3 "$HERE/sv_annotate.py" --out "$S/sv_candidates.tsv" --style-chr "$style_chr" "${dep[@]}" 2>"$S/sv_annotate.log" || true
    log "03 SV: no SV file; depth-only exon CNV candidates: $(( $(wc -l <"$S/sv_candidates.tsv") - 1 ))"
else
    log "03 SV: skipped (no --sv)"
fi

# =============================================================================
log "04 rank"
exo=(); for e in ${EXOS[@]+"${EXOS[@]}"}; do exo+=(--exomiser-genes "$e"); done
python3 "$HERE/rank_variants.py" --table "$OUT/01_table/annotated.tsv" --out "$OUT/04_rank" --max-af "$MAX_AF" \
    --cadd-min "$CADD_MIN" --top-discovery "$TOPDISC" ${GENES:+--genes "$GENES"} ${exo[@]+"${exo[@]}"} \
    $( [[ -s "$S/sv_candidates.tsv" ]] && echo --sv-table "$S/sv_candidates.tsv" || true ) \
    $( [[ -n "$FAILED" && -f "$FAILED" ]] && echo --failed-qc "$FAILED" || true ) \
    ${AG:+--alphagenome "$AG"} >"$OUT/04_rank/rank.log"
sed 's/^/  /' "$OUT/04_rank/rank.log" >&2
cp "$OUT/04_rank"/{ranked_genes.tsv,ranked_variants.tsv,targets.tsv,targets.bed,top_mva_variants.tsv} "$OUT/"

# =============================================================================
T="$OUT/05_targets"
cut -f2 "$OUT/targets.tsv" | tail -n +2 | sort -u >"$T/target_genes.txt"
nt=$(wc -l <"$T/target_genes.txt" | tr -d ' ')
log "05 targets: $nt gene(s)"
if [[ "$nt" -gt 0 ]]; then
    printf 'SOURCE\tCHROM\tPOS\tREF\tALT\tFILTER\tGT\tPS\tDP\tVAF\tSYMBOL\tConsequence\tIMPACT\tHGVSc\tHGVSp\tMAX_AF\tExisting_variation\n' >"$T/target_gene_variants.tsv"
    srcs=("$L1|$VCF1"); [[ -n "$VCF2" ]] && srcs+=("$L2|$VCF2")
    [[ -n "$ALLVCF" && -f "$ALLVCF" ]] && srcs+=("$L1(incl.filtered)|$ALLVCF")
    for sv in "${srcs[@]}"; do
        lab="${sv%%|*}" f="${sv#*|}"
        hdr=$(bcftools view -h "$f"); fl=$(bcftools +split-vep -l "$f" | cut -f2)
        fmt=""; for x in GT PS DP VAF; do grep -c "^##FORMAT=<ID=$x," <<<"$hdr" >/dev/null && fmt="$fmt\t%$x" || fmt="$fmt\t."; done
        csq=""; for x in SYMBOL Consequence IMPACT HGVSc HGVSp MAX_AF Existing_variation; do grep -qx "$x" <<<"$fl" && csq="$csq\t%$x" || csq="$csq\t."; done
        bcftools +split-vep "$f" -d -f "$lab\t%CHROM\t%POS\t%REF\t%ALT\t%FILTER[$fmt]$csq\n" 2>/dev/null |
            awk -F'\t' 'NR==FNR{g[$1]=1; next} ($11 in g) && !seen[$2 FS $3 FS $4 FS $5 FS $11]++' "$T/target_genes.txt" - >>"$T/target_gene_variants.tsv" || true
    done
    log "  $(($(wc -l <"$T/target_gene_variants.tsv") - 1)) small-variant rows in target genes (from ${#srcs[@]} source(s))"
    if [[ -s "$S/sv_candidates.tsv" ]]; then
        awk -F'\t' 'NR==FNR{g[$1]=1; next} FNR==1 || ($1 in g)' "$T/target_genes.txt" "$S/sv_candidates.tsv" >"$T/target_gene_svs.tsv"
        log "  $(($(wc -l <"$T/target_gene_svs.tsv") - 1)) SV/CNV rows in target genes"
    fi
fi

# =============================================================================
tier_meaning() {
    case "$1" in
        1) echo "MVA/panel gene, biallelic, ≥1 ClinVar P/LP" ;;
        2) echo "MVA/panel gene, biallelic, both alleles strong (SNV/indel/SV)" ;;
        3) echo "MVA/panel gene biallelic weaker, or emerging gene" ;;
        4) echo "MVA/panel gene, one strong allele: find the second hit" ;;
        5) echo "other gene, biallelic strong (discovery)" ;;
        *) echo "other" ;;
    esac
}
{
    echo "# Candidate ranking – summary"
    echo
    echo "- VCF $L1: \`$VCF1\`$( [[ -n "$VCF2" ]] && echo "  ·  VCF $L2: \`$VCF2\`  ·  combine: **$COMBINE**" || true)"
    [[ -n "$VCF2" ]] && echo "- Merge: $(sed 's/^\[merge\] //' "$OUT/01_table/merge.log" 2>/dev/null | tail -1 || true)"
    echo "- SV input: $( [[ ${#SVS[@]} -gt 0 ]] && printf '%s ' "${SVS[@]##*/}" || echo none)$( [[ -s "$S/sv_candidates.tsv" ]] && echo "· $(( $(wc -l <"$S/sv_candidates.tsv") - 1 )) SV/CNV-gene rows" || true)"
    [[ -s "$S/mva_gene_svs_all.tsv" ]] && echo "- SVs near MVA genes (input): $(cut -f1 "$S/mva_gene_svs_all.tsv" | tail -n +2 | sort -u | wc -l | tr -d ' ') SV(s); $(awk -F'\t' 'NR>1 && $10=="yes"' "$S/mva_gene_svs_all.tsv" | wc -l | tr -d ' ') SV-gene pair(s) used in the ranking; fate of each: \`03_sv/mva_gene_svs_all.tsv\`"
    echo "- AlphaGenome: $( [[ -n "$AG" ]] && echo "\`${AG##*/}\` · $(grep -h '^\[rank\] AlphaGenome' "$OUT/04_rank/rank.log" | sed 's/^\[rank\] AlphaGenome: //' | paste -sd'|' - | sed 's/|/ · /g')" || echo none)"
    echo "- Exomiser: ${#EXOS[@]} file(s)  ·  coverage: $( [[ -f "$C/exon_depth.tsv" ]] && echo done || echo "not run")  ·  module v$VERSION  ·  $(date)"
    echo
    echo "## Genes per tier"; echo
    echo '| tier | meaning | genes |'; echo '|---|---|---:|'
    for t in 1 2 3 4 5 6; do
        echo "| $t | $(tier_meaning $t) | $(awk -F'\t' -v t=$t 'NR>1 && $1==t' "$OUT/ranked_genes.tsv" | wc -l | tr -d ' ') |"
    done
    echo
    echo "## Targets (\`targets.tsv\`)"; echo
    echo '| tier | gene | genotype | phase | Exomiser | alleles | note |'; echo '|---|---|---|---|---|---|---|'
    awk -F'\t' 'NR>1{v=$9; gsub(/\|/, "\\|", v); printf "| %s | %s | %s | %s | %s | %s | %s |\n", $1, $2, $4, $5, $7, v, $10}' "$OUT/targets.tsv"
    echo
    if [[ -f "$C/exon_depth.tsv" ]] && awk -F'\t' 'NR>1 && $7 ~ /deletion|duplication/' "$C/exon_depth.tsv" | grep -c . >/dev/null; then
        echo "## Exons with CNV-like depth"; echo
        echo '| exon | depth | ratio | flag |'; echo '|---|---:|---:|---|'
        awk -F'\t' 'NR>1 && $7 ~ /deletion|duplication/ {printf "| %s | %s | %s | %s |\n", $4, $5, $6, $7}' "$C/exon_depth.tsv" | head -40; echo
    fi
    if [[ -f "$C/karyotype_screen.tsv" ]] && awk -F'\t' 'NR>1 && $5 ~ /mosaic/' "$C/karyotype_screen.tsv" | grep -c . >/dev/null; then
        echo "## Aneuploidy screen (bulk depth; normal does not exclude MVA)"; echo
        awk -F'\t' 'NR>1 && $5 ~ /mosaic/ {printf "- %s: ratio %s\n", $1, $3}' "$C/karyotype_screen.tsv"; echo
    fi
    echo "## Next steps"; echo
    echo "1. Tiers 1–2: confirm in IGV (\`targets.bed\`), Sanger/MLPA as appropriate; segregation in parents; ACMG-classify."
    echo "2. Alleles found in only one VCF, 'unphased' pairs and 'SV + small variant' pairs: confirm before use."
    echo "3. Tier 4: look for a second hit in \`05_targets/\` (non-coding, filtered, SV/CNV), RNA-seq."
    echo
    echo "Scores order candidates for review; they are not ACMG classifications."
    [[ -s "$OUT/logs/warnings.txt" ]] && { echo; echo "## Warnings"; echo; cat "$OUT/logs/warnings.txt"; }
} >"$OUT/summary.md"
log "done: $OUT/summary.md  ·  targets: $OUT/targets.tsv ($(($(wc -l <"$OUT/targets.tsv") - 1)) genes)"
