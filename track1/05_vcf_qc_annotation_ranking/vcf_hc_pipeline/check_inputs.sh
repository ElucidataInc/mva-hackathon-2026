#!/usr/bin/env bash
# Pre-flight check for vcf_hc_vep.sh: verifies tools and every input file you plan to pass,
# and that they agree with each other (genome build, contig names, VEP cache version).
# Read-only except that it offers the exact command to fix each problem.
#
#   ./check_inputs.sh -i sample.vcf.gz -r GRCh38.fa --vep-cache ~/.vep [--clinvar clinvar.vcf.gz]
#                     [--exclude-bed x.bed] [--include-bed targets.bed] [--assembly GRCh38]
set -uo pipefail

VCF="" REF="" VEP_CACHE="" CLINVAR="" EXCLUDE_BED="" INCLUDE_BED="" ASSEMBLY="GRCh38"
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--vcf) VCF="$2"; shift ;;
        -r|--ref) REF="$2"; shift ;;
        --vep-cache) VEP_CACHE="$2"; shift ;;
        --clinvar) CLINVAR="$2"; shift ;;
        --exclude-bed) EXCLUDE_BED="$2"; shift ;;
        --include-bed) INCLUDE_BED="$2"; shift ;;
        --assembly) ASSEMBLY="$2"; shift ;;
        -h|--help) sed -n '2,8p' "$0"; exit 0 ;;
        *) shift; continue ;;     # ignore pipeline-only options so the same command line works
    esac
    shift
done

OK=0 BAD=0 WARN=0
ok()   { printf '  [ OK ] %s\n' "$*"; OK=$((OK+1)); }
bad()  { printf '  [FAIL] %s\n' "$*"; BAD=$((BAD+1)); }
wrn()  { printf '  [WARN] %s\n' "$*"; WARN=$((WARN+1)); }
fix()  { printf '         fix: %s\n' "$*"; }
hr()   { printf '\n== %s ==\n' "$1"; }
have() { command -v "$1" >/dev/null 2>&1; }
build_of() {   # chr1 length -> build name
    case "$1" in 248956422) echo GRCh38 ;; 249250621) echo GRCh37 ;; 248387328) echo T2T-CHM13 ;; "") echo unknown ;; *) echo "unknown($1)" ;; esac
}
style_of() { grep -qE '^chr' <<<"$1" && echo "chr-style (chr1)" || echo "number-style (1)"; }

# ------------------------------------------------------------------ tools
hr "Tools"
for t in bcftools bgzip tabix samtools vep; do
    if have "$t"; then ok "$t: $(command -v "$t")"; else
        [[ $t == samtools ]] && { wrn "samtools not found (only needed to index a FASTA without .fai)"; fix "conda install -c bioconda samtools"; continue; }
        [[ $t == vep ]] && { bad "vep not found"; fix "conda install -c bioconda ensembl-vep   (or run the pipeline with --skip-vep)"; continue; }
        bad "$t not found"; fix "conda install -c bioconda bcftools htslib"
    fi
done
if have bcftools; then
    BV=$(bcftools --version | head -1 | awk '{print $2}')
    awk -v v="$BV" 'BEGIN{split(v,a,"."); exit !(a[1]>1 || (a[1]==1 && a[2]>=13))}' && ok "bcftools version $BV (>= 1.13)" ||
        { bad "bcftools $BV is too old (need >= 1.13)"; fix "conda install -c bioconda 'bcftools>=1.17'"; }
else
    echo; echo "bcftools is required for the remaining checks."; exit 1
fi

# ------------------------------------------------------------------ VCF
hr "Input VCF (-i)"
VCF_LEN1="" VCF_STYLE=""
if [[ -z "$VCF" ]]; then bad "no -i given"
elif [[ ! -f "$VCF" ]]; then bad "not found: $VCF"
elif ! bcftools view -h "$VCF" >/dev/null 2>&1; then bad "cannot read as VCF/BCF: $VCF"
else
    ok "readable: $VCF ($(du -h "$VCF" | cut -f1))"
    if [[ -f "$VCF.tbi" || -f "$VCF.csi" ]]; then ok "index present"
    else wrn "no index (.tbi/.csi) – not required, but useful"; fix "bcftools index -t $VCF   (if it fails: bcftools sort $VCF -Oz -o sorted.vcf.gz && bcftools index -t sorted.vcf.gz)"; fi
    ns=$(bcftools query -l "$VCF" | grep -c .)
    case "$ns" in
        0) bad "no samples (sites-only VCF) – a genotyped VCF is needed" ;;
        1) ok "1 sample: $(bcftools query -l "$VCF")" ;;
        *) wrn "$ns samples – pass --sample NAME to the pipeline (first: $(bcftools query -l "$VCF" | head -3 | paste -sd, -))" ;;
    esac
    h=$(bcftools view -h "$VCF")
    grep -q '<NON_REF>' <<<"$h" && wrn "looks like a gVCF – genotype it first (gatk GenotypeGVCFs) for meaningful QUAL/filters"
    first=$(sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' <<<"$h" | head -1)
    [[ -z "$first" ]] && first=$(bcftools query -f '%CHROM\n' "$VCF" 2>/dev/null | head -1)
    VCF_STYLE=$(style_of "$first")
    VCF_LEN1=$(sed -n 's/^##contig=<ID=\(chr\)\{0,1\}1,.*length=\([0-9]*\).*/\2/p' <<<"$h" | head -1)
    ok "contig names: $VCF_STYLE; build from header contig lengths: $(build_of "$VCF_LEN1")"
    ref_line=$(grep -m1 '^##reference=' <<<"$h"); [[ -n "$ref_line" ]] && ok "header ${ref_line#\#\#} (the FASTA it was called on)"
    callers=$(grep -oE 'HaplotypeCaller|GenotypeGVCFs|DeepVariant|DRAGEN|Strelka|FreeBayes|bcftools_call|Clair3|Sentieon|Mutect2' <<<"$h" | sort -u | paste -sd, -)
    [[ -n "$callers" ]] && ok "caller(s) in header: $callers"
    grep -q 'Mutect2' <<<"$h" && wrn "Mutect2 = somatic VCF; vcf_hc_vep.sh is designed for germline calls"
    for fld in GT AD DP GQ; do
        grep -q "^##FORMAT=<ID=$fld," <<<"$h" && ok "FORMAT/$fld present" || wrn "FORMAT/$fld missing – the filters that use it are skipped"
    done
fi

# ------------------------------------------------------------------ FASTA
hr "Reference FASTA (-r)"
REF_LEN1="" REF_STYLE=""
if [[ -z "$REF" ]]; then bad "no -r given"
elif [[ ! -f "$REF" ]]; then bad "not found: $REF"
else
    ok "found: $REF ($(du -h "$REF" | cut -f1))"
    if [[ "$REF" == *.gz ]] && ! (have bgzip && bgzip -t "$REF" 2>/dev/null); then
        bad "gzip-compressed but not BGZF – cannot be indexed"; fix "gunzip -c $REF > ${REF%.gz} && samtools faidx ${REF%.gz}"
    fi
    if [[ -f "$REF.fai" ]]; then
        ok ".fai index present"
        REF_STYLE=$(style_of "$(head -1 "$REF.fai" | cut -f1)")
        REF_LEN1=$(awk '$1=="1"||$1=="chr1"{print $2; exit}' "$REF.fai")
        ok "contig names: $REF_STYLE; build from chr1 length: $(build_of "$REF_LEN1"); $(wc -l <"$REF.fai" | tr -d ' ') contigs"
    else
        bad "no $REF.fai"; fix "samtools faidx $REF"
    fi
fi

# ------------------------------------------------------------------ VCF vs FASTA
hr "VCF vs FASTA consistency"
if [[ -n "$VCF_LEN1" && -n "$REF_LEN1" ]]; then
    [[ "$VCF_LEN1" == "$REF_LEN1" ]] && ok "same build ($(build_of "$REF_LEN1"))" ||
        { bad "different builds: VCF $(build_of "$VCF_LEN1") vs FASTA $(build_of "$REF_LEN1")"; fix "use the FASTA the VCF was called on (see ##reference above), or liftover the VCF (e.g. bcftools +liftover / Picard LiftoverVcf)"; }
elif [[ -n "$REF_LEN1" ]]; then wrn "VCF header has no contig lengths – build compared by REF alleles only"; fi
[[ "$ASSEMBLY" == "$(build_of "$REF_LEN1")" || -z "$REF_LEN1" ]] || wrn "--assembly $ASSEMBLY but FASTA looks like $(build_of "$REF_LEN1")"
if [[ -n "$VCF_STYLE" && -n "$REF_STYLE" && "$VCF_STYLE" != "$REF_STYLE" ]]; then
    ok "contig naming differs (VCF $VCF_STYLE, FASTA $REF_STYLE) – the pipeline renames automatically"
fi
if [[ -f "$VCF" && -f "$REF.fai" ]]; then
    # REF allele spot-check on the first 20,000 records (contigs renamed to the FASTA's style)
    tmpmap=$(mktemp); awk '{c=$1; if (c ~ /^chr/) print substr(c,4)"\t"c; else print "chr"c"\t"c}' "$REF.fai" >"$tmpmap"
    res=$( { bcftools view -h "$VCF"; bcftools view -H "$VCF" 2>/dev/null | head -20000; } |
           bcftools annotate --rename-chrs "$tmpmap" 2>/dev/null |
           bcftools norm -f "$REF" -c w -Ou -o /dev/null 2>&1 >/dev/null | grep -c REF_MISMATCH)
    rm -f "$tmpmap"
    if [[ "${res:-0}" -eq 0 ]]; then ok "REF alleles match the FASTA (first 20,000 records)"
    elif [[ "$res" -lt 200 ]]; then wrn "$res REF mismatches in first 20,000 records (<1%; removed by the pipeline)"
    else bad "$res REF mismatches in first 20,000 records – wrong reference or build"; fix "use the FASTA the VCF was called on"; fi
fi

# ------------------------------------------------------------------ VEP cache
hr "VEP cache (--vep-cache)"
if [[ -z "$VEP_CACHE" ]]; then wrn "no --vep-cache given (pipeline needs --vep-cache, --vep-gtf or --skip-vep)"
elif [[ ! -d "$VEP_CACHE" ]]; then bad "directory not found: $VEP_CACHE"
else
    found=$(ls -d "$VEP_CACHE"/homo_sapiens*/*_"$ASSEMBLY" 2>/dev/null | xargs -n1 basename 2>/dev/null | paste -sd' ' -)
    [[ -n "$found" ]] && ok "cache(s) for $ASSEMBLY: $found" || bad "no homo_sapiens cache for $ASSEMBLY in $VEP_CACHE"
    if have vep; then
        vv=$(vep --help 2>/dev/null | sed -n 's/.*ensembl-vep *: *\([0-9]*\).*/\1/p' | head -1)
        if [[ -n "$vv" ]]; then
            if ls -d "$VEP_CACHE"/homo_sapiens*/"${vv}_${ASSEMBLY}" >/dev/null 2>&1; then ok "cache matches VEP version $vv"
            else
                bad "VEP is version $vv but no ${vv}_${ASSEMBLY} cache"
                fix "vep_install -a cf -s homo_sapiens -y $ASSEMBLY -c $VEP_CACHE"
                [[ "$ASSEMBLY" == GRCh37 ]] && base="https://ftp.ensembl.org/pub/grch37/release-${vv}" || base="https://ftp.ensembl.org/pub/release-${vv}"
                fix "or: curl -O ${base}/variation/indexed_vep_cache/homo_sapiens_vep_${vv}_${ASSEMBLY}.tar.gz && tar -xzf homo_sapiens_vep_${vv}_${ASSEMBLY}.tar.gz -C $VEP_CACHE"
            fi
        else wrn "could not read the VEP version (does 'vep --help' run?)"; fi
    fi
fi

# ------------------------------------------------------------------ ClinVar
hr "ClinVar (--clinvar, optional)"
if [[ -z "$CLINVAR" ]]; then wrn "not given (optional)"
    fix "curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${ASSEMBLY}/clinvar.vcf.gz && curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${ASSEMBLY}/clinvar.vcf.gz.tbi"
elif [[ ! -f "$CLINVAR" ]]; then bad "not found: $CLINVAR"
elif ! bcftools view -h "$CLINVAR" >/dev/null 2>&1; then bad "cannot read: $CLINVAR"
else
    ch=$(bcftools view -h "$CLINVAR")
    fd=$(grep -m1 '^##fileDate' <<<"$ch" | sed 's/^##fileDate=//')
    ok "readable${fd:+ (release date $fd)}"
    cref=$(grep -m1 '^##reference=' <<<"$ch" | sed 's/^##reference=//')
    [[ -n "$cref" ]] && { [[ "$cref" == *"${ASSEMBLY#GRCh}"* ]] && ok "ClinVar build: $cref" || bad "ClinVar build $cref does not match $ASSEMBLY"; }
    [[ -f "$CLINVAR.tbi" || -f "$CLINVAR.csi" ]] && ok "index present" || { bad "no index"; fix "tabix -p vcf $CLINVAR"; }
    grep -q '^##INFO=<ID=CLNSIG,' <<<"$ch" && ok "CLNSIG present" || bad "no INFO/CLNSIG – is this the ClinVar VCF?"
fi

# ------------------------------------------------------------------ BEDs
for kind in EXCLUDE_BED INCLUDE_BED; do
    f="${!kind}"
    hr "$( [[ $kind == EXCLUDE_BED ]] && echo 'Exclude BED (--exclude-bed, optional)' || echo 'Include BED (--include-bed, optional)')"
    if [[ -z "$f" ]]; then
        wrn "not given (optional)"
        [[ $kind == EXCLUDE_BED ]] && fix "ENCODE blacklist ($ASSEMBLY): curl -LO https://github.com/Boyle-Lab/Blacklist/raw/master/lists/$( [[ $ASSEMBLY == GRCh37 ]] && echo hg19 || echo hg38 )-blacklist.v2.bed.gz"
        [[ $kind == INCLUDE_BED ]] && fix "exome/panel: use the capture kit's target BED from the vendor (Twist / IDT / Agilent) for the same build"
        continue
    fi
    [[ -f "$f" ]] || { bad "not found: $f"; continue; }
    first=$( (gzip -cdf "$f" 2>/dev/null || cat "$f") | grep -vE '^(#|track|browser)' | head -1)
    nf=$(awk -F'\t' '{print NF}' <<<"$first")
    [[ "${nf:-0}" -ge 3 ]] && ok "BED with $nf columns; first contig $(cut -f1 <<<"$first") ($(style_of "$(cut -f1 <<<"$first")"); renamed automatically)" ||
        bad "first data line is not tab-separated BED (chrom start end): $first"
    n=$( (gzip -cdf "$f" 2>/dev/null || cat "$f") | grep -vcE '^(#|track|browser)')
    sz=$( (gzip -cdf "$f" 2>/dev/null || cat "$f") | grep -vE '^(#|track|browser)' | awk -F'\t' '{s+=$3-$2} END{printf "%.1f", s/1e6}')
    ok "$n intervals, $sz Mb"
done

hr "Summary"
printf '  %d OK, %d warnings, %d failures\n' "$OK" "$WARN" "$BAD"
[[ $BAD -eq 0 ]] && echo "  Ready to run vcf_hc_vep.sh with these inputs." || echo "  Fix the [FAIL] items first."
exit $(( BAD > 0 ))
