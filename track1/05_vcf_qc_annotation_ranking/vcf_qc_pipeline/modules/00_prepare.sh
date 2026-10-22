#!/usr/bin/env bash
# Module 00 – prepare: tool inventory, bgzip + index the input, sorted copy if needed.
# Exports VCF_IN (bgzipped, indexed when bcftools is available) and PROFILE_VCF.

mod_prepare() {
    local d="${OUTDIR}/00_prepare"; mkdir -p "$d"
    {
        echo "# tool versions ($(date))"
        for t in bcftools tabix bgzip vcftools plot-vcfstats gatk somalier peddy vep snpEff hap.py rtg python3; do
            if have "$t"; then
                case $t in
                    bcftools|tabix|bgzip) v=$($t --version 2>&1 | head -1) ;;
                    vcftools) v=$(vcftools --version 2>&1 | head -1) ;;
                    gatk) v=$(gatk --version 2>/dev/null | grep -m1 -i 'gatk\|version') ;;
                    vep) v=$(vep --help 2>/dev/null | grep -m1 -i 'ensembl-vep' | xargs) ;;
                    snpEff) v=$(snpEff -version 2>&1 | head -1) ;;
                    somalier) v=$(somalier 2>&1 | grep -m1 -i version) ;;
                    peddy) v=$(peddy --version 2>&1 | head -1) ;;
                    hap.py) v="hap.py ($(command -v hap.py))" ;;
                    rtg) v=$(rtg version 2>/dev/null | head -1) ;;
                    python3) v=$(python3 --version 2>&1) ;;
                    *) v="available" ;;
                esac
                printf '%-14s %s\n' "$t" "${v:-available}"
            else
                printf '%-14s %s\n' "$t" "NOT FOUND"
            fi
        done
    } >"${LOGDIR}/tool_versions.txt"
    cat "${LOGDIR}/tool_versions.txt"

    if ! have bcftools; then
        mnote "bcftools not found: bcftools-based modules will be skipped; Python profiler used on raw input"
        return 0
    fi
    local src="$VCF_ORIG" out="$d/input.vcf.gz" abs
    abs="$(cd "$(dirname "$src")" && pwd)/$(basename "$src")"
    rm -f "$d"/input.*
    if [[ "$src" == *.vcf.gz && ( -f "$src.tbi" || -f "$src.csi" ) ]]; then
        ln -sf "$abs" "$out"
        for ix in tbi csi; do [[ -f "$src.$ix" ]] && ln -sf "$abs.$ix" "$out.$ix"; done
    else
        if [[ "$src" == *.vcf.gz ]]; then ln -sf "$abs" "$out"; fi
        if [[ ! -e "$out" ]] || ! bcftools index -t "$out" 2>"$d/index.err"; then
            rm -f "$out"
            log "re-compressing input with BGZF -> $out"
            bcftools view --threads "$THREADS" -Oz -o "$out" "$src" || return 1
            if ! bcftools index -t "$out" 2>"$d/index.err"; then
                mnote "input not indexable (unsorted?) - tool modules use a sorted copy"
                bcftools sort -T "$d/sort_tmp" -Oz -o "$d/input.sorted.vcf.gz" "$out" || return 1
                bcftools index -t "$d/input.sorted.vcf.gz" || return 1
                ln -sf input.sorted.vcf.gz "$d/input.use.vcf.gz"
                ln -sf input.sorted.vcf.gz.tbi "$d/input.use.vcf.gz.tbi"
            fi
        fi
    fi
    if [[ ! -e "$d/input.use.vcf.gz" ]]; then
        ln -sf input.vcf.gz "$d/input.use.vcf.gz"
        for ix in tbi csi; do [[ -e "$out.$ix" ]] && ln -sf "input.vcf.gz.$ix" "$d/input.use.vcf.gz.$ix"; done
    fi
    return 0
}
