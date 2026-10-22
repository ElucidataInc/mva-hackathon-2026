#!/usr/bin/env bash
# Module 03 – bcftools stats -s - (all samples) [+ -F ref for indel context], PASS-only stats,
# and plot-vcfstats.

mod_stats() {
    local d="${OUTDIR}/03_bcftools_stats"; mkdir -p "$d"
    have bcftools || { mnote "bcftools not found"; return 3; }
    local extra=()
    [[ -n "$REF" ]] && extra+=(-F "$REF")
    [[ -n "$TARGETS_BED" ]] && extra+=(-R "$TARGETS_BED")
    # -s - : per-sample stats for every sample (PSC/PSI sections)
    local src="$VCF_IN" attempt err
    for attempt in 1 2 3; do
        if bcftools stats --threads "$THREADS" -s - -d 0,500,1 ${extra[@]+"${extra[@]}"} "$src" \
               >"$d/bcftools.stats.txt" 2>"$d/stats.err"; then
            break
        fi
        err=$(grep -m1 -E 'differs|Incorrect AN/AC|E::' "$d/stats.err")
        if [[ -n "$REF" && " ${extra[*]-} " == *" -F "* ]] && \
           ! grep -q 'Incorrect AN/AC' "$d/stats.err"; then
            # -F aborts on the first REF mismatch, on contig names absent from the FASTA (chr1 vs 1)
            # or on a missing/unreadable .fai; REF problems are reported by module 02 / profiler
            err=$(grep -m1 -v '^\s*$' "$d/stats.err")
            mnote "bcftools stats -F failed (${err}); rerun without the reference"
            extra=(); [[ -n "$TARGETS_BED" ]] && extra+=(-R "$TARGETS_BED")
        elif grep -q 'Incorrect AN/AC' "$d/stats.err" && [[ "$src" == "$VCF_IN" ]]; then
            # INFO/AC,AN inconsistent with genotypes (e.g. after subsetting samples)
            mnote "bcftools stats failed (${err}); AC/AN recomputed with +fill-tags"
            bcftools +fill-tags "$VCF_IN" -Oz -o "$d/input.filltags.vcf.gz" -- -t AC,AN && \
                bcftools index -t -f "$d/input.filltags.vcf.gz" || return 1
            src="$d/input.filltags.vcf.gz"
        else
            cat "$d/stats.err"; return 1
        fi
        [[ $attempt -eq 3 ]] && { cat "$d/stats.err"; return 1; }
    done
    bcftools stats --threads "$THREADS" -s - -f PASS,. ${extra[@]+"${extra[@]}"} "$src" >"$d/bcftools.stats.PASS.txt" || true
    if have plot-vcfstats; then
        # -P: skip the PDF (no LaTeX needed); PNGs + python scripts are kept
        plot-vcfstats -P -p "$d/plots/" "$d/bcftools.stats.txt" >"$d/plot-vcfstats.log" 2>&1 || \
            mnote "plot-vcfstats failed (needs python matplotlib); figures re-drawn in report"
    else
        mnote "plot-vcfstats not found"
    fi
    return 0
}
