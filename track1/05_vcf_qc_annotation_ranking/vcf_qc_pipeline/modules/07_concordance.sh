#!/usr/bin/env bash
# Module 07 – truth-set concordance with hap.py and/or RTG vcfeval (e.g. GIAB HG002).

mod_concordance() {
    local d="${OUTDIR}/07_concordance"; mkdir -p "$d"
    [[ -n "$TRUTH" ]] || { mnote "no --truth VCF"; return 3; }
    [[ -n "$REF" ]] || { mnote "--ref required"; return 3; }
    local q="$VCF_IN"
    if [[ -n "$QUERY_SAMPLE" ]] && have bcftools; then
        bcftools view -s "$QUERY_SAMPLE" -c1 -Oz -o "$d/query.vcf.gz" "$VCF_IN" && bcftools index -t "$d/query.vcf.gz" && q="$d/query.vcf.gz"
    fi
    local ran=0
    if have hap.py; then
        mkdir -p "$d/happy"
        hap.py "$TRUTH" "$q" -r "$REF" ${TRUTH_BED:+-f "$TRUTH_BED"} ${TARGETS_BED:+-T "$TARGETS_BED"} \
            --threads "$THREADS" -o "$d/happy/happy" >"$d/happy/happy.log" 2>&1 && ran=1 || mnote "hap.py failed"
    fi
    if have rtg && [[ -n "$RTG_SDF" ]]; then
        rm -rf "$d/vcfeval"
        rtg vcfeval -b "$TRUTH" -c "$q" -t "$RTG_SDF" ${TRUTH_BED:+-e "$TRUTH_BED"} --threads "$THREADS" \
            -o "$d/vcfeval" >"$d/vcfeval.log" 2>&1 && ran=1 || mnote "rtg vcfeval failed"
    fi
    [[ $ran -eq 1 ]] || { mnote "neither hap.py nor rtg (+ --rtg-sdf) available"; return 3; }
    return 0
}
