#!/usr/bin/env bash
# Module 06 – annotation: dbSNP IDs + gnomAD AF (bcftools annotate) for known-vs-novel rates,
# functional consequences with VEP (preferred) or SnpEff, then annotation_summary.py.

mod_annotation() {
    local d="${OUTDIR}/06_annotation"; mkdir -p "$d"
    if [[ -z "$DBSNP$GNOMAD$VEP_CACHE$SNPEFF_DB" ]]; then
        mnote "no --dbsnp/--gnomad/--vep-cache/--snpeff-db given"; return 3
    fi
    have bcftools || { mnote "bcftools not found"; return 3; }
    local cur="$VCF_IN"
    if [[ -n "$DBSNP" ]]; then
        bcftools annotate --threads "$THREADS" -a "$DBSNP" -c ID -Oz -o "$d/step_dbsnp.vcf.gz" "$cur" && \
            bcftools index -t "$d/step_dbsnp.vcf.gz" && cur="$d/step_dbsnp.vcf.gz" || mnote "dbSNP annotation failed"
    fi
    if [[ -n "$GNOMAD" ]]; then
        bcftools annotate --threads "$THREADS" -a "$GNOMAD" -c "INFO/gnomAD_AF:=INFO/${GNOMAD_AF_TAG}" \
            -Oz -o "$d/step_gnomad.vcf.gz" "$cur" && bcftools index -t "$d/step_gnomad.vcf.gz" && \
            cur="$d/step_gnomad.vcf.gz" || mnote "gnomAD annotation failed"
    fi
    if [[ -n "$VEP_CACHE" ]] && have vep; then
        vep --input_file "$cur" --output_file "$d/annotated.vep.vcf.gz" --vcf --compress_output bgzip \
            --offline --cache --dir_cache "$VEP_CACHE" --assembly "$VEP_ASSEMBLY" ${REF:+--fasta "$REF"} \
            --fork "$THREADS" --symbol --biotype --canonical --check_existing --af_gnomadg --pick_allele \
            --stats_file "$d/vep_summary.html" --force_overwrite >"$d/vep.log" 2>&1 && \
            cur="$d/annotated.vep.vcf.gz" || mnote "VEP failed (see vep.log)"
    elif [[ -n "$SNPEFF_DB" ]] && have snpEff; then
        snpEff -Xmx8g -stats "$d/snpeff_summary.html" "$SNPEFF_DB" "$cur" 2>"$d/snpeff.log" | \
            bgzip -c >"$d/annotated.snpeff.vcf.gz" && cur="$d/annotated.snpeff.vcf.gz" || mnote "SnpEff failed"
    elif [[ -n "$VEP_CACHE$SNPEFF_DB" ]]; then
        mnote "vep/snpEff executable not found"
    fi
    ln -sf "$(basename "$cur")" "$d/annotated.vcf.gz" 2>/dev/null || cp "$cur" "$d/annotated.vcf.gz"
    "$PYTHON" "${SCRIPTS}/annotation_summary.py" --vcf "$cur" --out "$d/annotation_summary.json" --gnomad-tag gnomAD_AF
}
