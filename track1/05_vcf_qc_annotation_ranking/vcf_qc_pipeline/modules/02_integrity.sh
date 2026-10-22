#!/usr/bin/env bash
# Module 02 – integrity: index/contig check, REF concordance (bcftools norm -c w),
# normalisation summary and GATK ValidateVariants.
# (Sort order and exact duplicates are measured by the profiler in module 04.)

mod_integrity() {
    local d="${OUTDIR}/02_integrity"; mkdir -p "$d"
    have bcftools || { mnote "bcftools not found"; return 3; }
    bcftools index -s "$VCF_IN" >"$d/index_check.txt" 2>&1 || true
    if [[ -n "$REF" ]]; then
        # -c w: warn on REF mismatch; stderr also gets the Lines total/split/realigned/skipped summary
        # contig-name compatibility between VCF and FASTA (chr1 vs 1 is the usual problem)
        local fai="${REF}.fai"
        [[ -f "$fai" ]] || { have samtools && samtools faidx "$REF"; }
        if [[ -f "$fai" ]]; then
            bcftools index -s "$VCF_IN" 2>/dev/null | cut -f1 | sort >"$d/vcf_contigs.txt"
            cut -f1 "$fai" | sort >"$d/fasta_contigs.txt"
            local nv nm
            nv=$(wc -l <"$d/vcf_contigs.txt" | tr -d ' ')
            nm=$(comm -23 "$d/vcf_contigs.txt" "$d/fasta_contigs.txt" | tee "$d/contigs_missing_in_fasta.txt" | wc -l | tr -d ' ')
            [[ "$nm" -gt 0 ]] && mnote "$nm of $nv VCF contigs absent from FASTA (e.g. $(head -2 "$d/contigs_missing_in_fasta.txt" | paste -sd, -) vs FASTA $(head -2 "$d/fasta_contigs.txt" | paste -sd, -))"
        else
            mnote "no ${REF}.fai and samtools not found: run 'samtools faidx $REF'"
        fi
        bcftools norm --threads "$THREADS" -f "$REF" -c w -Ou -o /dev/null "$VCF_IN" 2>"$d/norm_checkref.log" || \
            mnote "bcftools norm failed: $(grep -m1 -iE 'error|fail|not' "$d/norm_checkref.log" | cut -c1-160)"
        # how many records would change on left-alignment / splitting
        bcftools norm --threads "$THREADS" -f "$REF" -c x -m -any -Ou -o /dev/null "$VCF_IN" 2>"$d/norm_split_leftalign.log" || true
    else
        mnote "no --ref: REF concordance and left-alignment checks skipped"
    fi
    if have gatk && [[ -n "$REF" ]]; then
        local dict="${REF%.*}.dict"; [[ "$REF" == *.gz ]] && dict="${REF%.*.*}.dict"
        if [[ -f "$dict" ]]; then
            gatk ValidateVariants -R "$REF" -V "$VCF_IN" --warn-on-errors \
                ${DBSNP:+--dbsnp "$DBSNP"} >"$d/gatk_ValidateVariants.log" 2>&1 || \
                mnote "ValidateVariants reported errors"
        else
            mnote "ValidateVariants skipped: sequence dictionary $dict missing (gatk CreateSequenceDictionary)"
        fi
    else
        mnote "GATK ValidateVariants skipped (needs gatk + --ref)"
    fi
    return 0
}
