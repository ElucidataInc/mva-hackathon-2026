#!/usr/bin/env bash
# Module 01 – inspect: size, checksum, header, samples, contig record counts.

mod_inspect() {
    local d="${OUTDIR}/01_inspect"; mkdir -p "$d"
    {
        echo "file:        $VCF_ORIG"
        echo "size_bytes:  $(wc -c <"$VCF_ORIG" | tr -d ' ')"
        echo "size_human:  $(du -h "$VCF_ORIG" | cut -f1)"
        echo "modified:    $(date -r "$VCF_ORIG" 2>/dev/null)"
        if have md5sum; then echo "md5:         $(md5sum "$VCF_ORIG" | cut -d' ' -f1)"
        elif have md5; then echo "md5:         $(md5 -q "$VCF_ORIG")"; fi
        echo "samples:     $(n_samples)"
    } | tee "$d/file_info.txt"
    bcf_header >"$d/header.txt"
    grep -E '^##(source|reference|.*Command|.*[Vv]ersion|fileDate)' "$d/header.txt" >"$d/provenance_lines.txt" || true
    if have bcftools; then
        bcftools query -l "$VCF_IN" >"$d/samples.txt" 2>/dev/null || true
        bcftools index -s "$VCF_IN" >"$d/records_per_contig.tsv" 2>/dev/null || true
    else
        awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++)print $i; exit}' "$d/header.txt" >"$d/samples.txt"
    fi
    return 0
}
