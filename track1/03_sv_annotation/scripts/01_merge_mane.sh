#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
Q='%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\n'
SIZE="INFO/SVLEN>=$SVMIN || INFO/SVLEN<=-$SVMIN"
bcftools query -i "$SIZE" -f "$Q" "$IN/manta/manta.pass.vcf.gz" | awk -v OFS='\t' '{if ($5 < 0) $5 = -$5} 1' > "$OUT/manta_ge4kb.tsv"
bcftools query -i "$SIZE" -f "$Q" "$IN/cue/reports/svs.vcf" > "$OUT/cue_ge4kb.tsv"
{ awk -v OFS='\t' '{print $0, "manta"}' "$OUT/manta_ge4kb.tsv"; awk -v OFS='\t' '{print $0, "cue"}' "$OUT/cue_ge4kb.tsv"; } \
    | sort -k1,1V -k2,2n -k3,3n > "$OUT/combined_ge4kb.tsv"

{
    printf '%s\n' '##fileformat=VCFv4.2' '##FILTER=<ID=PASS,Description="All filters passed">'
    bcftools view -h "$IN/manta/manta.pass.vcf.gz" | grep '^##contig'
    printf '%s\n' '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">' \
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">' \
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">' \
        '##INFO=<ID=SRC,Number=1,Type=String,Description="Caller">' \
        '##INFO=<ID=CALLER_TYPE,Number=1,Type=String,Description="Original type label">'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    awk -v OFS='\t' '{t = $4 == "IDUP" ? "DUP" : $4; print $1, $2, $6 "_" NR, "N", "<" t ">", ".", "PASS", "SVTYPE=" t ";END=" $3 ";SVLEN=" $5 ";SRC=" $6 ";CALLER_TYPE=" $4}' "$OUT/combined_ge4kb.tsv"
} | bcftools sort -Oz -o "$OUT/combined.vcf.gz" -
bcftools index -ft "$OUT/combined.vcf.gz"

zcat "$SVREF/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz" \
    | awk -F'\t' -v OFS='\t' '$3 == "transcript" && $9 ~ /tag "MANE_Select"/ {match($9, /gene_name "[^"]+"/); print $1, $4 - 1, $5, substr($9, RSTART + 11, RLENGTH - 12)}' \
    | sort -k1,1 -k2,2n > "$OUT/mane_genes.bed"
bcftools query -f '%CHROM\t%POS0\t%INFO/END\t%ID\n' "$OUT/combined.vcf.gz" | sort -k1,1 -k2,2n > "$OUT/sv.bed"
bedtools map -a "$OUT/sv.bed" -b "$OUT/mane_genes.bed" -c 4 -o distinct \
    | awk -v OFS='\t' '$5 != "." {print $1, $2 + 1, $5}' | bgzip > "$OUT/mane.tab.gz"
tabix -f -s1 -b2 -e2 "$OUT/mane.tab.gz"
echo '##INFO=<ID=MANE_GENE,Number=.,Type=String,Description="MANE Select genes overlapped">' > "$OUT/mane.hdr"
bcftools annotate -a "$OUT/mane.tab.gz" -h "$OUT/mane.hdr" -c CHROM,POS,INFO/MANE_GENE -Oz -o "$OUT/combined.mane.vcf.gz" "$OUT/combined.vcf.gz"
