#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
bcftools sort -Oz -o "$CR/svs.sorted.vcf.gz" "$CR/svs.vcf"
bcftools index -f "$CR/svs.sorted.vcf.gz"

READ_LEN=$(samtools view -F 0x904 "$FINAL_BAM" chr1:10000000-10100000 \
    | awk '{n[length($10)]++} END {for (l in n) if (n[l] > m) {m = n[l]; L = l}; print L}')
samtools idxstats "$FINAL_BAM" \
    | awk -v L="$READ_LEN" '$1 ~ /^chr([0-9]+|X|Y)$/ {printf "%s\t%.2f\n", $1, $3 * L / $2}' > "$CR/chrdepth.txt"

F='%CHROM\t%POS0\t%INFO/END\t%INFO/SVTYPE\t%QUAL\t[%GT]\n'
bcftools query -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP"' -f "$F" "$CR/svs.sorted.vcf.gz" > "$CR/dd.bed"
bcftools query -i 'INFO/SVTYPE="IDUP"' -f "$F" "$CR/svs.sorted.vcf.gz" > "$CR/idup.bed"
for b in dd idup; do
    samtools bedcov -Q 20 "$CR/$b.bed" "$FINAL_BAM" \
        | awk 'NR == FNR {d[$1] = $2; next} {printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n", $1, $2, $3, $4, $5, $6, $7 / ($3 - $2) / d[$1]}' "$CR/chrdepth.txt" - \
        > "$CR/$b.ratio.tsv"
done

awk -v OFS='\t' '($4 == "DEL" && $7 <= 0.7) || ($4 == "DUP" && $7 >= 1.3) {print $1, $2 + 1}' "$CR/dd.ratio.tsv" > "$CR/pass.pos"
bcftools view -T "$CR/pass.pos" -Oz -o "$CR/svs.depthvalidated.vcf.gz" "$CR/svs.sorted.vcf.gz"
bcftools view -i 'INFO/SVTYPE="INV"' -Oz -o "$C/cue.inv.vcf.gz" "$CR/svs.sorted.vcf.gz"
for v in "$CR/svs.depthvalidated.vcf.gz" "$C/cue.inv.vcf.gz"; do bcftools index -ft "$v"; done
