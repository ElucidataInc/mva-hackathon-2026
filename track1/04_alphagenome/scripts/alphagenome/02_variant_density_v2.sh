#!/usr/bin/env bash
# Patient variants in marker-gene regions -> per-gene counts and 5 kb window density.
set -euo pipefail
ROOT=/home/ubuntu/mva-hackathon-2026
STEP=$(cd "$(dirname "$0")/../.." && pwd)
export PATH=$ROOT/tools/pixi/.pixi/envs/qc/bin:$PATH
OUT=$STEP/results/alphagenome/v2_team11
VCF=$ROOT/data/variants/WGS_EX2312012_HGWCNDSX7.vcf.gz
BED=$OUT/marker_genes_flank50kb.bed
WIN=${WIN:-5000}
VAR=$OUT/marker_region_variants.tsv          # variant-level: stays on server
DEN=$OUT/variant_density_${WIN}bp.tsv
GEN=$OUT/variant_counts_per_gene.tsv

printf "gene\tchrom\tpos\tref\talt\ttype\tfilter\tgt\tregion_part\n" > "$VAR"
while IFS=$'\t' read -r chr s e gene gid strand gs ge; do
  bcftools view -r "$chr:$((s+1))-$e" -i 'GT="alt"' "$VCF" -Ou \
  | bcftools norm -m -any -Ou \
  | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%FILTER\t[%GT]\n" \
  | awk -F'\t' -v OFS='\t' -v g="$gene" -v gs="$gs" -v ge="$ge" -v st="$strand" '{
      if($2<=gs) part=(st=="+")?"upstream":"downstream";
      else if($2>ge) part=(st=="+")?"downstream":"upstream";
      else part="gene_body"; print g,$0,part}' >> "$VAR"
done < "$BED"

# 5 kb windows
printf "gene\tchrom\twin_start\twin_end\tn_total\tn_snv\tn_indel_other\tn_pass\n" > "$DEN"
awk -F'\t' -v OFS='\t' -v W="$WIN" 'NR>1{w=int(($3-1)/W)*W; k=$1 OFS $2 OFS w OFS w+W
  n[k]++; if($6=="SNP") s[k]++; else o[k]++; if($7=="PASS") p[k]++}
  END{for(k in n) print k,n[k],s[k]+0,o[k]+0,p[k]+0}' "$VAR" | sort -t$'\t' -k1,1 -k5,5nr >> "$DEN"

# per-gene summary (normalised per 10 kb of region)
printf "gene\tregion_kb\tn_total\tn_snv\tn_indel_other\tn_pass\tupstream\tgene_body\tdownstream\tper_10kb\n" > "$GEN"
awk -F'\t' -v OFS='\t' 'FNR==NR{len[$4]=($3-$2)/1000; next}
  FNR>1{n[$1]++; if($6=="SNP")s[$1]++; else o[$1]++; if($7=="PASS")p[$1]++; part[$1,$9]++}
  END{for(g in len) printf "%s\t%.0f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", g,len[g],n[g],s[g],o[g],p[g],
      part[g,"upstream"],part[g,"gene_body"],part[g,"downstream"],(n[g]/len[g])*10}' "$BED" "$VAR" \
  | sort -t$'\t' -k10,10nr >> "$GEN"

echo "== per-gene counts =="; column -t -s$'\t' "$GEN"
echo; echo "== top 3 windows per gene (counts only) =="
awk -F'\t' 'NR>1 && c[$1]++<3' "$DEN" | cut -f1,3-8 | column -t
echo; echo "variant rows written: $(($(wc -l < "$VAR")-1))  (file kept on server, not for sharing)"
