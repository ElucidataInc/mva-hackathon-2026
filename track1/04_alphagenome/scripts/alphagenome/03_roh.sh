#!/usr/bin/env bash
# Genome-wide runs of homozygosity (bcftools roh); report overlap with marker genes.
set -euo pipefail
ROOT=/home/ubuntu/mva-hackathon-2026
STEP=$(cd "$(dirname "$0")/../.." && pwd)
export PATH=$ROOT/tools/pixi/.pixi/envs/qc/bin:$PATH
OUT=$STEP/results/alphagenome
VCF=$ROOT/data/variants/WGS_EX2312012_HGWCNDSX7.vcf.gz
BED=$OUT/marker_genes_flank50kb.bed
ROH=$OUT/roh_regions.tsv

# autosomes only, PASS, biallelic SNVs; default alt-AF 0.4 (no population AF tag in VCF)
bcftools view -f PASS -v snps -m2 -M2 -t ^X,Y,MT "$VCF" -Ou \
 | bcftools roh -G30 --AF-dflt 0.4 -O r -o "$OUT/roh_raw.txt"
printf "chrom\tstart\tend\tlength_bp\tn_markers\tquality\n" > "$ROH"
awk -F'\t' -v OFS='\t' '$1=="RG"{print $3,$4,$5,$6,$7,$8}' "$OUT/roh_raw.txt" >> "$ROH"

echo "== genome-wide ROH summary =="
awk -F'\t' 'NR>1{n++; t+=$4; if($4>=1e6){n1++; t1+=$4}; if($4>=5e6) n5++}
  END{printf "segments: %d | total: %.1f Mb | >=1 Mb: %d segments, %.1f Mb | >=5 Mb: %d\n", n,t/1e6,n1,t1/1e6,n5}' "$ROH"
echo "== ROH >=1 Mb per chromosome (Mb) =="
awk -F'\t' 'NR>1 && $4>=1e6{s[$1]+=$4} END{for(c in s) printf "%s\t%.1f\n",c,s[c]/1e6}' "$ROH" | sort -k1,1V | column -t
echo "== marker genes inside ROH (segment length) =="
awk -F'\t' 'FNR==NR{if(FNR>1){c[++n]=$1; s[n]=$2; e[n]=$3; l[n]=$4}; next}
  {hit="none"; for(i=1;i<=n;i++) if(c[i]==$1 && s[i]<=$8 && e[i]>=$7+1){hit=sprintf("%.2f Mb", l[i]/1e6)}
   print $4"\t"$1"\t"hit}' "$ROH" "$BED" | column -t
