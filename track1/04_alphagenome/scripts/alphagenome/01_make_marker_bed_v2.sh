#!/usr/bin/env bash
# v2: marker genes from team list (01_alignment_qc_calling/scripts/marker_genes.tsv), gene body +/- FLANK, Ensembl 116.
set -euo pipefail
ROOT=/home/ubuntu/mva-hackathon-2026
STEP=$(cd "$(dirname "$0")/../.." && pwd)
OUT=$STEP/results/alphagenome/v2_team11
GTF=$ROOT/references/annotation/ensembl_116/Homo_sapiens.GRCh38.116.gtf.gz
GENEFILE=$STEP/../01_alignment_qc_calling/scripts/marker_genes.tsv
FLANK=${FLANK:-50000}
mkdir -p "$OUT/logs"
GENES=$(awk -F'\t' 'NR>1 && $1!=""{print $1}' "$GENEFILE" | tr '\n' ' ')
N=$(echo $GENES | wc -w)
echo "GENEFILE=$GENEFILE  ($N genes): $GENES"
zcat "$GTF" | awk -F'\t' -v genes="$GENES" -v f="$FLANK" '
  BEGIN{n=split(genes,a," "); for(i=1;i<=n;i++) g[a[i]]=1}
  $3=="gene" && match($9,/gene_name "[^"]+"/){
    nm=substr($9,RSTART+11,RLENGTH-12)
    if(nm in g){ match($9,/gene_id "[^"]+"/); id=substr($9,RSTART+9,RLENGTH-10)
      s=$4-1-f; if(s<0) s=0
      print $1"\t"s"\t"$5+f"\t"nm"\t"id"\t"$7"\t"$4-1"\t"$5 } }' \
  | sort -k1,1V -k2,2n > "$OUT/marker_genes_flank50kb.bed"
awk -v OFS='\t' '{$1="chr"$1; print}' "$OUT/marker_genes_flank50kb.bed" > "$OUT/marker_genes_flank50kb.chr.bed"
cat "$OUT/marker_genes_flank50kb.bed"
echo "# genes found: $(cut -f4 "$OUT/marker_genes_flank50kb.bed" | sort -u | wc -l) / $N expected"
for g in $GENES; do grep -qw "$g" "$OUT/marker_genes_flank50kb.bed" || echo "MISSING: $g"; done
