#!/usr/bin/env bash
# Build BED of literature marker genes (gene body +/- FLANK) from Ensembl 116 GTF.
set -euo pipefail
ROOT=/home/ubuntu/mva-hackathon-2026
STEP=$(cd "$(dirname "$0")/../.." && pwd)
OUT=$STEP/results/alphagenome
GTF=$ROOT/references/annotation/ensembl_116/Homo_sapiens.GRCh38.116.gtf.gz
FLANK=${FLANK:-50000}
GENES="BUB1B CEP57 TRIP13 CENATAC CCDC84 SLF2 SMC5 MAD1L1 BUB1 MAD2L1BP CEP192"
mkdir -p "$OUT"
echo "GTF=$GTF  FLANK=$FLANK"
zcat "$GTF" | awk -F'\t' -v genes="$GENES" -v f="$FLANK" '
  BEGIN{n=split(genes,a," "); for(i=1;i<=n;i++) g[a[i]]=1}
  $3=="gene" && match($9,/gene_name "[^"]+"/){
    nm=substr($9,RSTART+11,RLENGTH-12)
    if(nm in g){ match($9,/gene_id "[^"]+"/); id=substr($9,RSTART+9,RLENGTH-10)
      s=$4-1-f; if(s<0) s=0
      print $1"\t"s"\t"$5+f"\t"nm"\t"id"\t"$7"\t"$4-1"\t"$5 } }' \
  | sort -k1,1V -k2,2n > "$OUT/marker_genes_flank50kb.bed"
echo "# BED cols: chrom start end gene gene_id strand gene_start0 gene_end"
cat "$OUT/marker_genes_flank50kb.bed"
echo "# genes found: $(cut -f4 "$OUT/marker_genes_flank50kb.bed" | sort -u | wc -l) / 10 expected"
for g in $GENES; do grep -qw "$g" "$OUT/marker_genes_flank50kb.bed" || echo "MISSING: $g"; done
