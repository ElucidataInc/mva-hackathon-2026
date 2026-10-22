#!/usr/bin/env bash
# Marker-region variants from the step 05 VEP output -> one row per variant (marker-gene CSQ),
# rarity flag (MAX_AF), high-confidence membership; prints counts only.
# NOTE: raw VCF uses "15"; the annotated VCFs use "chr15" -> chr-prefixed BED for those.
set -euo pipefail
ROOT=/home/ubuntu/mva-hackathon-2026
STEP=$(cd "$(dirname "$0")/../.." && pwd)
export PATH=$ROOT/tools/pixi/.pixi/envs/qc/bin:$PATH
OUT=$STEP/results/alphagenome
VEPDIR=$STEP/../05_vcf_qc_annotation_ranking/results/hc_out/04_vep
RAW=$ROOT/data/variants/WGS_EX2312012_HGWCNDSX7.vcf.gz
ALL=$VEPDIR/annotated.vcf.gz
HC=$VEPDIR/annotated.high_confidence.vcf.gz
BED=$OUT/marker_genes_flank50kb.bed
BEDC=$OUT/marker_genes_flank50kb.chr.bed
TSV=$OUT/marker_region_vep.tsv          # variant-level: stays on server
MAXAF=${MAXAF:-0.01}

awk -v OFS='\t' '{$1="chr"$1; print}' "$BED" > "$BEDC"

echo "== records with non-ref GT in marker regions =="
printf "%-40s %s\n" "$(basename "$RAW")" "$(bcftools view -H -R "$BED"  -i 'GT="alt"' "$RAW" | wc -l)"
printf "%-40s %s\n" "$(basename "$ALL")" "$(bcftools view -H -R "$BEDC" -i 'GT="alt"' "$ALL" | wc -l)"
printf "%-40s %s\n" "$(basename "$HC")"  "$(bcftools view -H -R "$BEDC" -i 'GT="alt"' "$HC"  | wc -l)"

bcftools query -R "$BEDC" -f '%CHROM:%POS:%REF:%ALT\n' "$HC" > "$OUT/tmp_hc_keys.txt"

FMT='%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT]\t%Allele\t%SYMBOL\t%Gene\t%Consequence\t%IMPACT\t%MANE_SELECT\t%PICK\t%HGVSc\t%HGVSp\t%gnomADg_AF\t%gnomADe_AF\t%MAX_AF\t%CLIN_SIG\t%am_class\t%am_pathogenicity\t%EXON\t%INTRON\t%DISTANCE\n'
bcftools view -R "$BEDC" -i 'GT="alt"' "$ALL" -Ou \
 | bcftools +split-vep -d -f "$FMT" > "$OUT/tmp_csq.tsv"

printf "marker_gene\tregion_part\tchrom\tpos\tref\talt\tfilter\tgt\tzygosity\tsymbol\tconsequence\timpact\tmane\thgvsc\thgvsp\tgnomADg_AF\tgnomADe_AF\tMAX_AF\trare\tclin_sig\tam_class\tam_path\texon\tintron\tdistance\tin_high_conf\n" > "$TSV"
awk -F'\t' -v OFS='\t' -v M="$MAXAF" '
  FILENAME==ARGV[1]{c[++n]=$1; s[n]=$2; e[n]=$3; g[n]=$4; st[n]=$6; gs[n]=$7; ge[n]=$8; next}
  FILENAME==ARGV[2]{hc[$1]=1; next}
  { mg=""; for(i=1;i<=n;i++) if(c[i]==$1 && $2>s[i] && $2<=e[i]){mg=g[i]; k=i}
    if(mg=="") next
    key=$1":"$2":"$3":"$4
    sc=($8==mg && $12!=".")?3:($8==mg)?2:($13=="1")?1:0
    if(!(key in best) || sc>best[key]){ best[key]=sc
      part=($2<=gs[k])?(st[k]=="+"?"upstream":"downstream"):($2>ge[k])?(st[k]=="+"?"downstream":"upstream"):"gene_body"
      gt=$6; gsub(/\|/,"/",gt); zyg=(gt ~ /^([1-9])\/\1$/)?"hom":"het"
      rare=($18=="." || $18+0<M)?"yes":"no"
      row[key]=mg OFS part OFS $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS zyg OFS $8 OFS $10 OFS $11 OFS $12 OFS $14 OFS $15 OFS $16 OFS $17 OFS $18 OFS rare OFS $19 OFS $20 OFS $21 OFS $22 OFS $23 OFS $24 OFS ((key in hc)?"yes":"no") } }
  END{for(k in row) print row[k]}' "$BEDC" "$OUT/tmp_hc_keys.txt" "$OUT/tmp_csq.tsv" \
 | sort -t$'\t' -k3,3V -k4,4n >> "$TSV"
rm -f "$OUT/tmp_csq.tsv" "$OUT/tmp_hc_keys.txt"

echo; echo "== per gene (rare = MAX_AF < $MAXAF or absent) =="
awk -F'\t' 'NR>1{g=$1; n[g]++; if($19=="yes"){r[g]++; im[g,$12]++; if($9=="hom")h[g]++; if($2=="gene_body")b[g]++;
  if($18==".")nov[g]++; if($26=="no")nohc[g]++; if($20~/athogenic/ && $20!~/onflicting/)p[g]++; if($21=="likely_pathogenic")am[g]++}}
  END{printf "gene\ttotal\trare\tnovel\tHIGH\tMOD\tLOW\tMODIF\trare_hom\trare_body\tclinvar_PLP\tAM_LP\trare_not_in_HC\n";
  for(g in n) printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",g,n[g],r[g]+0,nov[g]+0,im[g,"HIGH"]+0,im[g,"MODERATE"]+0,im[g,"LOW"]+0,im[g,"MODIFIER"]+0,h[g]+0,b[g]+0,p[g]+0,am[g]+0,nohc[g]+0}' "$TSV" | column -t
echo; echo "rows: $(($(wc -l < "$TSV")-1))  -> $TSV (server only)"
