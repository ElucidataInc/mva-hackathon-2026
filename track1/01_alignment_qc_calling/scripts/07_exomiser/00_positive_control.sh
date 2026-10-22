#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen --no-install -e qc)"
EX=$PROJECT_ROOT/tools/exomiser-cli-15.1.0
T=$RESULTS/exomiser/positive_control
CLINVAR=$(echo "$PROJECT_ROOT"/references/clinvar/*/*.vcf.gz)
EXYML=$(echo "$EX"/examples/*_38.yml)
[[ -f $EXYML && -f $CLINVAR ]] || { echo "expected one GRCh38-named example and one ClinVar VCF" >&2; exit 1; }
mkdir -p "$T"

IFS=_ read -r _ _ _ _ GENE _ POS37 _ <<< "$(basename "$EXYML" .yml)"
read -r R37 A37 < <(zcat "$EX/$(awk '$1 == "vcf:" {print $2}' "$EXYML")" | awk -v p="$POS37" '!/^#/ && $2 == p {print $4, $5}')
bcftools query -i "INFO/GENEINFO~\"^$GENE:\" && INFO/CLNSIG==\"Pathogenic\" && REF==\"$R37\" && ALT==\"$A37\"" \
    -f '%CHROM\t%POS\t%REF\t%ALT\n' "$CLINVAR" > "$T/control_variant.tsv"
[[ $(wc -l < "$T/control_variant.tsv") -eq 1 ]] || { echo "expected one ClinVar record for the control" >&2; exit 1; }
read -r CHR POS R A < "$T/control_variant.tsv"
CHR=chr${CHR#chr}

[[ $(samtools faidx "$REF" "$CHR:$POS-$((POS + ${#R} - 1))" | tail -n +2 | tr -d '\n') == "$R" ]] || { echo "reference base mismatch" >&2; exit 1; }
printf '##fileformat=VCFv4.2\n##contig=<ID=%s,length=%s>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCONTROL\n%s\t%s\t.\t%s\t%s\t100\tPASS\t.\tGT\t1/1\n' \
    "$CHR" "$(awk -v c="$CHR" '$1 == c {print $2}' "$REF.fai")" "$CHR" "$POS" "$R" "$A" > "$T/control.vcf"

sed -e 's|^  genomeAssembly: hg19$|  genomeAssembly: hg38|' -e "s|^  vcf: .*|  vcf: $T/control.vcf|" -e 's|^  outputFileName: .*|  outputFileName: control|' "$EXYML" > "$T/control.yml"
[[ $(grep -cE '^  (genomeAssembly: hg38|vcf: /.*|outputFileName: control)$' "$T/control.yml") -eq 3 ]] || { echo "yml edit failed" >&2; exit 1; }

run rank java -Xmx8g -Dspring.config.location="$EX/application.properties" -jar "$EX/exomiser-cli-15.1.0.jar" \
    analyse --analysis "$T/control.yml" --output-directory "$T" > "$RUN_DIR/exomiser.log" 2>&1
grep -q 'references/exomiser/2602' "$RUN_DIR/exomiser.log" || { echo "2602 data directory not found in log" >&2; exit 1; }
[[ $(awk -F'\t' 'NR == 2 {print $3}' "$T/control.genes.tsv") == "$GENE" ]] || { echo "control gene not ranked first" >&2; exit 1; }
