#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
EX=$PROJECT_ROOT/tools/exomiser-cli-15.1.0
X=$RESULTS/exomiser
DOCX=$(echo "$PROJECT_ROOT"/data/phenotype/*.docx)
IN=$X/$SAMPLE.deepvariant.primary.pass.vcf.gz
[[ -f $DOCX ]] || { echo "expected one phenotype document" >&2; exit 1; }
mkdir -p "$X"

unzip -p "$DOCX" word/document.xml | sed 's|</w:p>|\n|g; s|<[^>]*>||g' | awk '
    f && NF {print id "\t" (tolower($0) ~ /parental history/ ? "parental" : "proband"); f = 0}
    match($0, /HP:[0-9]{7}/) && $0 ~ /^[[:space:]]*HP:[0-9]{7}[[:space:]]*$/ {id = substr($0, RSTART, RLENGTH); f = 1}' > "$X/hpo_terms.tsv"
grep -q proband "$X/hpo_terms.tsv" || { echo "no proband HPO terms parsed" >&2; exit 1; }

run qc bcftools view -f PASS -t "$(awk '$1 ~ /^chr([0-9]+|X|Y|M)$/ {printf "%s%s", s, $1; s = ","}' "$REF.fai")" \
    -Oz -o "$IN" "$RESULTS/deepvariant/$SAMPLE.deepvariant.vcf.gz"
run qc bcftools index -t "$IN"
PROBAND=$(run qc bcftools query -l "$IN")

for set in proband_terms all_terms; do
    hpo=$(awk -v s=$set 's == "all_terms" || $2 == "proband" {printf "%s%s", c, $1; c = ","}' "$X/hpo_terms.tsv")
    mkdir -p "$X/$set"
    {
        printf -- '---\nanalysis:\n  genomeAssembly: hg38\n  vcf: %s\n  ped:\n  proband: %s\n  hpoIds: [ %s ]\n' "$IN" "$PROBAND" "$(sed "s/[^,]*/'&'/g; s/,/, /g" <<< "$hpo")"
        grep -vE '^[[:space:]]*(#|---|$)' "$EX/examples/preset-genome-analysis.yml" | sed 's/^/  /'
        printf 'outputOptions:\n  outputContributingVariantsOnly: false\n  numGenes: 0\n  outputFileName: %s\n  outputFormats: [ HTML, JSON, TSV_GENE, TSV_VARIANT, VCF ]\n' "$SAMPLE.$set"
    } > "$X/$set/$SAMPLE.$set.yml"
    run rank java -Xmx16g -Dspring.config.location="$EX/application.properties" -jar "$EX/exomiser-cli-15.1.0.jar" \
        analyse --analysis "$X/$set/$SAMPLE.$set.yml" --output-directory "$X/$set" > "$RUN_DIR/$set.exomiser.log" 2>&1
    grep -q 'references/exomiser/2602' "$RUN_DIR/$set.exomiser.log" || { echo "2602 data directory not found in log" >&2; exit 1; }
done
