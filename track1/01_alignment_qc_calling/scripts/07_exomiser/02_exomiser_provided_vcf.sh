#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
EX=$PROJECT_ROOT/tools/exomiser-cli-15.1.0
X=$RESULTS/exomiser_provided_vcf
MAP=$RESULTS/deepvariant/compare/chr_map.txt
HPO_TSV=$RESULTS/exomiser/hpo_terms.tsv
RENAMED=$X/$SAMPLE.chr.vcf.gz
IN=$X/$SAMPLE.provided.primary.pass.vcf.gz
[[ -f $MAP ]] || { echo "chr map not found: $MAP" >&2; exit 1; }
[[ -f $HPO_TSV ]] || { echo "hpo_terms.tsv not found, run 01_exomiser.sh first" >&2; exit 1; }
mkdir -p "$X"

# 1. rename contigs to chr form
run qc bcftools annotate --rename-chrs "$MAP" -Oz -o "$RENAMED" "$VCF"
run qc bcftools index -t "$RENAMED"

# 2. same PASS + primary-contig filter as 01_exomiser.sh
run qc bcftools view -f PASS -t "$(awk '$1 ~ /^chr([0-9]+|X|Y|M)$/ {printf "%s%s", s, $1; s = ","}' "$REF.fai")" \
    -Oz -o "$IN" "$RENAMED"
run qc bcftools index -t "$IN"
echo "records per contig after filtering:"; run qc bcftools index -s "$IN"
[[ $(run qc bcftools index -n "$IN") -gt 0 ]] || { echo "no variants left after filtering, check chr_map" >&2; exit 1; }
PROBAND=$(run qc bcftools query -l "$IN")

# 3. Exomiser, both term sets
for set in proband_terms all_terms; do
    hpo=$(awk -v s=$set 's == "all_terms" || $2 == "proband" {printf "%s%s", c, $1; c = ","}' "$HPO_TSV")
    mkdir -p "$X/$set"
    {
        printf -- '---\nanalysis:\n  genomeAssembly: hg38\n  vcf: %s\n  ped:\n  proband: %s\n  hpoIds: [ %s ]\n' "$IN" "$PROBAND" "$(sed "s/[^,]*/'&'/g; s/,/, /g" <<< "$hpo")"
        grep -vE '^[[:space:]]*(#|---|$)' "$EX/examples/preset-genome-analysis.yml" | sed 's/^/  /'
        printf 'outputOptions:\n  outputContributingVariantsOnly: false\n  numGenes: 0\n  outputFileName: %s\n  outputFormats: [ HTML, JSON, TSV_GENE, TSV_VARIANT, VCF ]\n' "$SAMPLE.$set"
    } > "$X/$set/$SAMPLE.$set.yml"
    run rank java -Xmx16g -Dspring.config.location="$EX/application.properties" -jar "$EX/exomiser-cli-15.1.0.jar" \
        analyse --analysis "$X/$set/$SAMPLE.$set.yml" --output-directory "$X/$set" 2>&1 | tee "$RUN_DIR/$set.exomiser.log"
    grep -q 'references/exomiser/2602' "$RUN_DIR/$set.exomiser.log" || { echo "2602 data directory not found in log" >&2; exit 1; }
done
