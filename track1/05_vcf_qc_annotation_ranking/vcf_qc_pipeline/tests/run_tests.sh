#!/usr/bin/env bash
# Build the synthetic test data set and run the pipeline on the germline cohort and the
# somatic tumor/normal VCF. Usage: tests/run_tests.sh [workdir]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="${1:-${HERE}/test_run}"
mkdir -p "$W"
python3 "${HERE}/make_test_data.py" "$W/data"
D="$W/data"
if command -v bcftools >/dev/null; then
    bgzip -f "$D/cohort.vcf"            # deliberately unsorted -> pipeline must create a sorted copy
    bcftools view -Oz -o "$D/known_sites.vcf.gz" "$D/known_sites.vcf" && bcftools index -t -f "$D/known_sites.vcf.gz"
    bgzip -f "$D/somatic.vcf"
    KNOWN=(--dbsnp "$D/known_sites.vcf.gz" --gnomad "$D/known_sites.vcf.gz")
else
    gzip -f "$D/cohort.vcf"; gzip -f "$D/somatic.vcf"; KNOWN=()
fi
"${HERE}/../run_vcf_qc.sh" -i "$D/cohort.vcf.gz" -o "$W/germline" -r "$D/ref.fa" --ped "$D/family.ped" \
    --sex "$D/sex.tsv" ${KNOWN[@]+"${KNOWN[@]}"} --bin-size 200000 --title "Synthetic germline cohort – QC"
"${HERE}/../run_vcf_qc.sh" -i "$D/somatic.vcf.gz" -o "$W/somatic" -r "$D/ref.fa" --bin-size 200000 \
    --title "Synthetic tumor/normal (Mutect2-style) – QC"
echo "Reports:"; ls -1 "$W"/*/vcf_qc_report.html
