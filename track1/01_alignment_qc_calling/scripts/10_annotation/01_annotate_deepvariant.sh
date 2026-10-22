#!/usr/bin/env bash
source "$(dirname "$0")/../01_alignment/common.sh"
PIPE=$TRACK1/05_vcf_qc_annotation_ranking/vcf_hc_pipeline
O=$RESULTS/annotation/deepvariant
mkdir -p "$(dirname "$O")"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen --no-install -e annotate)"
source "$PROJECT_ROOT/tools/py3/bin/activate"

bash "$PIPE/vcf_hc_vep.sh" -i "$RESULTS/deepvariant/$SAMPLE.deepvariant.vcf.gz" -r "$REF" -o "$O" \
    --vep-cache "$PROJECT_ROOT/references/vep/116" \
    --clinvar "$(echo "$PROJECT_ROOT"/references/clinvar/*/*.vcf.gz)" \
    --alphamissense "$(echo "$PROJECT_ROOT"/references/alphamissense/*.tsv.gz)" \
    --vep-plugins-dir "$(echo "$PROJECT_ROOT"/tools/pixi/.pixi/envs/annotate/share/ensembl-vep-*)" \
    --annotate-all -t $((THREADS / 4)) --cleanup
