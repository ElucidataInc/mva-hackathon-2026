#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
diff <(samtools view -H "$FINAL_BAM" | awk -F'\t' '$1 == "@SQ" {print substr($2, 4) "\t" substr($3, 4)}' | sort) \
     <(cut -f1,2 "$REF.fai" | sort) > /dev/null || { echo "BAM and reference contigs differ" >&2; exit 1; }
printf 'bam: %s\nfai: %s\nchr_names: null\n' "$FINAL_BAM" "$REF.fai" > "$C/data.yaml"
printf 'model_path: %s\ngpu_ids: []\nn_cpus: %s\n' "$CUE_MODEL" "${CUE_CPUS:-10}" > "$C/model.yaml"
with "$CUE_BIN" env PYTHONPATH="$CUE_SRC" python "$CUE_SRC/engine/call.py" --data_config "$C/data.yaml" --model_config "$C/model.yaml"
