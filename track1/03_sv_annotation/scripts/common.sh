source "$(dirname "${BASH_SOURCE[0]}")/../../01_alignment_qc_calling/scripts/01_alignment/common.sh"
ANNOT_PIXI=${ANNOT_PIXI:-$(dirname "$TRACK1")/env/sv_annotation/pixi.toml}
eval "$("$PIXI" shell-hook --manifest-path "$ANNOT_PIXI" --frozen)"
STEP=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
IN=${SV_OUT:-$TRACK1/02_sv_calling/results}
OUT=${SVANN_OUT:-$STEP/results}
SVREF=${SVREF:-$PROJECT_ROOT/references/sv_annotation}
SVMIN=4000
AS=$OUT/annotsv_out
mkdir -p "$OUT"
