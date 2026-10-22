source "$(dirname "${BASH_SOURCE[0]}")/../../01_alignment_qc_calling/scripts/01_alignment/common.sh"
eval "$("$PIXI" shell-hook --manifest-path "$PIXI_MANIFEST" --frozen --no-install -e qc)"
STEP=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
OUT=${SV_OUT:-$STEP/results}
ENVS=$(dirname "$TRACK1")/env
envbin() { echo "$(dirname "$1")/.pixi/envs/default/bin"; }
with() { local b=$1; shift; PATH="$b:$PATH" "$@"; }
CUE_BIN=$(envbin "${CUE_PIXI:-$ENVS/cue/pixi.toml}")
MANTA_BIN=$(envbin "${MANTA_PIXI:-$ENVS/manta/pixi.toml}")
TRUVARI_BIN=$(envbin "${TRUVARI_PIXI:-$ENVS/truvari/pixi.toml}")
CUE_SRC=${CUE_SRC:-$ENVS/cue/cue}
CUE_MODEL=${CUE_MODEL:-$CUE_SRC/data/models/cue.v2.pt}
SVMIN=4000

C=$OUT/cue
CR=$C/reports
MD=$OUT/manta
TV=$OUT/truvari
mkdir -p "$CR" "$MD" "$TV"
