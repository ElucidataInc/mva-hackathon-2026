#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
rm -rf "$MD/manta_run"
with "$MANTA_BIN" configManta.py --bam "$FINAL_BAM" --referenceFasta "$REF" --runDir "$MD/manta_run"
with "$MANTA_BIN" python2 "$MD/manta_run/runWorkflow.py" -m local -j "${MANTA_JOBS:-$((THREADS / 2))}"
