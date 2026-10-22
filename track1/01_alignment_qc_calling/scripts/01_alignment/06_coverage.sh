#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
mkdir -p "$QC/coverage"
run align mosdepth -t 4 -n -Q 20 -F 3844 --by 1000000 -T 5,10,20,30 "$QC/coverage/$SAMPLE" "$FINAL_BAM"
