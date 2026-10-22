#!/usr/bin/env bash
source "$(dirname "$0")/common.sh"
mkdir -p "$QC/fastqc"
run qc fastqc --threads 8 --outdir "$QC/fastqc" "$FQ_DIR"/*.fastq.gz
