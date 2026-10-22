# Structural variant annotation

Merges the Manta and Cue calls of 4 kb or more and annotates them with MANE Select genes, gnomAD-SV v4.1 frequencies and AnnotSV 3.5. The final table is the `--sv` input of `05_vcf_qc_annotation_ranking`.

## Setup

```bash
pixi install --manifest-path env/sv_annotation/pixi.toml
bash track1/03_sv_annotation/scripts/00_references.sh <run_folder>
```

`00_references.sh` downloads the MANE v1.5 GTF, the SVAFotate gnomAD-SV v4.1 BED and the AnnotSV 3.5 annotations (about 5.8 GB to download) into `references/sv_annotation`. Set `SVREF` to use another folder.

Input: `manta/manta.pass.vcf.gz` and `cue/reports/svs.vcf` from `02_sv_calling`.

## Run

Each script takes a run folder as its first argument and writes to `results/`.

| Script | Step | Output |
|---|---|---|
| `01_merge_mane.sh` | Manta PASS and Cue calls with absolute SVLEN >= 4 kb merged into one VCF (IDUP recorded as DUP); MANE Select genes overlapped by each SV | `combined.mane.vcf.gz` |
| `02_gnomad_annotsv.sh` | SVAFotate best gnomAD-SV match at 50% overlap; AnnotSV full annotation; MANE and gnomAD fields extracted into columns (`extract_all_annotations.py`) | `annotsv_out/combined.MANE_gnomAD_AnnotSV.tsv` |

Rerun on the original calls, both scripts reproduce the original files exactly.
