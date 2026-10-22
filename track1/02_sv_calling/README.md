# Structural variant calling

Calls large structural variants (SVs) from the final BAM with Cue and Manta, tests Cue deletions and duplications against read depth, and keeps the Cue calls that Manta confirms.

## Setup

```bash
for e in cue manta truvari; do pixi install --manifest-path env/$e/pixi.toml; done
git clone https://github.com/PopicLab/cue.git env/cue/cue
git -C env/cue/cue checkout 551cb72f7b9f4177c9ac71743f8ab7cf3d7f28dc
env/cue/.pixi/envs/default/bin/pip install --no-cache-dir -r env/cue/cue/install/requirements.txt
wget --directory-prefix=env/cue/cue/data/models/ https://storage.googleapis.com/cue-models/latest/cue.v2.pt
```

The model file `cue.v2.pt` has MD5 `3b5f0a5af6079ed270368fd6f540ad3d`. `env/cue/pip_freeze.txt` lists the Python packages of the environment used for the analysis.

Input: the final BAM and the reference from `01_alignment_qc_calling`.

## Run

Each script takes a run folder as its first argument and writes to `results/`.

| Script | Step | Output |
|---|---|---|
| `01_cue_call.sh` | Cue on all chromosomes, CPU only | `cue/reports/svs.vcf` |
| `02_cue_depth_filter.sh` | Mean depth per chromosome from `samtools idxstats`; depth ratio of each DEL and DUP from `samtools bedcov` at MAPQ >= 20; DEL kept at ratio <= 0.7, DUP at >= 1.3; inversions split out | `cue/reports/svs.depthvalidated.vcf.gz`, `cue/cue.inv.vcf.gz` |
| `03_manta_call.sh` | Manta 1.6.0 germline workflow | `manta/manta_run/results/variants/diploidSV.vcf.gz` |
| `04_manta_filter.sh` | PASS calls; breakend pairs converted to inversions | `manta/manta.pass.vcf.gz`, `manta/manta.pass.inv.vcf.gz` |
| `05_truvari.sh` | truvari bench of Cue against Manta: DEL and DUP with type ignored, INV with type matched; 4 kb to 3 Mb, 50% reciprocal overlap, 1 kb reference distance | `cue.manta_confirmed.sorted.vcf.gz`, `igv_review.tsv` |

`igv_review.tsv` lists the depth-validated Cue DEL and DUP calls without a Manta match, for review in IGV.

Optional settings: `CUE_CPUS` (default 10), `MANTA_JOBS` (default half the cores), `SV_OUT` (output folder).

Steps 02, 04 and 05, rerun on the original Cue and Manta calls, reproduce the original files exactly.
