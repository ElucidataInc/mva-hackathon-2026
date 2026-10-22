# MVA Hackathon 2026

Code and methods for [Rare Disease, Real Kid: The MVA Hackathon 2026](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026), run by Sage Bionetworks with the MVA Society, Hugging Face and BEACON.

- **Track 1, variant prediction:** rank the variants most likely to cause the child's mosaic variegated aneuploidy (MVA), using whole-genome sequencing and HPO terms.
- **Track 2, drug repurposing:** describe the variant mechanism and propose approved drugs for follow-up.

Drug candidates here are hypotheses for further study, not evidence that a medicine works, and nothing here is medical advice.

## Repository layout

| Path | Contents |
|---|---|
| [`track1/`](track1/) | Variant prediction: shared preprocessing and six prioritisation approaches |
| [`track2/`](track2/) | Mechanism analysis and drug scoring |
| [`env/`](env/) | pixi workspaces with lock files: `pixi` (shared tools), `cue`, `manta`, `truvari`, `sv_annotation` |
| [`docs/`](docs/) | [`MVA_full_flowchart.html`](docs/MVA_full_flowchart.html) (every Track 1 step with inputs, outputs and run times) and background research |
| [`methods/`](methods/) | Methods reports for each track |
| [`submissions/`](submissions/) | Track 1 predictions file as submitted |

## Data

The dataset is gated and cannot be redistributed, so it is not in this repository.

1. Request access on the [dataset page](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data).
2. Place it under the project root set as `PROJECT_ROOT` in `track1/01_alignment_qc_calling/scripts/01_alignment/common.sh` (and `ROOT` in `track1/04_alphagenome/scripts/alphagenome/*.sh`): `data/fastq/`, `data/variants/`, `data/phenotype/`.

It contains one child's whole-genome sequencing (8 paired FASTQ files over 4 lanes), a VCF with index, and HPO phenotype terms: about 85 GB. The VCF is on GRCh38 (no-alt analysis set plus hs38d1 decoy, with GRC exclusion masking v2, no `chr` prefix).

References go under `references/` in the same root: the NCBI GRCh38 no-alt plus hs38d1 analysis set and the GIAB masked set (`genome/`), Ensembl 116 (`annotation/ensembl_116/`), ClinVar 2026-09-23, AlphaMissense, ReMM 0.4, Exomiser data 2602, the VEP 116 cache, UCSC hg38 centromeres, and the structural variant references from `track1/03_sv_annotation/scripts/00_references.sh`.

## Reproducing the results

1. Install pixi 0.81.0. The shared workspace in `env/pixi` is installed at `$PROJECT_ROOT/tools/pixi`; the step workspaces in `env/<name>` are installed in place with `pixi install --manifest-path env/<name>/pixi.toml`.
2. Install the tools outside pixi: GATK 4.7.0.0 with OpenJDK 17, Exomiser CLI 15.1.0, DeepVariant 1.10.0 (Docker image `google/deepvariant:1.10.0`) and a Python 3.12 environment at `tools/py3` with `alphagenome` 0.9.0.
3. Run the Track 1 steps in the order given in [`track1/README.md`](track1/README.md).

## Compute

All analyses ran on a single server: Intel Core i9-14900K (32 logical CPUs), 125 GiB RAM, 1.8 TB disk, Ubuntu 24.04.5 LTS, no GPU. Run one after another, Track 1 takes about 25 to 27 hours of compute; the flowchart gives the time of each step. Track 2 needs little compute.

## Contributing

- Never commit raw or individual-level data: FASTQ, BAM/CRAM, VCF, genotype tables or the phenotype file. `.gitignore` blocks the common formats; check `git status` before pushing.
- Keep anything that could identify the child or family out of code, comments, issues and commit messages.

## License

Code, reports and results are released under [CC BY 4.0](LICENSE), as the hackathon rules require.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

Please cite the hackathon dataset using the reference on the hackathon's [Synapse page](https://www.synapse.org/Synapse:syn76251147).
