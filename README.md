# MVA Hackathon 2026

Code and methods for [Rare Disease, Real Kid: The MVA Hackathon 2026](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026), run by Sage Bionetworks with the MVA Society, Hugging Face and BEACON.

- **Track 1, variant prediction:** rank the variants most likely to cause the child's mosaic variegated aneuploidy (MVA), using whole-genome sequencing and HPO terms.
- **Track 2, drug repurposing:** describe the variant mechanism and propose approved drugs for follow-up.

Drug candidates here are hypotheses for further study, not evidence that a medicine works, and nothing here is medical advice.

## Repository layout

| Path | Contents |
|---|---|
| [`docs/research/`](docs/research/) | Background research: disease, literature, *BUB1B* variants, biology, drugs, rules, data and tools |
| [`track1/`](track1/) | Variant prediction pipeline |
| [`track2/`](track2/) | Mechanism analysis and drug scoring |
| [`methods/`](methods/) | Methods reports for each track |
| [`submissions/`](submissions/) | Submitted files, added after the hackathon closes |
| [`env/`](env/) | Container and environment definitions with pinned versions |

## Data

The dataset is gated and cannot be redistributed, so it is not in this repository.

1. Request access on the [dataset page](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data).
2. Download it outside the repository and point the config at it (see [`track1/config/config.example.yaml`](track1/config/config.example.yaml)).

It contains one child's whole-genome sequencing (8 paired FASTQ files over 4 lanes), a VCF with index, and HPO phenotype terms: about 85 GB. The VCF is on GRCh38 (no-alt analysis set plus hs38d1 decoy, no `chr` prefix).

## Reproducing the results

Each track's README gives the exact commands. In outline:

1. Build the environment from [`env/`](env/).
2. Copy `track1/config/config.example.yaml` to `config.yaml` and set the data paths.
3. Run the Track 1 workflow to produce the ranked variant table.
4. Run the Track 2 scripts on that table to produce the scored drug list.

## Compute

Analysis runs on a single server (PNAP s3.c3.medium: Intel i9-14900K, 32 threads, 128 GB RAM (125 GiB usable), 2 × 2 TB NVMe) that holds the data, references, work files and results. A full Track 1 run from FASTQ is estimated at 7-9 hours on 48 vCPU, so expect longer on this machine. Track 2 needs little compute and can run on laptops. See [`docs/research/07_data_and_tools.md`](docs/research/07_data_and_tools.md), section 5.

## Contributing

- Never commit raw or individual-level data: FASTQ, BAM/CRAM, VCF, genotype tables or the phenotype file. `.gitignore` blocks the common formats; check `git status` before pushing.
- Keep anything that could identify the child or family out of code, comments, issues and commit messages.

## License

Code, reports and results are released under [CC BY 4.0](LICENSE), as the hackathon rules require.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

Please cite the hackathon dataset using the reference on the hackathon's [Synapse page](https://www.synapse.org/Synapse:syn76251147).
