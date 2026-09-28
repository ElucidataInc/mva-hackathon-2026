# Track 1: Variant prediction

Ranks candidate causal variants from the WGS data and HPO terms.

## Layout

- `workflow/`: pipeline definition
- `scripts/`: annotation, scoring and output helpers
- `config/`: copy `config.example.yaml` to `config.yaml` (git-ignored) and set local paths

## Run

Commands will be added as the pipeline is built.

## Output

The submission CSV follows the organisers' [template](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/resolve/main/static/templates/track1_submission_template.csv):

```
proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes
```

- `proband_id` is `PROBAND01`; GRCh38 coordinates; up to 10 rows.
- Each row is one variant or one compound heterozygous pair (leave the `_2` fields blank for a single variant).
- `epcr` is in (0, 1]; `finding_type` is `primary` or `secondary`.
