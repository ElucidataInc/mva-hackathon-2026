# AlphaGenome Atlas ranking

Scores the DeepVariant high-confidence variants in flanked windows around the marker genes with precomputed AlphaGenome Atlas (AVI) scores, ranks them, and keeps the rare ones.

## Setup

Python environment `tools/py3` with `alphagenome` 0.9.0 (including `alphagenome.atlas`), numpy and pandas. The Atlas client needs an AlphaGenome API key, passed with `--api-key` or the `AG_API` environment variable.

Inputs: `marker_genes.tsv` in this folder, the Ensembl 116 GTF, and `results/annotation/deepvariant/final.high_confidence.vcf.gz` from `01_alignment_qc_calling`.

## Run

From the repository root, with `PROJECT_ROOT` set:

```bash
A=track1/06_alphagenome_atlas
V=track1/01_alignment_qc_calling/results/annotation/deepvariant/final.high_confidence.vcf.gz
python3 $A/alphagenome_atlas_gene_scores.py --genes $A/marker_genes.tsv \
    --gtf $PROJECT_ROOT/references/annotation/ensembl_116/Homo_sapiens.GRCh38.116.gtf.gz \
    --variants $V --out-prefix $A/results/atlas_scores_all_hc --fetch-mode targeted --no-long
python3 $A/rank_atlas_variants.py --prefix $A/results/atlas_scores_all_hc
python3 $A/filter_rare_classify.py --ranked $A/results/atlas_scores_all_hc.ranked_input_variants.tsv --vcf $V --max-af 0.01
```

| Script | Step | Output |
|---|---|---|
| `alphagenome_atlas_gene_scores.py` | Gene windows from the GTF; Atlas scores for the variants in them | `atlas_scores_all_hc.variants_atlas_summary.tsv`, `atlas_scores_all_hc.position_alleles.tsv` |
| `rank_atlas_variants.py` | One row per input variant, ranked by AVI PHRED, then by strongest per-track effect | `atlas_scores_all_hc.ranked_input_variants.tsv` |
| `filter_rare_classify.py` | Variants with max population AF below 1% or absent, classified as coding, synonymous or non-coding | `atlas_scores_all_hc.ranked_input_variants.rare.tsv` |

The rare table is list (d) in `07_consensus`.
