# ACMG point-based scoring

Scores small variants in the marker-gene regions with ACMG/AMP points (Tavtigian et al. 2020) plus a gene-relevance weight G for the 13 panel genes below, and scores structural variants near the panel genes on four evidence components.

| G | Meaning | Genes |
|---|---|---|
| 3 | Biallelic MVA reported in several unrelated families (peer-reviewed) | BUB1B, CEP57, TRIP13 |
| 2 | Biallelic aneuploidy or chromosome-instability disease, fewer families | CENATAC, MAD2L1BP, MAD1L1, BUB1, SMC5, SLF2 |
| 1 | Preprint-only evidence, or mitotic or centrosome function without a verified MVA disease | CEP192, TRIM37, BUB3, CDC20 |
| 0 | Any other gene, or a variant in a neighbouring gene | none |

Small variants: PVS1 for null variants in genes with G >= 2 (reduced to moderate in the last exon), PM2/BS1/BA1 from the maximum population frequency, PP3/BP4 from the raw AlphaMissense score; ClinVar is flagged, not scored. Total = ACMG points + G.

Structural variants: gene content (S1), disruption of a panel gene (S2), number of genes affected (S3) and distance to the nearest panel gene (S4), with caller concordance as confidence.

## Run

Python environment from `env/sv_annotation`. Inputs: `marker_region_vep.tsv` from the first marker-gene pass of `04_alphagenome` (`01_make_marker_bed.sh`, `04_marker_rare_variants.sh`), and `combined.MANE_gnomAD_AnnotSV.tsv` and `mane_genes.bed` from `03_sv_annotation`.

```bash
python3 track1/08_acmg_scoring/score_mva_variants_v2.py \
    --snv track1/04_alphagenome/results/alphagenome/marker_region_vep.tsv \
    --sv track1/03_sv_annotation/results/annotsv_out/combined.MANE_gnomAD_AnnotSV.tsv \
    --panel-bed track1/03_sv_annotation/results/mane_genes.bed \
    --outdir track1/08_acmg_scoring/results
```

Outputs: `MVA_scored_SNVs_v2.tsv`, `MVA_scored_SVs_v2.tsv` (panel-linked or within 1 Mb) and `all_SVs_scored_v2.tsv`.

List (f) in `07_consensus` was produced by an earlier version of this scoring, which used AlphaMissense classes instead of raw scores; that version is not available.
