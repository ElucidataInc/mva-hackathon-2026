# Track 1: Variant prediction

Six prioritisation approaches applied to one whole-genome sample, built on shared alignment, variant calling, phasing and annotation. [`docs/MVA_full_flowchart.html`](../docs/MVA_full_flowchart.html) shows every step with its inputs, outputs and run time.

## Steps

| Folder | Contents |
|---|---|
| [`01_alignment_qc_calling`](01_alignment_qc_calling/) | Reference build, alignment, QC, DeepVariant calling, comparison with the provided calls, phasing, runs of homozygosity, marker-gene coverage, Exomiser on both call sets, candidate review |
| [`02_sv_calling`](02_sv_calling/) | Cue and Manta structural variant calls, read-depth test, truvari confirmation |
| [`03_sv_annotation`](03_sv_annotation/) | MANE, gnomAD-SV and AnnotSV annotation of structural variants of 4 kb or more |
| [`04_alphagenome`](04_alphagenome/) | AlphaGenome scoring of rare variants within 50 kb of the 11 marker genes |
| [`05_vcf_qc_annotation_ranking`](05_vcf_qc_annotation_ranking/) | VCF QC, normalisation and VEP annotation, tiered rule-based ranking |
| [`06_alphagenome_atlas`](06_alphagenome_atlas/) | Rare variants ranked by AlphaGenome Atlas (AVI) score |
| [`07_consensus`](07_consensus/) | Comparison of the six ranked lists and evidence flags for review |
| [`08_acmg_scoring`](08_acmg_scoring/) | ACMG point-based scoring of marker-gene variants and nearby structural variants |

The six ranked lists come from `01_alignment_qc_calling` (Exomiser on each call set), `04_alphagenome`, `05_vcf_qc_annotation_ranking`, `06_alphagenome_atlas` and `08_acmg_scoring`; `07_consensus` compares them.

## Run order

Scripts in `01_alignment_qc_calling`, `02_sv_calling` and `03_sv_annotation` take a run folder as their first argument; `00_reference/build_masked_reference.sh` takes none.

1. `01_alignment_qc_calling`: `00_reference` to `04_variant_calling`.
2. `05_vcf_qc_annotation_ranking`: `vcf_qc_pipeline`, then `vcf_hc_pipeline` on the provided VCF, with output in `05_vcf_qc_annotation_ranking/results/hc_out`.
3. `01_alignment_qc_calling`: `10_annotation`, `03_phasing` with `provided_hc` and with `deepvariant_hc`, then `05_homozygosity` to `09_marker_genes`.
4. `02_sv_calling`, then `03_sv_annotation`.
5. `04_alphagenome`: `01_make_marker_bed_v2.sh`, `02_variant_density_v2.sh`, `03_roh.sh`, `04_marker_rare_variants_v2.sh`, then `06` to `09` (Python, in `tools/py3` with `ALPHAGENOME_API_KEY` set). `01_make_marker_bed.sh` and `04_marker_rare_variants.sh` build the first-pass table used by `08_acmg_scoring`.
6. `06_alphagenome_atlas`: commands in [its README](06_alphagenome_atlas/README.md).
7. `05_vcf_qc_annotation_ranking`: ranking, run from the repository root with `PROJECT_ROOT` set (command below).
8. `08_acmg_scoring`: commands in [its README](08_acmg_scoring/README.md).
9. `07_consensus`: comparison of the six lists.

```bash
T=track1
S=WGS_EX2312012
bash $T/05_vcf_qc_annotation_ranking/vcf_rank_pipeline/run_ranking.sh \
    --vcf1 $T/01_alignment_qc_calling/results/phasing/provided_hc/$S.phased.vcf.gz \
    --vcf2 $T/01_alignment_qc_calling/results/phasing/deepvariant_hc/$S.phased.vcf.gz \
    --combine union \
    --sv $T/03_sv_annotation/results/annotsv_out/combined.MANE_gnomAD_AnnotSV.tsv \
    --gene-gtf $PROJECT_ROOT/references/sv_annotation/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
    --bam $T/01_alignment_qc_calling/results/alignment/$S.md.bam \
    --exomiser-genes $T/01_alignment_qc_calling/results/exomiser_provided_vcf/proband_terms/$S.proband_terms.genes.tsv \
    --exomiser-genes $T/01_alignment_qc_calling/results/exomiser/proband_terms/$S.proband_terms.genes.tsv \
    --alphagenome $T/04_alphagenome/results/alphagenome/v2_team11/alphagenome_model/FINAL_alphagenome_variant_table.tsv \
    -i $T/05_vcf_qc_annotation_ranking/results/hc_out \
    -o $T/05_vcf_qc_annotation_ranking/results/rank_out_v38 -t 8
```

## Submission

The submission file follows the organisers' [template](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/resolve/main/static/templates/track1_submission_template.csv): `PROBAND01`, GRCh38 coordinates, up to 10 rows, one variant or one compound heterozygous pair per row, `epcr` in (0, 1], `finding_type` primary or secondary.
