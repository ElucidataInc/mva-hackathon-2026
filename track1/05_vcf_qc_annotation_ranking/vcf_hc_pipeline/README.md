# High-confidence variants: QC → normalisation → filtering → VEP

`vcf_hc_vep.sh` takes a **single-sample human VCF** and produces a **high-confidence, normalised,
VEP-annotated variant set**, with a QC summary that explains every variant that was removed.

```bash
conda env create -f environment.yml && conda activate vcf-hc
# one-off: VEP cache matching your VEP version and genome build (~25 GB)
vep_install -a cf -s homo_sapiens -y GRCh38 -c ~/.vep

./vcf_hc_vep.sh -i sample.vcf.gz -r GRCh38.fa -o hc_out --vep-cache ~/.vep -t 8
```

Common additions:
```bash
  --exclude-bed GRCh38-blacklist.bed      # ENCODE blacklist / low-complexity / segdups
  --include-bed exome_targets.bed --exome --bed-padding 50   # exome/panel: padded targets, no auto max-depth
  --clinvar clinvar.vcf.gz                # ClinVar, matched by exact allele after normalisation
  --resume                                # rerun after a failure without repeating an unchanged VEP run
  --sample NAME                           # if the VCF has several samples
```

## What it does

| Step | Best-practice action | Why |
|---|---|---|
| **01 prepare** | select the sample; rename contigs to match the FASTA (`1`↔`chr1`, `MT`↔`chrM`); **stop if VCF and FASTA contig lengths differ** (different builds); keep chr1–22, X, Y, M; sort and index; remove any old `CSQ` | stops build/naming mismatches from silently dropping or mis-checking variants |
| **02 normalize** | count REF mismatches vs the FASTA (warns if >1%, which suggests the wrong build); `bcftools norm -m -any` splits multi-allelics; `-f ref` left-aligns and trims indels; drop REF-mismatched records; `norm -d exact` removes duplicates; drop `*`, symbolic (`<NON_REF>`, SVs) and records where the sample has no ALT | gives one consistent representation per variant, so annotation and database matching work |
| **03 filter** | **soft filters** add a label in FILTER, so every reason is kept: caller FILTER (VQSR, DeepVariant, …), `QUAL<20`, `GQ<20`, `DP<10`, `DP>3×median` (WGS), het allele fraction outside 0.2–0.8, hom-alt fraction <0.9, GATK hard filters (QD/FS/SOR/MQ/rank-sum, only when the caller did no VQSR/CNN filtering), exclude/include BEDs. **chrM is exempt from depth and allele-balance filters** (~100× nuclear depth, heteroplasmy). Then **keep PASS** | removes low-support calls, mapping/strand artefacts, collapsed-repeat and contamination-like calls |
| **04 ClinVar** (optional) | ClinVar is renamed, split and left-aligned exactly like the sample, then matched on CHROM/POS/REF/ALT with `bcftools annotate`; P/LP variants that **failed** QC are listed separately | a pathogenic entry for a *different* allele at the same position never transfers, and pathogenic calls are never silently filtered away |
| **04 vep** | VEP offline cache with `--everything` (SYMBOL, MANE, canonical, HGVS, SIFT/PolyPhen, gnomAD/1000G AFs, rsIDs, ClinVar significance, …), `--flag_pick_allele_gene`; optional CADD, SpliceAI, REVEL, AlphaMissense plugins | functional and population context for downstream work |
| **05 tables** | `high_confidence.annotated.tsv`: one row per variant per gene (VEP's PICK = MANE/canonical transcript), with GT, GQ, DP, VAF, AD. `prioritised.tsv` (with a `priority_reason` column): rare (max population AF < 1% or absent) **and** HIGH/MODERATE impact; **or** ClinVar P/LP with ≥1 review star and no conflicts; **or** VEP-cache CLIN_SIG pathogenic | ready for filtering, review or loading elsewhere |
| **06 qc** | counts at every step, Ts/Tv, Het/HomAlt, known-variant (rsID) rate, per-filter breakdown, bcftools stats before/after → `summary.md` | confirms the final set looks like a real human germline sample |

## Outputs (`-o hc_out`)
```
final.high_confidence.vcf.gz(.tbi)     <- the deliverable: PASS, normalised, VEP-annotated
summary.md                             <- QC summary: step counts, metrics, filter reasons, sanity checks
05_tables/high_confidence.annotated.tsv
05_tables/prioritised.tsv
05_tables/clinvar_plp_failed_qc.tsv    <- ClinVar P/LP variants that FAILED a filter: confirm before discarding
03_filter/soft_filtered.vcf.gz         <- every normalised variant with its FILTER labels (audit trail)
04_vep/vep_summary.html                <- VEP's own statistics
06_qc/                                 <- step_counts.tsv, metrics.tsv, filter_breakdown.tsv, bcftools stats
logs/                                  <- run_info.txt (parameters + versions), pipeline.log, warnings.txt
```

## Tuning
All thresholds are flags (`--min-gq`, `--min-dp`, `--max-dp`, `--het-ab-min/max`, `--min-qual`,
`--gatk-hard-filters yes|no|auto`, `--max-af`). Check `summary.md` → *Why variants were filtered*:
if one filter removes a surprising share of variants (especially as the *only* reason), adjust it.

Expected for a good germline sample: Ts/Tv ≈ 2.0–2.1 (WGS) or 2.8–3.3 (exome); Het/HomAlt ≈ 1.4–2.0;
>95–99% of SNVs known in dbSNP. Use `--qc-report` to also run the full HTML QC report
(`../vcf_qc_pipeline`) on the final VCF.

## Notes
* **gVCF input**: reference blocks are dropped, but QUAL and filters only mean something after
  joint genotyping (`gatk GenotypeGVCFs` / glnexus). Genotype first if you can.
* **Somatic (tumour) VCFs** need different filters (VAF, TLOD, panel of normals). This script is for germline.
* **Large VCFs**: VEP is the slow step; use `-t`/`--vep-fork`. Only PASS variants are annotated
  unless you add `--annotate-all`. After a failure, rerun with `--resume`: VEP is skipped if its input
  and options are unchanged (checked by content fingerprint). `--cleanup` removes large intermediates.
* **chrM**: kept and annotated, but nuclear-DNA filters are not appropriate for it; for real
  mitochondrial analysis use a dedicated caller (GATK Mutect2 `--mitochondria-mode`).
* **Sex chromosomes**: heterozygous calls on chrX outside the PARs in a male sample are artefacts;
  check the Het/HomAlt on chrX if the sample's sex is known.
* No VEP cache? `--vep-gtf genes.gtf.gz` (sorted, bgzipped, tabix-indexed GTF) gives basic
  gene/consequence/HGVS annotation, without population frequencies or SIFT/PolyPhen.
