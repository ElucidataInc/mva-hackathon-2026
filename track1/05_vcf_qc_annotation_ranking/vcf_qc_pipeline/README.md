# VCF QC pipeline

A modular shell workflow that profiles a VCF and writes a self-contained **HTML report** with
statistics, tables and figures. It uses the standard tools (bcftools, vcftools, GATK, somalier,
peddy, hap.py/RTG, VEP/SnpEff) where they apply, and Python scripts for everything else. Only
`python3` with numpy and matplotlib is strictly required. Any module whose tool or input is
missing is skipped, and the report says why.

```bash
conda env create -f environment.yml && conda activate vcf-qc

./run_vcf_qc.sh -i cohort.vcf.gz -o qc_out -r GRCh38.fa \
    --ped family.ped --sex sex.tsv \
    --dbsnp dbsnp.vcf.gz --gnomad gnomad.sites.vcf.gz --vep-cache ~/.vep \
    --somalier-sites sites.hg38.vcf.gz --somalier-labels ancestry-labels-1kg.tsv --somalier-ref-dir 1kg-somalier/ \
    --truth HG002_GRCh38.vcf.gz --truth-bed HG002.bed --rtg-sdf GRCh38.sdf --query-sample HG002 \
    -t 8
# -> qc_out/vcf_qc_report.html
```

Somatic VCFs (Mutect2 and similar) are detected automatically. Add `--tumor T --normal N` if the
header does not declare them. Run `./run_vcf_qc.sh -h` for every option, and use `--modules profile,report`
to run only some modules.

## Modules

| # | module | tools | output dir |
|---|--------|-------|-----------|
| 00 | `prepare` | bcftools/bgzip/tabix | `00_prepare/` – BGZF + index. If the input is unsorted, a sorted copy is made for the tools (the original is still profiled, so the sort problem is reported) |
| 01 | `inspect` | bcftools | `01_inspect/` – size, md5, header, provenance lines, samples, records per contig |
| 02 | `integrity` | `bcftools norm -c w`, GATK `ValidateVariants` | `02_integrity/` – REF mismatches, left-align/split counts, validation log |
| 03 | `stats` | `bcftools stats -s - [-F ref]`, `plot-vcfstats` | `03_bcftools_stats/` – all-sample and PASS-only stats. Retries automatically after bad AC/AN or REF mismatches |
| 04 | `profile` | `scripts/vcf_profile.py` | `04_profile/` – profile.json, per_chrom.tsv, per_sample.tsv, density_bins.tsv, gt_matrix.npz |
| 05 | `sampleqc` | vcftools, `+guess-ploidy`, `+check-ploidy`, `roh`, `+mendelian2`, somalier, peddy, `scripts/sample_qc.py` | `05_sample_qc/` |
| 06 | `annotation` | `bcftools annotate` (dbSNP IDs, gnomAD AF), VEP or SnpEff, `scripts/annotation_summary.py` | `06_annotation/` |
| 07 | `concordance` | hap.py, RTG vcfeval | `07_concordance/` |
| 08 | `gatk` | GATK `CollectVariantCallingMetrics` | `08_gatk_metrics/` |
| 09 | `report` | `scripts/parse_tool_outputs.py`, `scripts/build_report.py` | `vcf_qc_report.html` |

Every run writes `logs/module_status.tsv` (OK / SKIPPED / FAILED plus the reason), `logs/tool_versions.txt`,
`logs/params.txt` and one log per module.

## Where each metric comes from

| Metric | Source |
|---|---|
| Size, compression, index, checksum | inspect, profile |
| How it was produced (caller, command lines, build, reference) | header parsing (callers detected: GATK HC/GenotypeGVCFs/Mutect2/VQSR, DeepVariant, GLnexus, DRAGEN, Strelka2, FreeBayes, bcftools, Clair3, Sentieon, …). Build inferred from contig lengths |
| Sort order, duplicates | profile (positions going backwards, contigs that reappear, exact duplicate records, records sharing a position) |
| REF concordance | profile (indexed FASTA reader) + `bcftools norm -c w` |
| Filters, QUAL, DP, GATK INFO annotations (QD/FS/SOR/MQ/…) | profile + bcftools stats |
| Variant kinds and classes (SNV/MNV/ins/del/complex/symbolic SV/BND), Ts/Tv, substitution spectrum | profile + bcftools stats |
| Biallelic vs multiallelic sites | profile |
| Indel length distribution | profile + bcftools stats IDD |
| Counts, types, het/hom-alt, missingness per chromosome | profile → per_chrom.tsv |
| Density and depth chromosome maps (Manhattan-style + heat-map ideogram) | profile → density_bins.tsv (`--bin-size`) |
| Hom-ref/het/hom-alt counts, phasing (phased hets, PS blocks) | profile, bcftools stats PSC |
| Genotype missingness plots (per sample, site, chromosome) | profile + vcftools `--missing-indv/--missing-site` |
| Known vs novel (rsID / DB flag / gnomAD AF), Ts/Tv known vs novel | profile, annotation module |
| Functional annotation (consequence, impact, genes, ClinVar, gnomAD AF bins) | existing CSQ/ANN/BCSQ/CLNSIG in the VCF, or VEP/SnpEff in module 06 |
| Per-sample Ts/Tv, Het/HomAlt, missingness, mean DP/GQ, singletons, private variants | profile, bcftools stats PSC, vcftools |
| Sex check | chrX non-PAR het ratio + chrY call rate vs declared sex (sample_qc.py), `+guess-ploidy`, peddy, somalier |
| Relatedness / duplicates | KING-robust kinship + IBS0 (sample_qc.py), somalier relate, peddy |
| Ancestry | within-cohort PCA (sample_qc.py). Projection onto 1000G via `somalier ancestry` |
| Contamination | het-rate outliers and het allele-balance skew (VCF-level proxy; confirm with VerifyBamID2 on the BAM) |
| Runs of homozygosity (consanguinity, UPD hints) | `bcftools roh`, or a fallback scan of the genotype matrix |
| Mendelian errors | profile (trios from PED) + `bcftools +mendelian2` |
| HWE, inbreeding F | vcftools `--hardy`, `--het` |
| Somatic: VAF, TLOD/NLOD/NALOD/GERMQ/POPAF, tumor/normal depth, F1R2/F2R1 orientation bias, SB strand bias, PON/germline/artefact filters | profile (somatic mode) |
| Truth concordance | hap.py / RTG vcfeval |

## Tests

```bash
tests/run_tests.sh            # synthetic 12-sample cohort + Mutect2-style tumor/normal
```

The synthetic cohort contains known problems the report should catch: a duplicate sample pair,
a trio, a contaminated sample, a sample whose declared sex is wrong, a sample with a long ROH,
a high-missingness sample, REF mismatches, a duplicate record and an unsorted record.
`example_reports/` holds the two reports generated from it.

## Notes and limits
* The profiler reads the file in a single pass with pure Python (about 20k records/s for a few samples).
  For very large cohorts, use `--max-records` for a quick look, or run it on a subset or region
  (`bcftools view -r`). For biobank-scale data, use Hail.
* Kinship, PCA and the fallback ROH scan use a random sample of autosomal PASS biallelic SNVs
  (`--max-gt-sites`, default 100k).
* PAR coordinates are GRCh38 or GRCh37, chosen from the inferred build.
* gVCFs are recognised (`<NON_REF>` reference blocks are counted separately). It is still better to
  QC the genotyped VCF.
