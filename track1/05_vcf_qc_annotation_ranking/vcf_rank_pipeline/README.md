# Candidate ranking module (v3.8)

This module looks for the genetic cause of a **recessive** disorder in one patient. The default focus is
**mosaic variegated aneuploidy (MVA)**. A recessive disorder needs **two damaging hits in the same gene**, one on
each copy. The module collects all hits from:
- small variants called by two variant callers;
- structural variants (SVs);
- read depth from the BAM;
- AlphaGenome predictions.

It then decides, for each gene, whether a two-hit genotype is plausible, ranks the genes, and writes a short list
of genes to confirm in the lab.

**Scores order the candidates for review. They are not ACMG classifications and are not a diagnosis.**

```bash
bash run_ranking.sh \
    --vcf1 provided_hc/WGS_EX2312012.phased.vcf.gz --vcf2 deepvariant_hc/WGS_EX2312012.phased.vcf.gz --combine union \
    --sv combined.MANE_gnomAD_AnnotSV.tsv --gene-gtf MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
    --bam WGS_EX2312012.md.bam \
    --exomiser-genes <provided>.proband_terms.genes.tsv --exomiser-genes <deepvariant>.proband_terms.genes.tsv \
    --alphagenome FINAL_alphagenome_variant_table.tsv \
    -i 01_VCF_Annotate/hc_out -o rank_out -t 8
```

**Environment:** run it in the pixi `align` environment. That environment has bcftools 1.19 (with `+split-vep`),
samtools 1.24 and mosdepth.

---

## 1. The ranking steps and why each one matters

| # | Step | What it does | Why it matters |
|---|---|---|---|
| 0 | **Re-annotate** (optional, `reannotate_plugins.sh`) | Rebuilds the VEP annotation of each phased VCF with REVEL, SpliceAI and CADD. Genotypes, phase and ClinVar are kept. | Without these scores, missense and splicing evidence is weaker (see §2). |
| 1 | **Variant table** (`01_table`) | Reads both VCFs and keeps, per variant and gene, the VEP preferred (MANE) transcript. Removes common variants (population AF ≥1%, unless ClinVar pathogenic), variants with no gene, and variants with no qualifying evidence. Adds ClinVar P/LP calls that **failed** QC, flagged. | Most of the genome is common or harmless. Removing it early leaves a few hundred candidates and keeps memory low (~20 MB for a whole genome). Failed-QC pathogenic calls are kept because a known pathogenic variant should be confirmed, not discarded. |
| 1b | **Merge the two callers** | Matches variants on position, alleles and gene. `union` keeps variants from either caller, labelled `both`, `A` or `B`; `consensus` keeps only shared ones. An indel written differently by the two callers counts as one variant. | A variant found by both callers is more likely real. Counting one indel twice would create a **false compound heterozygote**. |
| 2 | **Coverage** (`02_coverage`, needs `--bam`) | Exon depth of the 10 MVA genes vs the median exon. Median depth of 1-Mb windows per chromosome, ignoring gap windows. | Finds exon deletions or duplications the SV callers missed. Screens for **mosaic aneuploidy**, the cellular hallmark of MVA, directly in the patient's DNA. |
| 3 | **Structural variants** (`03_sv`) | Turns AnnotSV output into one candidate per SV and gene. Assigns an effect (coding loss, in-frame, splice, regulatory, intronic). Removes common SVs. Rates breakpoint reliability. **Measures depth inside every deletion and duplication** in the BAM to get its real genotype. | An SV can be one of the two hits. AnnotSV writes a placeholder genotype (`0/1`), so without depth a **homozygous deletion would be missed**. Depth also rejects SV calls with normal coverage (likely false). |
| 4 | **Evidence class per allele** | Gives each variant one class and score (§2, table A). AlphaGenome can raise the class of a non-coding variant. | Puts very different evidence (ClinVar, loss of function, missense predictors, SVs) on one scale. |
| 5 | **Genotype per gene** | Homozygous or hemizygous → two hits. Two heterozygous hits → checks **phase** (same or opposite copies) within each caller's phase sets. Chooses the best pair that is not on the same copy. A SNV that looks homozygous inside a heterozygous deletion → hemizygous, on the opposite copy to the deletion. | Two hits on the **same** copy (cis) leave one working copy, so they don't cause a recessive disorder. Phase is what separates a real compound heterozygote from a false one. |
| 6 | **Tier and order** | Tiers 1–6 (§2, table B). Within a tier: genes supported by both callers first, then Exomiser phenotype score, then evidence score. | MVA genes with strong two-hit evidence come first. Exomiser brings in the patient's HPO phenotype. |
| 7 | **Targets** (`05_targets`) | Lists every variant (both callers plus the filtered calls) and SV in the target genes, with an IGV BED file. | The reviewer sees everything in each target gene, not only what passed the scoring. A missed second hit is often among the filtered or low-scored calls. |

---

## 2. What evidence the ranking uses

### Evidence sources, and their status for this case (WGS_EX2312012)

The numbers come from the last full run (log of 15 Oct 2026, module v3.5). They will change after the re-run with
v3.7.

| Source | What it tells us | Strength | Status for this case |
|---|---|---|---|
| ClinVar (exact allele, ≥1 star, not conflicting) | Variant already reported as pathogenic | strongest | ✅ in both VCFs. 1 failed-QC P/LP call added from `hc_out` (unphased). |
| VEP consequence and impact (MANE transcript) | Loss of function: stop, frameshift, canonical splice, start loss | strong | ✅ |
| AlphaMissense | Missense pathogenicity | moderate | ✅ |
| SIFT + PolyPhen | Missense damage, older predictors | supporting at most | ✅ used when REVEL is missing |
| REVEL | Missense pathogenicity, ClinGen-calibrated grades | strong–supporting | ❌ missing. Add it with `reannotate_plugins.sh`. |
| SpliceAI | Splice effects outside the canonical ±2 bp | moderate | ❌ missing. **This is the main gap** for finding a hidden second hit. |
| CADD | Non-coding harm score | weak | ❌ missing (small effect on the ranking) |
| Population frequency (gnomAD, 1000 Genomes via VEP `MAX_AF`) | Common = unlikely to cause a rare disorder | filter | ✅ |
| Two callers (provided + DeepVariant) | Call reliability | confidence | ✅ 470 variant-gene rows in both, 63 only in A, 14 only in B |
| Phasing (WhatsHap, read-based) | Same or opposite copy | essential for compound hets | ✅ 85% (A) and 83% (B) of heterozygous calls phased. Pairs across phase blocks stay "unphased". |
| AnnotSV (Manta + cue) | SVs, their effect, frequency, known pathogenic CNVs | strong for coding deletions | ⚠️ the file has **only full lines** (no per-gene effect columns). Effects are computed from gene models: use `--gene-gtf` (MANE). Without it, only the 10 MVA genes are evaluated. |
| BAM read depth | Exon CNVs, SV genotypes, mosaic aneuploidy | strong for deletions | ✅ (42× genome mean). In v3.5 the exon and chromosome screens had a reference artefact. Fixed in v3.6: re-run. |
| Exomiser (proband HPO terms only) | Phenotype match | orders genes within a tier | ✅ 238 genes, max of the 2 files |
| AlphaGenome (46 variants in marker-gene regions) | Splicing, expression and regulatory effect of non-coding variants | supporting (not clinically calibrated) | ✅ 8 of 46 with a class, 2 failed QC, 13 without signal |
| Parents, RNA, long reads | Phase, de novo status, splice and expression proof | confirmatory | ❌ not available yet. These are the confirmation steps. |

### A. Evidence classes (one per variant)

| Class | Score | Rule |
|---|---:|---|
| CLINVAR_PLP | 8 | ClinVar pathogenic/likely pathogenic, ≥1 star, not conflicting |
| LOF / SV_LOF | 6 | stop, frameshift, canonical splice, start loss / SV deleting coding sequence (frameshift or ≥50% of CDS) or the whole gene |
| MISSENSE_STRONG | 5 | REVEL ≥ 0.932 |
| MISSENSE_MODERATE | 4 | REVEL ≥ 0.773, or AlphaMissense likely pathogenic |
| SPLICE_PREDICTED / SV_LOF_POSSIBLE / AG_SPLICE (both AlphaGenome methods) | 4 | SpliceAI ≥ 0.2 / in-frame exon deletion, partial duplication, breakpoint in the gene / AlphaGenome splice change confirmed by both methods |
| MISSENSE_SUPPORTING | 3 | REVEL ≥ 0.644; **without REVEL:** SIFT deleterious **and** PolyPhen probably damaging |
| CNV_DEPTH / SV_SPLICE / AG_SPLICE / AG_EXPRESSION | 3 | exon depth CNV with no SV call / intronic SV ≤50 bp from a splice site / AlphaGenome splice change (one method) / reduced expression |
| MISSENSE_VUS / SV_REGULATORY / AG_REGULATORY | 2 | other missense / SV in promoter, UTR or a regulatory element / AlphaGenome regulatory change |
| MISSENSE_PRED_BENIGN | 1 | **without REVEL:** AlphaMissense likely benign (kept visible, can't be a strong hit) |
| SPLICE_REGION / NONCODING_CADD / SV_INTRONIC / SV_DOSAGE | 1 | splice region ±3–8 bp / CADD ≥20 / intronic SV / whole-gene duplication |

**Removed:**
- population AF ≥1% (unless ClinVar P/LP)
- REVEL ≤0.290
- no qualifying consequence
- failed QC (unless ClinVar P/LP)

**SV score adjustments:** −1 for low breakpoint confidence, for a benign-CNV overlap of all its genes, and for no depth support (≥1 kb). An allele counts as **strong** only with a strong class **and** a score ≥4.

### B. Tiers (one per gene)

| Tier | Rule | What it means |
|---|---|---|
| 1 | MVA/panel gene, biallelic, ≥1 ClinVar P/LP allele, other allele ≥3 | Best candidate: confirm first |
| 2 | Established MVA/panel gene, biallelic, both alleles strong (≥4) | Strong candidate |
| 3 | MVA/panel gene, biallelic with weaker evidence, or an emerging MVA gene | Plausible: one allele needs more support |
| 4 | MVA/panel gene, one strong allele only | **Look for the second hit** (SV, depth, splicing, RNA) |
| 5 | Other gene, biallelic, both alleles strong | Discovery candidate outside the MVA list |
| 6 | Everything else | Not followed up |

**MVA genes:**
- established: BUB1B, CEP57, TRIP13, CENATAC, MAD1L1, SLF2, SMC5
- emerging: BUB1, MAD2L1BP, TRIM37

Add your own panel with `--genes`.

---

## 3. What needs to be reviewed

### Case-level checks (before reading the gene list)

| Check | Where | Action |
|---|---|---|
| All steps ran, nothing skipped | `logs/ranking.log`, `logs/warnings.txt` | Every warning should be understood. `02 coverage` must say `(mosdepth)` or `(samtools)`, never "skipped". |
| SV effects are per gene | `03_sv/sv_annotate.log` | "only a full line" means use `--gene-gtf` or an AnnotSV file with split lines. The count on the "gene without a model" line should be small. |
| Mosaic aneuploidy screen | `02_coverage/karyotype_screen.tsv` | A real mosaic gain or loss changes depth consistently across many windows (`windows_used`). One flagged autosome is a key MVA finding: compare with the karyotype or cytogenetics. Check the sex chromosomes against the recorded sex. |
| Exon CNV screen | `02_coverage/exon_depth.tsv`, `summary.md` | Several adjacent exons with the same change are more convincing than a single exon. GC-rich first exons are often low. |
| Missing evidence | log `note: no REVEL ...` | Know which classes couldn't be assigned (REVEL, SpliceAI, CADD), and add them if possible. |
| Phenotype input | `--exomiser-genes` | Use the **proband_terms** files. Use `all_terms` only as a sensitivity check. |

### Per target gene (`targets.tsv`, then IGV with `targets.bed`)

1. **Read the `note` column first.** Each note is a specific thing to confirm:

   | Note | Meaning | How to confirm |
   |---|---|---|
   | `confirm in trans` / `unphased` / `SV-small variant phase unknown` | Two hits, but it is unknown whether they are on opposite copies | Parental testing (one allele from each parent), or long reads |
   | `found in only one VCF` | Only one caller saw the allele | IGV: reads, strand, mapping quality; Sanger |
   | `genotype differs between the two VCFs` | Callers disagree (e.g. 0/1 vs 1/1) | IGV allele balance |
   | `FAILED QC` | Known pathogenic call that failed a filter | IGV, then Sanger |
   | `homozygous deletion called from read depth` | Depth ≈0 inside the SV | MLPA / qPCR, breakpoint PCR |
   | `SV not supported by read depth` | SV call with normal depth: probably false | IGV split and discordant reads; usually discard |
   | `SV breakpoint in a gap/blacklist/segdup region` | Unreliable SV region | IGV, long reads |
   | `depth-only CNV` | Exon depth change with no SV call | MLPA / qPCR |
   | `allele supported only by AlphaGenome prediction` | Non-coding or splice hit from a model | RNA (RT-PCR or RNA-seq from blood or fibroblasts), reporter assay |
   | `second hit?` (tier 4) | Only one strong allele | Look at `05_targets/target_gene_variants.tsv` and `target_gene_svs.tsv` for weaker or filtered calls; RNA-seq |
   | `dosage-sensitive gene` | One SV copy may be enough | Consider a dominant model |
   | `mitochondrial` | Heteroplasmy, not biallelic logic | Interpret by allele fraction |

2. **Check each allele in IGV:** read support, allele balance (~50% heterozygous, ~100% homozygous; mosaic if clearly lower), strand bias, repeats, and the phase (reads linking the two variants).
3. **Check the variant against the literature and ClinVar:** HGVS on the MANE transcript (`hgvsc`, `hgvsp` in `ranked_variants.tsv`), and the disease mechanism of the gene (loss of function for MVA).
4. **Check the phenotype fit:** the Exomiser score and OMIM/GenCC disease (`gene_info` for SV genes). Ask whether the patient's features match.
5. **Decide the follow-up:** phase (parents or long reads), CNV (MLPA/qPCR), splicing and expression (RNA), then ACMG classification of each confirmed allele by a clinical scientist.

### Known limits

- **Phase only within one caller's phase blocks.** Hits far apart, or one SV and one small variant, stay "unphased" until confirmed.
- **Mosaic variants:** a low allele fraction is not modelled as a separate genotype. Check the allele balance in IGV.
- **Predictions** (AlphaMissense, SIFT/PolyPhen, AlphaGenome) support, but don't prove, an effect.
- **De novo variants and dominant models** are not the focus. Dosage-sensitive genes are only noted.

---

## 4. Reference

### Files

| File | Role |
|---|---|
| `run_ranking.sh` | driver (steps above) |
| `reannotate_plugins.sh` | optional, before ranking: rebuilds VEP annotation with REVEL / SpliceAI / CADD / AlphaMissense (GT, PS, ClinVar kept) |
| `prefilter_table.py` | streams each table and keeps possible candidates |
| `merge_tables.py` | union / consensus merge of the two VCFs |
| `append_aligned.py` | adds failed-QC ClinVar rows, matching columns by name |
| `window_medians.py` | per-chromosome median of 1-Mb window depths (gaps ignored) |
| `sv_annotate.py` | SV candidates per gene (AnnotSV TSV or SV VCF), depth genotypes from the BAM |
| `rank_variants.py` | evidence classes, genotypes, tiers |
| `mva_genes_GRCh38.bed`, `mva_exons_GRCh38.bed` | bundled MVA gene spans and exons (MANE v1.5) |

### Main options

| Option | Use |
|---|---|
| `--vcf1`, `--vcf2`, `--label1/2`, `--combine union\|consensus` | the two phased, VEP-annotated VCFs |
| `-i hc_out` | `vcf_hc_vep.sh` output: failed-QC ClinVar calls and filtered variants |
| `--sv` (repeatable), `--sv-gt auto\|placeholder\|measured`, `--gene-gtf`, `--sv-pop` | SV input; gene models for full-line-only AnnotSV files or SV VCFs |
| `--bam`, `--bam-ref`, `--mosdepth`, `--exon-bed` | coverage checks and SV depth genotypes |
| `--exomiser-genes` (repeatable), `--alphagenome`, `--genes` | phenotype score, AlphaGenome table, extra gene panel |
| `--max-af`, `--cadd-min`, `--top-discovery`, `-t` | thresholds, number of discovery genes, threads |

### How the AnnotSV columns are used

| Columns | Use |
|---|---|
| `Annotation_mode` | Split lines give per-gene effects. SV-wide fields are copied from the full line. Full-line-only SVs → effect from gene models (`--gene-gtf`, else MVA genes only). |
| `SAMPLE` (GT) | A placeholder if all SVs share one GT-only value. With `--bam`, depth inside the SV vs 2–10 kb flanks: **DEL** <0.15× homozygous, 0.15–0.70× heterozygous, 0.70–0.85× partial (mosaic?), ≥0.85× not supported. **DUP** ≥1.8× homozygous/amplified, 1.3–1.8× heterozygous, 1.15–1.3× partial, <1.15× not supported. |
| `Location(2)`, `Overlapped_CDS_length/_percent`, `Frameshift` | effect class (frameshift or ≥50% CDS → SV_LOF; in-frame internal exons → SV_LOF_POSSIBLE) |
| `Dist_nearest_SS`, `Nearest_SS_type` | intronic SV ≤50 bp from a splice site → SV_SPLICE |
| `RE_gene` | regulatory element of a gene → SV_REGULATORY (SVs ≤1 Mb only) |
| `Tx` + MANE IDs in `INFO` | MANE transcript line preferred |
| `B_*_AFmax`, gnomAD AF columns and `INFO` keys | SVs with AF ≥ `--max-af` removed |
| `po_B_*_allG_source` / `P_*`, `po_P_*` | benign overlap of all genes (−1) / known pathogenic CNV (note) |
| `Gap_*`, `ENCODE_blacklist_*`, `SegDup_*`, `Repeat_type_*` | breakpoint confidence high / medium / low (low: −1) |
| `ACMG_class`, `AnnotSV_ranking_score/_criteria` | whole-SV class: raises a coding effect only (class 5 DEL → SV_LOF) |
| `HI`, `TS`, `GnomAD_pLI`, `LOEUF_bin`, OMIM / GenCC | `gene_info`, `gene_flag` (dosage-sensitive note) |
| `SUPP`/`SUPP_VEC`, `FILTER`, `IMPRECISE`, `Gene_count`, `CytoBand` | support and context notes |
| overlapping records of one event | merged (no false compound het) |
| coordinates vs the MVA genes | An SV overlapping an MVA gene (or ≤10 kb upstream) is a candidate even if AnnotSV doesn't name the gene. Old symbols are mapped (CCDC84→CENATAC, FAM178A→SLF2, …). Every SV within 50 kb of an MVA gene is listed in `mva_gene_svs_all.tsv` with its fate. |

### How the AlphaGenome table is used

| Signal | Class | Score |
|---|---|---:|
| `ag_splicing` and `pv_splice_max_absdiff` ≥0.2 | AG_SPLICE | 4 |
| `ag_splicing` | AG_SPLICE | 3 |
| `ag_expression` with log2FC <0 | AG_EXPRESSION | 3 |
| `ag_regulatory`, or evidence moderate/strong, or ≥2 cell types, not "up" | AG_REGULATORY | 2 |
| increase only, weak or no signal | note only | – |

- **Matching:** on chrom, pos, ref, alt and `marker_gene`.
- **Pre-filter:** the listed sites always pass the pre-filter, so they keep the VCF genotype and phase.
- **Quality:** a site with `quality_ok` False is excluded, unless the allele is ClinVar P/LP.
- **Not used for scoring:** `final_rank`, `track1_candidate`, `rank_key` and `final_tier`.

### Outputs (`-o`)

```
summary.md                          merge statistics, tiers, targets, next steps
targets.tsv / targets.bed           follow-up genes (tiers 1-4 + top discovery) and IGV positions
top_mva_variants.tsv                top variants (small + SV) in MVA genes (--top-mva, default 10)
ranked_genes.tsv                    tier gene gene_set genotype phase gene_score exomiser n_qualifying variants note
ranked_variants.tsv                 tier gene gene_set chrom pos id ref alt gt class score consequence impact hgvsc hgvsp
                                    max_af clinvar dp vaf evidence
01_table/annotated.tsv              merged table (SOURCE, GT_NOTE); *.prefilter.log, merge.log: counts
02_coverage/exon_depth.tsv          chrom start end exon mean_depth ratio_to_reference flag
02_coverage/karyotype_screen.tsv    chrom median_1Mb_window_depth ratio_to_autosomal_median windows_used flag
03_sv/mva_gene_svs_all.tsv          EVERY input SV within 50 kb of an MVA gene, with its fate (kept + class, or why removed)
03_sv/sv_candidates.tsv             gene chrom start end svtype svlen gt gt_source depth_ratio callers filter qual pr sr
                                    pop_af confidence class score effect exons acmg cytoband gene_info gene_flag notes ids
05_targets/target_gene_variants.tsv all small variants in target genes (both VCFs + filtered); target_gene_svs.tsv
logs/ranking.log, logs/warnings.txt
```
