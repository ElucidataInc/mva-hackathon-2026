# 07. Data, tools, licences and compute

The hackathon dataset, the reference resources a Track 1 pipeline needs, tool licences and versions, recommended compute and storage, and a data-deletion map. The dataset manifest uses public metadata only.

## Summary

- The gated dataset is 85.0 GB: 8 FASTQ.gz files (4 lanes, paired), a VCF with index (317.5 MB), an HPO phenotype docx and a README.
- The working reference set is about 350 GB. The VEP 116 cache already carries gnomAD v4.1 frequencies and ClinVar, so the 761 GB of gnomAD VCFs are optional.
- A full run from FASTQ needs about 1-1.3 TB of storage in total, including alignments and pipeline scratch.
- Recommended machine: 48 vCPU and 128 GB RAM, plus a separate 2 TB persistent volume. Provisioned machine: PNAP s3.c3.medium: Intel i9-14900K, 32 threads, 128 GB RAM (125 GiB usable), 2 × 2 TB NVMe. A full run from FASTQ needs about 1 TB of storage and takes about 7-9 hours; Track 2 needs little compute and can run on laptops.
- Non-commercial only: MELT, xTea, Manta, LIRICAL, CADD (scores and, read literally, scripts), SpliceAI scores and models, and the dbNSFP academic branch.
- Permissive: GATK, DeepVariant, WhatsHap, bwa-mem2, CNVkit, gatk-sv, slivar, Talos, nf-cavalier (Apache-2.0, BSD-3 or MIT). TIDDIT is GPL-3.0. Exomiser is AGPL-3.0, which matters only if a modified version is offered as a service.
- AlphaMissense predictions are CC BY 4.0, so it can replace CADD and SpliceAI for missense scoring.

## 1. Dataset

From the Hugging Face API (`/api/datasets/SageBio/mva-hackathon-2026-data` and `/tree/main`). The file listing is public; the file contents are gated.

| File | Bytes | Size | Contents |
|---|---|---|---|
| `.gitattributes` | 2,575 | 2.5 KiB | Git LFS config |
| `README.md` | 5,076 | 5.0 KiB | Dataset card (gated) |
| `Challenge_Clinical_Phenotype_1.docx` | 16,865 | 16.5 KiB | HPO terms |
| `WGS_EX2312012_HGWCNDSX7.vcf.gz` | 315,153,971 | 300.6 MiB | Variant calls: GRCh38 no-alt + hs38d1 with GRC exclusion masking v2, Sentieon 202308.02 (VCF header) |
| `WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi` | 2,343,376 | 2.2 MiB | Tabix index |
| `..._S16_L001_R1_001.fastq.gz` | 10,467,073,561 | 9.75 GiB | Lane 1, read 1 |
| `..._S16_L001_R2_001.fastq.gz` | 11,035,873,766 | 10.28 GiB | Lane 1, read 2 |
| `..._S16_L002_R1_001.fastq.gz` | 10,375,801,767 | 9.66 GiB | Lane 2, read 1 |
| `..._S16_L002_R2_001.fastq.gz` | 10,911,120,034 | 10.16 GiB | Lane 2, read 2 |
| `..._S16_L003_R1_001.fastq.gz` | 10,195,738,985 | 9.50 GiB | Lane 3, read 1 |
| `..._S16_L003_R2_001.fastq.gz` | 10,621,354,344 | 9.89 GiB | Lane 3, read 2 |
| `..._S16_L004_R1_001.fastq.gz` | 10,306,364,082 | 9.60 GiB | Lane 4, read 1 |
| `..._S16_L004_R2_001.fastq.gz` | 10,755,107,565 | 10.02 GiB | Lane 4, read 2 |
| **Total** | **84,985,953,392** | **85.0 GB / 79.1 GiB** | |

FASTQ names share the prefix `WGS_EX2312012_HGWCNDSX7`.

## 2. Reference and annotation resources

| Resource | Version | Date | Size | Source |
|---|---|---|---|---|
| GRCh38 no-alt + hs38d1 FASTA | GCA_000001405.15 | 2013 assembly | 874.6 MB compressed | `ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz` |
| VEP indexed cache, GRCh38 (includes gnomAD v4.1 exome and genome frequencies, ClinVar 2025-09) | Ensembl 116 | 2026 | 27.6 GB | `ftp.ensembl.org/pub/current_variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz` |
| MANE | v1.5 | | small | `ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/` |
| gnomAD v4.1 exomes sites (optional) | v4.1 | 2024 | 198.4 GB (24 files) | `gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/` |
| gnomAD v4.1 genomes sites (optional) | v4.1 | 2024 | 563.0 GB (24 files) | `gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/` |
| ClinVar weekly VCF | weekly | | ~200 MB | `ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/` |
| CADD SNVs, GRCh38 | v1.7 | | 87.5 GB | `kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz` |
| CADD gnomAD r4.0 indels | v1.7 | | 1.26 GB | `kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz` |
| SpliceAI precomputed scores | Jaganathan 2019 | | not measured (BaseSpace login) | `basespace.illumina.com/s/otSPW8hnhaZR` |
| AlphaMissense hg38 | 2023 release | 2023-08-03 | 643.0 MB | `storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz` |
| REVEL | v1.3 | | 667.1 MB | `zenodo.org/records/7072866` (mirror) |
| HPO | 2026-09-01 | 2026-09-03 | 35.8 MB (`phenotype.hpoa`) | `purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa` |
| Exomiser CLI | 15.1.0 | 2026-06-09 | 170.4 MB | `github.com/exomiser/Exomiser/releases` |
| Exomiser data | 2512 (hg38 + phenotype) | 2025-12 | 23.6 GB + 13.1 GB zipped; about 80 GB unpacked | `data.monarchinitiative.org/exomiser/data/` |
| bwa-mem2 index, GRCh38 | built locally | | ~10 GB (building it needs ~85-90 GB RAM) | |
| REMM (Genomiser) | v0.4 hg38 | 2022-05-24 | 16.0 GB | `kircherlab.bihealth.org/download/ReMM/ReMM.v0.4.hg38.tsv.gz` |

Without the gnomAD VCFs, the working set (FASTA and index, VEP cache, Exomiser data, CADD, REMM, AlphaMissense, REVEL, HPO, containers) is about 350 GB. Adding the gnomAD VCFs brings it to about 1.1 TB.

## 3. Licences

Hackathon outputs are CC BY 4.0, and this pipeline may be reused commercially, so each licence is assessed for commercial use. This is a technical summary, not legal advice.

| Resource | Licence | Commercial use | Notes |
|---|---|---|---|
| CADD scores | "freely available for all non-commercial applications" (cadd.gs.washington.edu/faq) | No, without a UW CoMotion licence | Use AlphaMissense instead |
| CADD scripts | LICENSE grants rights "to all non-commercial users and licensees of CADD" | No, read literally | |
| SpliceAI code | PolyForm Strict 1.0.0 | No | |
| SpliceAI scores and models | Scores free for academic and non-profit use; models CC BY-NC 4.0 | No, without an Illumina licence | |
| AlphaMissense code | Apache-2.0 | Yes | Model weights not released |
| AlphaMissense predictions | CC BY 4.0 | Yes, with attribution | |
| REVEL | Not published on the download page | Unknown | Also bundled in the non-commercial dbNSFP academic branch |
| dbNSFP academic | CC BY-NC-ND 4.0 | No | |
| dbNSFP commercial | Paid (Genos Bioinformatics): $5,000/yr internal research, $10,000/yr product development | Yes, paid | Excludes CADD, VEST, M-CAP, MutScore, PolyPhen-2, PrimateAI and RGC Million Exome |
| Exomiser code | AGPL-3.0 | Yes | A modified version offered as a network service must publish its source |
| Exomiser data | No licence statement found | Unknown | |
| Talos, nf-cavalier, WhatsHap | MIT | Yes | |
| Manta | PolyForm Strict 1.0.0 | No | |
| TIDDIT | GPL-3.0 | Yes | Copyleft on distributed modifications |
| RetroSeq (MEI caller in nf-core/raredisease) | No licence file | Unknown | |
| DeepVariant | BSD-3-Clause | Yes | |
| GATK | Apache-2.0 (top-level LICENSE.TXT) | Yes | GitHub shows "Other" because files are licensed individually |
| VEP | Apache-2.0 | Yes | Plugins carry their own data licences |
| MELT | Academic non-commercial; bars "commercial, public, cloud, time-sharing or service bureau use" | No, without a University of Maryland licence | |
| xTea | Harvard academic licence; bars commercial and fee-for-service use | No | |
| SCRAMble | Unknown | Unknown | |
| gnomAD data | See gnomad.broadinstitute.org/terms | Unknown | Browser code is MIT |
| ClinVar | NCBI | Usually public domain | |
| OMIM API | Registration required | Usually needs a licence | Research API key; commercial use has needed a licence |
| DrugBank | DrugBank terms | Usually needs a licence | Usually needs a commercial licence beyond academic use |
| Open Targets | See platform-docs.opentargets.org/licence | Unknown | |
| ChEMBL | CC BY-SA 3.0 | Yes | Derived datasets must carry CC BY-SA 3.0 |

## 4. Tool versions

From the GitHub Releases API unless stated.

| Tool | Version | Resources |
|---|---|---|
| nf-core/raredisease | 3.1.2 (2026-07-06) | Starts from FASTQ. bwa-mem2 alignment, DeepVariant (default) for SNVs and indels, Manta and TIDDIT for SVs, RetroSeq for mobile elements, then VEP and ranking. Nextflow >= 25.10.4. Default resources suit the test data only and must be raised for a genome |
| nf-core/sarek | 3.10.0 (2026-08-12) | General germline/somatic pipeline |
| bwa-mem2 | 2.3 (2025-06-30) | Index ~10 GB on disk and in memory; building it needs about 28 GB per gigabase of reference (~85-90 GB for GRCh38) |
| DeepVariant | 1.10.0 (2026-03-05) | Whole genome in 69 min on 96 vCPU / 384 GiB: `make_examples` 46 min, `call_variants` 16 min, `postprocess_variants` 7 min. A GPU speeds up only `call_variants` |
| GATK | 4.7.0.0 (2026-08-18) | |
| Exomiser | 15.1.0 (2026-06-09) | 8 GB RAM minimum for a genome, 16 GB recommended; at least 100 GB free disk, SSD preferred |
| LIRICAL | 2.4.1 (2026-06-18) | |
| slivar | 0.3.4 (2026-05-02) | Single Go binary |
| WhatsHap | 2.8 (PyPI) | |
| CNVkit | 0.9.14 (2026-07-01) | |
| Manta | 1.6.0 (2019-07-09) | No release since 2019 |
| Talos | No tagged release | Cohort reanalysis against ClinVar and PanelApp Australia |
| nf-cavalier | No release; run from `main` | Nextflow >= 24.04. Annotates with vcfanno, SVAFotate and VEP |

nf-core/raredisease fits Track 1 better than sarek, which targets cancer and large cohorts. The Exomiser CLI and its data release must match; 2512 is the newest data release on the public mirror.

## 5. Compute and storage

### 5.1 Recommended machine

| Resource | Recommended | Why |
|---|---|---|
| CPU | 48 vCPU | With 128 GB RAM, alignment and DeepVariant can safely use about 48 threads. More cores would sit idle during those steps |
| RAM | 128 GB | bwa-mem2 has peaked at 30-110 GB with 36 threads on hg38, and each DeepVariant `make_examples` shard uses 1-2 GB, so both are capped at about 40-48 threads. Building the bwa-mem2 index needs ~85-90 GB once. 256 GB or more removes these caps |
| Storage | 2 TB persistent volume | Dataset, references, Nextflow `work/`, temp files, final CRAM and results. Keeping everything on the volume keeps patient data in one place |
| GPU | Not needed | Only DeepVariant `call_variants` benefits; `make_examples` is the larger stage and runs on CPU |

Settings for 128 GB: cap Nextflow at 120 GB and 48 CPUs; run bwa-mem2 with 32-40 threads and `samtools sort -@ 8 -m 2G`; run DeepVariant with `--num_shards 48`; build the bwa-mem2 index with nothing else running; never run alignment and DeepVariant at the same time. nf-core/raredisease defaults these steps to 12 CPUs and 72 GB, so raise them in the config.

Use a dedicated volume for the hackathon data rather than a shared one, so it can be wiped and deleted in one step by 23 Nov 2026. On the provisioned bare-metal server the data sit on the local NVMe, so releasing the server deletes them; keep it running until the work is finished and the deletion is confirmed.

### 5.2 Storage budget

| Item | Size |
|---|---|
| Hackathon dataset | 85 GB |
| Aligned reads: BAM during the run, CRAM kept | ~100 GB + ~45 GB |
| Pipeline scratch (sorting, `work/`) | ~200-300 GB at peak |
| GRCh38 FASTA and bwa-mem2 index | ~20 GB |
| VEP 116 cache | ~30 GB |
| Exomiser data 2512 | ~80 GB |
| CADD v1.7 SNVs and indels | ~89 GB |
| REMM | 16 GB |
| AlphaMissense, REVEL, HPO, ClinVar, chr15 phasing panel | ~5 GB |
| Containers and tools | ~30 GB |
| Track 2 databases (optional local ChEMBL) | 0-30 GB |
| **Total** | **~0.9-1.1 TB** |
| Optional gnomAD v4.1 VCFs | +761 GB |

### 5.3 Runtime

| Step | 48 vCPU, 128 GB |
|---|---|
| bwa-mem2 alignment and sorting | ~3-4 h |
| DeepVariant | ~2-2.5 h |
| SV, MEI and CNV calling, phasing | ~1 h |
| Annotation and Exomiser ranking | ~0.5-1 h |
| **Full Track 1 run from FASTQ** | **~7-9 h** |
| Re-ranking from existing calls | ~1-3 h |

These times are estimates extrapolated from the DeepVariant benchmark and bwa-mem2 scaling; they have not been measured on this dataset. Track 2 uses gene- and variant-level results only and can run on laptops.

### 5.4 Data-deletion map

Everything derived from the dataset must be deleted by 23 Nov 2026 ([document 06](06_hackathon_rules.md), section 4).

| Artefact | Action |
|---|---|
| FASTQ and VCF downloads | Delete |
| BAM/CRAM and any re-sorted or re-indexed copy | Delete |
| VCF splits, annotated VCFs, genotype-scale TSV/Parquet exports | Delete |
| Nextflow `work/`, `.nextflow.log`, cached intermediate data | Delete |
| Notebook checkpoints and kernel state with patient data | Delete |
| LLM prompts and chat logs containing variant rows or reads | Delete (provider-side logs are exempt if the provider meets the processor test in thread #2) |
| Any model trained on the genome or its intermediates | Delete |
| Public reference downloads | Keep (not patient data) |
| Submission CSV, reports, HPO term list, gene and pathway rankings, mechanism write-up, drug list, code, pitch | Keep |

## 6. Candidate tool stack

| Step | Tool | Version | Licence |
|---|---|---|---|
| Alignment | bwa-mem2 | 2.3 | MIT |
| Orchestration | nf-core/raredisease | 3.1.2 | MIT; containers carry their own licences |
| SNV/indel calling | GATK or DeepVariant | 4.7.0.0 / 1.10.0 | Apache-2.0 / BSD-3 |
| Structural variants | Manta, TIDDIT, CNVkit | 1.6.0 / 0.9.14 | PolyForm Strict (non-commercial) / GPL-3.0 / Apache-2.0 |
| Mobile-element insertions | RetroSeq (in nf-core/raredisease), MELT or xTea | | RetroSeq has no licence file; MELT and xTea are non-commercial |
| Phasing | WhatsHap | 2.8 | MIT |
| Filtering | slivar | 0.3.4 | MIT |
| Annotation | VEP | 116 | Apache-2.0 |
| Pathogenicity | AlphaMissense; REVEL once its licence is known | 2023-08-03 / 1.3 | CC BY 4.0 / unknown |
| Population frequency | gnomAD via the VEP cache | 4.1 | See gnomAD terms |
| Clinical database | ClinVar (also in the VEP cache) | weekly | Usually public domain |
| Phenotype ranking | Exomiser; LIRICAL as secondary | 15.1.0 / 2.4.1 | AGPL-3.0 / non-commercial |
| Reanalysis | Talos | | MIT |
| Reporting | nf-cavalier | `main` | MIT |
| Compute | PNAP s3.c3.medium: Intel i9-14900K, 32 threads, 128 GB RAM (125 GiB usable), 2 × 2 TB NVMe | | |

## Open questions

1. REVEL licence terms.
2. gnomAD, ClinVar and Open Targets data terms.
3. SCRAMble repository and licence, and RetroSeq's licence.
4. Measured runtimes for alignment and the full pipeline on this dataset.
5. Library insert size (thread #25), which affects read-backed phasing.

## References

1. Hugging Face dataset API, `SageBio/mva-hackathon-2026-data` card and file tree.
2. NCBI GRCh38 analysis-set FASTA (URL in section 2).
3. Ensembl REST `/info/data`; VEP 116 cache and licence (`github.com/Ensembl/ensembl-vep`, release/116).
4. NCBI MANE directory listing.
5. gnomAD v4.1 sites VCFs via the Google Cloud Storage JSON API; gnomAD terms page; gnomAD browser README.
6. NCBI ClinVar weekly VCF directory.
7. CADD downloads, FAQ and CADD-scripts LICENSE (`github.com/kircherlab/CADD-scripts`).
8. Illumina SpliceAI README (`github.com/Illumina/SpliceAI`).
9. Google DeepMind AlphaMissense README and LICENSE (`github.com/google-deepmind/alphamissense`); hg38 predictions file.
10. REVEL download page and Zenodo mirror (record 7072866).
11. dbNSFP licence page (`dbnsfp.org/license`).
12. HPO releases (`github.com/obophenotype/human-phenotype-ontology`).
13. Exomiser releases, CLI README and data mirror (`data.monarchinitiative.org/exomiser/data/`).
14. REMM downloads (Kircher lab).
15. GitHub licence API for gatk, deepvariant, whatshap, manta, gatk-sv, pbsv, talos, nf-cavalier, Exomiser, LIRICAL, slivar, cnvkit, bwa-mem2 and xTea.
16. MELT downloads and licence (`melt.igs.umaryland.edu/downloads.php`).
17. xTea LICENSE (`github.com/parklab/xTea`).
18. GATK LICENSE.TXT (`github.com/broadinstitute/gatk`).
19. WhatsHap PyPI metadata.
20. ChEMBL "About" page (data licence).
21. Open Targets licence page; DrugBank terms; OMIM API help.
22. nf-core/raredisease README and usage docs; raredisease and sarek release tags.
23. DeepVariant `docs/metrics.md` and `docs/deepvariant-details.md` (r1.10).
24. nf-cavalier and Talos READMEs.
25. bwa-mem2 README (`github.com/bwa-mem2/bwa-mem2`).
26. Exomiser installation docs (`exomiser.readthedocs.io`).
27. Ensembl VEP cache documentation (release 116 cache contents).
28. Manta LICENSE.txt (`github.com/Illumina/manta`); TIDDIT COPYING (`github.com/SciLifeLab/TIDDIT`); RetroSeq repository (`github.com/tk2/RetroSeq`).
29. [Document 06](06_hackathon_rules.md) in this folder.
