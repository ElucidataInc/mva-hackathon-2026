## Team name	

RareELleles

## Model number (fill out this form once per model/approach submitted; up to 6)

1

# RareELleles: Track 1 methods report

MVA Hackathon 2026, Track 1 (Variant Prediction). Model 1.

## Please describe your model/approach in detail.

Approach:

1. Alignment and QC: reads aligned to GRCh38 (no-alt + hs38d1, GIAB v2 mask, matching the provided VCF) with bwa mem 0.7.19; duplicates marked with GATK MarkDuplicates; mean depth about 40x.

2. Small variants: DeepVariant 1.10.0 calls plus the provided Sentieon VCF (genotype concordance 0.994), annotated with VEP 116.2 (gnomAD v4.1, AlphaMissense) and ClinVar, and phased with WhatsHap 2.8.

3. Structural variants: Cue and Manta, cross-confirmed with truvari and annotated with gnomAD-SV and AnnotSV 3.5; exon read depth screened over the MVA genes.

4. Ranking: six independent prioritisation approaches were run on the proband's WGS and their outputs integrated into a single ranked list: 
(a, b) Exomiser 15.1.0 on the DeepVariant and the provided call sets (two lists), with the HPO terms from the provided phenotype document; 
(c) a tiered rule-based ranking combining small variants, SVs, phasing, Exomiser and AlphaGenome; 
(d) a ranking of rare variants by AlphaGenome Atlas (AVI) score; 
(e) AlphaGenome 0.9.0 scoring of rare variants within 50 kb of 11 MVA marker genes; 
(f) ACMG point-based scoring of marker-gene variants and nearby SVs.

5. Consensus: the top 30 of each list were matched by variant and given six evidence flags (MVA gene, two or more independent methods, protein-altering or SV, rare, ClinVar/AlphaMissense pathogenic, fits recessive inheritance). All six lists ranked BUB1B p.Leu737Ter and p.Asn1002Lys at the top.

EPCR scoring: each row was assigned a value in (0, 1] by evidence tier, in descending order: (i) a biallelic candidate in an established recessive MVA gene with a ClinVar pathogenic/likely pathogenic allele (0.95; below 1 because phase is unconfirmed); (ii) each allele of that pair alone (0.40 for the ClinVar P/LP stop-gain, 0.10 for the missense); (iii) rare non-coding variants near MVA marker genes with weak AlphaGenome support (0.05, 0.04); (iv) rare variants in genes not established for MVA, ordered by strength of supporting evidence (0.03 to 0.01): a predicted-damaging missense in ACLY, a 5' UTR splice-donor deletion in DIDO1, whose loss destabilises BubR1 and causes aneuploidy in mouse cells, and a frameshift in ARL9; (v) incidental findings (0.005, 0.004). All values are distinct, so each row forms its own F-max threshold.

## (required) If commercially-available Generative AI was used, record the provider, the plan or tier, and the relevant setting in your methods description. A line is enough, for example, "Anthropic API, Claude Sonnet, commercial terms, no training on customer content."

Anthropic Claude (Sonnet and Opus models), Elucidata business plan (Team/Enterprise), commercial terms, no training on customer content; no feedback or ratings were given on outputs. Used for writing and checking analysis scripts, summarising and cross-checking result tables, and drafting documentation.

## Please state whether the submission file is the automated output of your computational approach (preferred), or if it has undergone downstream manual review and curation.

The submission file has undergone downstream manual review and curation. The candidate lists are the automated outputs of six pipelines; the submitted rows, their order and EPCR values were set by manual review of those outputs, as described below.

## Please describe any downstream manual review in detail.

- Compared the top 30 of each of the six rankings, including the SVs, to focus on the strongest-scoring candidates of each method (148 distinct variants, 126 of them in only one list).
- Reviewed an evidence table with six flags per variant and warnings for ClinVar benign, AF of 1% or more, cis pairs and annotation errors.
- Checked top candidates in the BAM (both call sets, allele counts, strand, MAPQ). For example, a ClinVar P/LP call in PRSS1 was set aside: DeepVariant did not call it and all alternate reads at MAPQ 20 or more were on one strand.
- Checked WhatsHap phase for same-gene pairs: the ADAMTS1 and FANCD2 pairs are in cis; the BUB1B pair (about 11 kb apart) could not be phased.
- Grouped candidate genes by shared KEGG pathways and GO biological processes. Assessed the likely mechanism of each candidate: loss of function fits the recessive BUB1B MVA model (stop-gain with predicted nonsense-mediated decay, missense in the last exon); the non-coding candidates near BUB1B and SMC5 are listed as possible regulatory alleles for follow-up.
- Assigned EPCR values in descending order of supporting evidence: the BUB1B compound-heterozygous pair, then each BUB1B variant alone, then rare non-coding variants near MVA marker genes, then rare predicted-damaging variants in other genes; incidental findings received the lowest values.

## Please state whether your approach only used publicly available data (preferred) or if proprietary data was also used. If proprietary data were used, where possible, please also include a submission based only on publicly available data.

Only publicly available data were used. Apart from the dataset provided by the organisers (WGS FASTQs, VCF and phenotype document), all resources were public: the GRCh38 reference and GIAB masks, Ensembl and MANE gene models, gnomAD (v4.1 and gnomAD-SV), ClinVar, AlphaMissense, ReMM and the Exomiser data release, UCSC centromere models, the AlphaGenome model and AlphaGenome Atlas, and published MVA gene literature, and KEGG pathway and Gene Ontology annotations. No proprietary data were used for Track 1.

## Please describe any public data used in detail (e.g., reference panels, population databases, literature).

Reference: GRCh38 NCBI no-alt analysis set with hs38d1 decoys and the GIAB v2 masked exclusions. Gene models: Ensembl release 116 GTF and MANE Select v1.5. Population frequencies: gnomAD v4.1 (VEP 116 cache) and gnomAD-SV v4.1 (via SVAFotate). Clinical and predictive annotation: ClinVar (23 September 2026 release, plus ClinVar 202509 in the VEP cache), AlphaMissense (hg38), ReMM v0.4, Exomiser data release 2602, AnnotSV 3.5 annotation sources, the AlphaGenome 0.9.0 model and the AlphaGenome Atlas of precomputed variant scores (queried by genomic position only). Other: UCSC hg38 centromere models (downloaded 6 October 2026) and published MVA literature, used to define the 11 marker genes and the MVA1 allele rs576524605. Pathway and functional annotation: Gene Ontology terms from the Polly biomedical knowledge graph, KEGG pathway and network annotations, and g:Profiler enrichment of the marker genes.

## Please describe any proprietary data used in detail.

None

## Is your approach able to output proposed pairs of compound heterozygous candidate variants, or only single candidate variants?

Yes. Compound-heterozygous pairs are proposed when a gene carries two rare, predicted-damaging heterozygous variants (Exomiser autosomal recessive model and the two-hit tiers of the rule-based ranking). WhatsHap phase sets are used to exclude pairs in cis. Pairs whose phase cannot be resolved from short reads, as for BUB1B, are proposed as one row with the trans configuration unconfirmed.

## How did you handle secondary or incidental findings, if any? (See the finding_type column in your submission.)

Two incidental findings are reported with finding_type secondary: LZTR1 p.(Tyr748Ter) (heterozygous stop-gain, ClinVar pathogenic/likely pathogenic) and GNRHR p.(Arg139His) (heterozygous, ClinVar pathogenic; carrier finding in an autosomal recessive gene). Neither gene is established for MVA; they were identified while reviewing ClinVar P/LP calls and are listed for clinical review only. In the KEGG and GO annotation of the candidate genes, GNRHR maps to GnRH receptor signalling (GnRH-GnRHR-PLCB-PKC) and LZTR1 to RAS ubiquitination by the CUL3 ubiquitin ligase complex.

## Please provide an estimate of run time and cost for your approach (if possible).

All analyses ran on a single cloud server: Intel Core i9-14900K (1 socket, 32 logical CPUs), 125 GiB RAM, 1.8 TB disk, Ubuntu 24.04.5 LTS, no GPU. Approximate compute time per stage: reference build, alignment, QC, DeepVariant calling, phasing, Exomiser and annotation about 14 hours (alignment 4.4 hours and DeepVariant 3.5 hours are the longest steps); structural-variant calling with Cue and Manta about 5 hours; structural-variant annotation about 1 hour; Exomiser on the provided calls about 40 minutes; AlphaGenome marker-gene scoring 1.5 to 3 hours; AlphaGenome Atlas ranking about 2 hours; functional annotation about 1 hour; final ranking 1 to 2 hours. Run one after another, this is about 26 to 29 hours of compute; several stages can run in parallel. Compute cost was about $236.

## Provide a method abstract (up to 500 words). Please include strengths and limitations of your method as well as methodological detail.

Candidate causal variants for the proband were prioritised from whole-genome sequencing by combining independent variant calling, phasing, structural-variant detection and ranking approaches based on variant type, MVA association, phase, annotation and mechanism, and the submitted rows were selected by manual review of cross-method agreement and supporting evidence.

Methods: Reads were aligned to GRCh38 (NCBI no-alt analysis set with hs38d1 decoys and the GIAB v2 mask, matching the provided VCF) with bwa mem and duplicate-marked with GATK; sample identity was confirmed against the provided VCF (concordance 0.994). Small variants were called with DeepVariant 1.10.0 and analysed alongside the provided Sentieon calls; both sets were normalised, filtered to PASS, annotated with Ensembl VEP 116.2 (gnomAD v4.1, AlphaMissense) and ClinVar, and read-backed phased with WhatsHap. Structural variants were called with Cue and Manta, cross-confirmed with truvari and annotated with gnomAD-SV and AnnotSV, and exon read depth was screened over the MVA genes. Six lists were produced: Exomiser on each call set with the provided HPO terms; a tiered rule-based ranking that integrates small variants, SVs, phase, Exomiser and AlphaGenome and puts biallelic genotypes in MVA genes first; a rare-variant ranking by AlphaGenome Atlas (AVI) score; AlphaGenome scoring of rare variants around 11 MVA genes; and ACMG point-based scoring of marker-gene variants and SVs. The top 30 of each list were matched by variant and annotated with six evidence flags.

Results: All six lists placed BUB1B c.2210T>G p.(Leu737Ter) and c.3006T>G p.(Asn1002Lys) at the top. Below these, agreement was low: 148 distinct variants, 126 of them in a single list. The submitted rows, their order and EPCR values were set by an evidence-tier scoring (variant type, MVA association, phase and annotation) and manual curation based on annotation, associated pathways and mechanism.

Strengths: two independent call sets; methods with different assumptions (phenotype-driven, rule-based, regulatory prediction, ACMG-style scoring); SV and exon CNV screening; read-level and phase review; functional annotation and mechanism view; and only public resources with versioned software environments.

Limitations: a single proband without parental samples, so the two BUB1B variants (about 11 kb apart) could not be phased and their trans configuration is unconfirmed; CADD was not available, and REVEL and SpliceAI scores were available only through Exomiser; depth-based exon CNVs and single-caller SVs are unconfirmed; and distal regulatory predictions are less reliable.