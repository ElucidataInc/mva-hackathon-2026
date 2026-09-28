# 03. *BUB1B* variant spectrum and pipeline failure modes

*BUB1B* gene facts, ClinVar counts for five MVA genes, the known pathogenic *BUB1B* variants with gnomAD frequencies, two published founder and noncoding alleles, and the ways a singleton short-read WGS pipeline can miss a causal compound heterozygote.

## Summary

- MANE Select transcript: ENST00000287598.11 / NM_001211.6 (NP_001202.5), 23 exons, 1050 aa, chr15:40,161,069-40,221,123 (GRCh38, + strand).
- *BUB1B* tolerates heterozygous loss of function (pLI ~0, LOEUF 0.75), as expected for a recessive gene.
- gnomAD lists 64 pathogenic or likely pathogenic *BUB1B* variants on the MANE transcript. 30 are absent from gnomAD v4 and the rest have allele counts of 1-6 in 0.7-1.6 million alleles.
- Matsuura's recurrent 1833delT is c.1833del at chr15:40,206,279. All 3 gnomAD carriers are East Asian (popmax AF 7.6e-5), so standard frequency filters keep it.
- The Ochiai 2014 upstream allele (ss802470619) has no resolved GRCh38 coordinate.
- Exomiser's exome preset drops upstream, intergenic and regulatory variants. Only the genome preset (Genomiser) keeps them.
- Read-backed phasing cannot link variants kilobases apart with 300-500 bp fragments. Phase for a *BUB1B* pair will usually rest on population data and published precedent.

## 1. *BUB1B* gene facts

| Field | Value | Source |
|---|---|---|
| Symbol / Entrez ID | *BUB1B* / 701 | NCBI Gene |
| Ensembl gene | ENSG00000156970, chr15:40,160,984-40,221,137 (+), 15q15.1 | Ensembl REST |
| MANE Select | ENST00000287598.11 / NM_001211.6, NP_001202.5, CCDS10053.1 | Ensembl REST, NCBI, gnomAD API |
| MANE span | chr15:40,161,069-40,221,123 | Ensembl REST |
| Exons / protein | 23 / 1050 aa | NCBI gene table |
| Other transcripts | At least 8 more protein-coding, 31 in total | Ensembl REST |

### 1a. gnomAD v4 constraint

| Metric | Value |
|---|---|
| pLI | 5.52e-17 |
| LOEUF | 0.747 |
| oe_lof | 0.624 (90% CI 0.525-0.747) |
| oe_mis | 0.914 (90% CI 0.870-0.960) |
| LoF Z / missense Z / synonymous Z | 3.74 / 1.30 / -0.06 |

From the gnomAD GraphQL API; the gene returned 4,974 v4 variants. Strongly constrained genes have LOEUF below about 0.35, so heterozygous *BUB1B* null carriers are expected to be unaffected and not rare.

### 1b. Paralogs and repeats

- *BUB1* (chr2) has a near-identical constraint profile. The two genes share a kinase domain but sit on different chromosomes, so segmental-duplication mismapping between them is not expected. Local misalignment within the kinase-domain exons was not tested.
- The UCSC `genomicSuperDups` track (hg38) has no entry overlapping *BUB1B* plus 1 kb flanks.

## 2. ClinVar counts for five MVA genes

NCBI `esearch` term counts. A variant can fall under several classification terms, so the columns overlap and sum to more than the total.

| Gene | Total | Pathogenic | Likely pathogenic | VUS | Likely benign | Benign | Conflicting |
|---|---|---|---|---|---|---|---|
| *BUB1B* | 2,501 | 141 | 32 | 1,682 | 958 | 1,068 | 140 |
| *CEP57* | 1,005 | 50 | 8 | 690 | 346 | 369 | 48 |
| *TRIP13* | 418 | 153 | 9 | 301 | 289 | 292 | 36 |
| *CENATAC* | 126 | 34 | 4 | 93 | 89 | 89 | 68 |
| *MAD1L1* | 261 | 53 | 6 | 190 | 179 | 181 | 122 |

*BUB1B* by consequence: missense 1,497, synonymous 614, intron 459, 5' UTR 121, nonsense 41, frameshift 41, splice acceptor 8, splice donor 7. These counts include non-MANE transcripts and benign variants.

## 3. Pathogenic and likely pathogenic *BUB1B* variants

From gnomAD's `clinvar_variants` resolver on the MANE transcript, joined with gnomAD v4 frequencies. Use this list rather than the raw NCBI count of 141, which double-counts. All 64 are frameshift or stop-gained except two splice variants (c.1402-5A>G, c.1402-1G>T) and c.1833del. Stars are ClinVar review status. AF is gnomAD v4 exome / genome; a dash means not observed.

| GRCh38 chr15 | HGVS.c | HGVS.p | ClinVar | Stars | AF exome / genome |
|---|---|---|---|---|---|
| 40165186 | c.169C>T | p.Gln57Ter | Pathogenic | 1 | - / - |
| 40170073 | c.192del | p.Tyr64Ter | Pathogenic | 1 | - / - |
| 40170081 | c.199C>T | p.Arg67Ter | Pathogenic | 2 | 1.30e-05 / 6.57e-06 |
| 40170539 | c.243del | p.Tyr81Ter | Pathogenic/Likely pathogenic | 2 | 6.84e-07 / - |
| 40170578 | c.282_283del | p.Ser96Ter | Pathogenic | 1 | 6.84e-07 / - |
| 40170604 | c.310del | p.Arg104GlufsTer3 | Pathogenic | 1 | - / - |
| 40170637 | c.340C>T | p.Arg114Ter | Pathogenic | 2 | 6.84e-06 / - |
| 40170655 | c.358C>T | p.Arg120Ter | Pathogenic | 2 | 2.05e-06 / - |
| 40176513 | c.422dup | p.Tyr141Ter | Pathogenic | 1 | - / - |
| 40176555 | c.464del | p.Tyr155LeufsTer32 | Pathogenic | 1 | - / - |
| 40176655 | c.564_571dup | p.Ser191TyrfsTer40 | Pathogenic | 1 | 1.37e-06 / - |
| 40176669 | c.578del | p.His193ProfsTer35 | Pathogenic | 1 | 6.84e-07 / - |
| 40176672 | c.580C>T | p.Arg194Ter | Pathogenic/Likely pathogenic | 2 | 6.16e-06 / 6.58e-06 |
| 40183727 | c.595C>T | p.Arg199Ter | Pathogenic/Likely pathogenic | 2 | 2.74e-06 / - |
| 40183825 | c.693_694insTT | p.Ser232LeufsTer15 | Pathogenic | 1 | - / - |
| 40183838 | c.709_712del | p.Thr237GlnfsTer8 | Pathogenic | 1 | - / 6.57e-06 |
| 40185212 | c.799C>T | p.Gln267Ter | Pathogenic | 2 | 6.84e-07 / 6.58e-06 |
| 40185212 | c.803dup | p.Asn268LysfsTer2 | Pathogenic | 1 | - / - |
| 40185283 | c.871del | p.Thr291GlnfsTer6 | Pathogenic | 1 | 5.47e-06 / - |
| 40185290 | c.877C>T | p.Gln293Ter | Pathogenic | 1 | 6.84e-07 / - |
| 40185304 | c.897dup | p.Met300HisfsTer8 | Pathogenic | 1 | 6.84e-07 / - |
| 40185304 | c.897del | p.Met300CysfsTer31 | Pathogenic | 1 | 1.37e-06 / 6.58e-06 |
| 40185326 | c.913G>T | p.Glu305Ter | Pathogenic | 1 | - / - |
| 40185338 | c.925C>T | p.Gln309Ter | Pathogenic | 1 | - / - |
| 40185628 | c.1045del | p.Arg349AspfsTer5 | Pathogenic | 1 | 1.37e-06 / - |
| 40196624 | c.1141del | p.Arg381GlyfsTer16 | Pathogenic | 1 | - / - |
| 40196633 | c.1147C>T | p.Gln383Ter | Pathogenic | 1 | - / - |
| 40196684 | c.1198A>T | p.Lys400Ter | Pathogenic | 1 | 4.79e-06 / - |
| 40199624 | c.1298T>A | p.Leu433Ter | Pathogenic | 1 | - / - |
| 40199626 | c.1302del | p.Ser435ValfsTer19 | Pathogenic | 1 | 6.84e-07 / - |
| 40199633 | c.1313_1317del | p.Lys438SerfsTer8 | Pathogenic | 1 | - / - |
| 40199704 | c.1378C>T | p.Gln460Ter | Pathogenic | 1 | 1.37e-06 / - |
| 40199713 | c.1387A>T | p.Arg463Ter | Pathogenic | 1 | - / - |
| 40200239 | c.1402-5A>G | (splice region) | Pathogenic | 0 | 2.06e-06 / - |
| 40200243 | c.1402-1G>T | (splice acceptor) | Pathogenic | 0 | - / - |
| 40200283 | c.1441C>T | p.Gln481Ter | Pathogenic | 1 | - / - |
| 40200306 | c.1466_1470dup | p.Gly491TyrfsTer4 | Pathogenic | 1 | - / - |
| 40200354 | c.1512T>A | p.Cys504Ter | Pathogenic | 1 | - / - |
| 40202608 | c.1648C>T | p.Arg550Ter | Pathogenic | 1 | 2.05e-06 / 1.97e-05 |
| 40202623 | c.1663C>T | p.Arg555Ter | Pathogenic | 1 | 4.10e-06 / - |
| 40202628 | c.1672del | p.Ala559GlnfsTer29 | Pathogenic | 1 | - / - |
| 40206189 | c.1743del | p.Phe581LeufsTer7 | Pathogenic | 1 | 2.05e-06 / - |
| 40206279 | c.1833del | p.Phe611LeufsTer15 | Pathogenic | 0 | 2.05e-06 / - |
| 40206354 | c.1906del | p.Glu636ArgfsTer9 | Pathogenic | 1 | - / - |
| 40206403 | c.1954C>T | p.Gln652Ter | Pathogenic | 1 | 6.84e-07 / - |
| 40209643 | c.2153del | p.Thr718IlefsTer14 | Pathogenic | 1 | 6.84e-06 / - |
| 40209697 | c.2208_2211dup | p.Ser738ValfsTer16 | Pathogenic | 1 | 1.71e-05 / - |
| 40209699 | c.2210del | p.Leu737Ter | Pathogenic | 1 | - / - |
| 40209701 | c.2210T>G | p.Leu737Ter | Pathogenic/Likely pathogenic | 2 | 7.87e-05 / 3.29e-05 |
| 40210123 | c.2298C>G | p.Tyr766Ter | Pathogenic | 1 | - / - |
| 40210133 | c.2308C>T | p.Arg770Ter | Pathogenic | 1 | 1.44e-05 / 6.58e-06 |
| 40210141 | c.2316C>G | p.Tyr772Ter | Pathogenic | 1 | - / - |
| 40210151 | c.2328dup | p.Asp777ArgfsTer41 | Pathogenic | 2 | - / 1.31e-05 |
| 40210159 | c.2334C>G | p.Tyr778Ter | Pathogenic | 1 | - / - |
| 40210185 | c.2363_2364del | p.Ser788CysfsTer29 | Pathogenic/Likely pathogenic | 2 | 4.81e-06 / 1.32e-05 |
| 40210188 | c.2364dup | p.Ala789CysfsTer29 | Pathogenic | 1 | - / - |
| 40212546 | c.2436_2439dup | p.Arg814GlyfsTer5 | Pathogenic | 1 | - / - |
| 40212591 | c.2478T>A | p.Tyr826Ter | Pathogenic | 1 | - / - |
| 40212591 | c.2478T>G | p.Tyr826Ter | Pathogenic | 1 | - / - |
| 40213359 | c.2566del | p.His856MetfsTer3 | Pathogenic | 2 | 2.74e-06 / - |
| 40213431 | c.2636dup | p.His880ProfsTer3 | Pathogenic | 1 | 6.84e-07 / - |
| 40217665 | c.2848C>T | p.Gln950Ter | Pathogenic | 1 | 2.05e-06 / - |
| 40218481 | c.2876T>G | p.Leu959Ter | Pathogenic | 1 | - / - |
| 40218501 | c.2896A>T | p.Lys966Ter | Pathogenic | 1 | - / - |

c.2386-11A>G (Rio Frio 2010) is at chr15:40,212,488 but is classified "Conflicting classifications of pathogenicity" (1 star), so it is not in the table.

## 4. Two published *BUB1B* alleles

### 4a. Matsuura 2006 "1833delT"

- Recurrent in 4 of 7 Japanese PCS families.
- NM_001211.6:c.1833del, p.(Phe611LeufsTer15). ClinVar VCV000006768 (Pathogenic, 0 stars; OMIM allele 602860.0011). dbSNP rs1566826570.
- GRCh38 chr15:40,206,279 (SPDI `NC_000015.10:40206279:TTT:TT`; ClinVar displays 40,206,280). GRCh37 chr15:40,498,481.
- gnomAD v4 exomes: AC 3 / AN 1,461,876 (AF 2.05e-6), all East Asian (AC 3 / AN 39,694, AF 7.56e-5). No genome data.
- A filter at AF or popmax AF < 0.001 keeps it. Only a very strict global threshold (for example AF < 1e-5 without popmax) or an "absent from gnomAD" rule would remove it.

### 4b. Ochiai 2014 upstream allele

- A G>A substitution at ss802470619, "44 kb upstream of a *BUB1B* transcription start site", validated by TALEN editing, cosegregating in the families studied and homozygous in one further Japanese infant (PMID 24344301).
- The ss number did not resolve to an rs number or GRCh38 coordinate through E-utilities or the NCBI Variation API. The dbSNP merge table or the paper's supplement (PMC3910577) should give it. Its gnomAD frequency is therefore unknown.
- Qu 2025 (PMID 40555658) reports the same pattern, a coding allele plus an upstream regulatory variant, in two fetuses. Its coordinates were not extracted.

## 5. Failure modes for a singleton WGS pipeline

Constraints: FASTQ only (no BAM/CRAM), GRCh38 no-alt + hs38d1, no `chr` prefix in the supplied VCF ([document 06](06_hackathon_rules.md)).

### 5a. Noncoding second allele filtered out

Default filters keep coding and splice consequences, so an upstream or deep intronic allele like Ochiai 2014, Qu 2025 or Rio Frio 2010 never reaches the candidate list.

Mitigation: call and annotate the whole locus plus the 44 kb upstream window (roughly chr15:40,117,000-40,265,000) with no consequence filter (VEP `--everything`), and rank noncoding candidates with regulatory annotation (ENCODE cCREs, GeneHancer).

### 5b. Mobile-element insertions

Kato 2017 reported an Alu insertion in *BUB1B*. SNV/indel callers (DeepVariant, GATK, Sentieon DNAscope, which produced the supplied VCF) do not call mobile-element insertions.

Mitigation: run a dedicated MEI caller on the realigned reads: RetroSeq (built into nf-core/raredisease), MELT (Gardner 2017), xTea (Chu 2021) or SCRAMble (benchmarked in Wijngaard 2024). MELT and xTea are non-commercial only; RetroSeq has no licence file and SCRAMble's licence is unconfirmed ([document 07](07_data_and_tools.md), section 3).

### 5c. Phase without parents

The organisers will not confirm phase (thread #24). WhatsHap links two variants only through overlapping reads or a chain of heterozygous sites (Patterson 2015). With 300-500 bp fragments (insert size not yet confirmed) and *BUB1B* spanning 60 kb, most candidate pairs cannot be phased from reads.

Mitigation: try WhatsHap and report if it fails; then use population phasing (SHAPEIT5, Hofmeister 2023; or Beagle) against a public panel such as All of Us (thread #22). Check whether the two alleles co-occur in cis in large panels, cite the Ochiai and Qu precedent, and report phase as a probability.

### 5d. Exon-level CNV

SNV callers miss single- or few-exon deletions, and depth-based WGS callers often miss them too.

Mitigation: run CNVkit (WGS mode), GATK-SV, TIDDIT or Manta over the *BUB1B* locus. Manta is non-commercial only. No *BUB1B* exon CNV has been published, so this is a completeness check.

### 5e. Low-VAF mosaic variants

A variant could appear at low allele fraction if it is itself mosaic or if aneuploid clones shift local copy number. Germline callers assume diploid genotypes and can miss it.

Mitigation: relax allele-balance, GQ and DP thresholds for a manual pileup review at *BUB1B* (Exomiser's UDN-derived defaults are 15-85% VAF, GQ >= 20, DP >= 10). Mutect2 in single-sample mode is a heavier option.

### 5f. Variant representation versus the answer key

The same indel can be written several ways. An unnormalised call may not string-match the answer key, and the `chr` prefix differs between the VCF and the template.

Mitigation: run `bcftools norm -f <reference> -m -both` before annotation and before writing the submission (Danecek 2021). Convert contig names separately.

### 5g. Transcript choice

*BUB1B* has at least 8 other protein-coding transcripts. Annotating against one of them can change HGVS and predicted severity.

Mitigation: report consequences on the MANE Select transcript (NM_001211.6 / ENST00000287598), as ClinVar and gnomAD do, and treat other transcripts as secondary.

### 5h. Exomiser exome preset drops noncoding variants

The example exome analysis file (`test-analysis-exome.yml`) removes 5' and 3' UTR, non-coding transcript, intronic, upstream, downstream, intergenic and regulatory variants before prioritisation. The genome preset (Genomiser) keeps them and scores them with REMM, CADD and SpliceAI, with `pathogenicityFilter` set for noncoding targets.

Mitigation: run `--preset genome` with REMM/CADD/SpliceAI enabled whenever a noncoding second allele is possible (Exomiser documentation; Cooperstein 2025).

## 6. Detecting mosaic aneuploidy from blood WGS

Per-chromosome depth shifts and B-allele-fraction deviation from 0.5 at heterozygous sites reveal chromosome-level mosaic gains, losses and copy-neutral LOH. MoChA uses these signals at population scale (Terao 2020; Loh 2020; Zekavat 2020). Detection drops sharply below about 10% cell fraction and depends on the number of informative SNPs; exact figures were not extracted from the papers.

A MoChA-style scan of the patient's WGS could provide direct molecular support for the MVA diagnosis alongside the *BUB1B* genotype. It is optional.

## Open questions

- GRCh38 coordinate and gnomAD frequency of the Ochiai 2014 allele.
- Coordinates of the Qu 2025 variants.
- Overlap of *BUB1B* with GIAB low-mappability regions.
- Whether *BUB1* reads mismap to *BUB1B* kinase-domain exons.
- Sensitivity of mosaic-aneuploidy detection at low clone fractions.

## References

1. Ensembl REST API, `/lookup/id` and `/overlap/id` for ENSG00000156970 and ENST00000287598, GRCh38.
2. NCBI Gene 701, esummary and gene_table. https://eutils.ncbi.nlm.nih.gov/entrez/eutils/.
3. gnomAD GraphQL API (`gnomad_r4`). https://gnomad.broadinstitute.org/api.
4. UCSC Genome Browser REST API, `genomicSuperDups`, hg38. https://api.genome.ucsc.edu/.
5. NCBI E-utilities, db=clinvar, for *BUB1B*, *CEP57*, *TRIP13*, *CENATAC* and *MAD1L1*; ClinVar UID 6768.
6. Ochiai H, et al. TALEN-mediated single-base-pair editing identification of an intergenic mutation upstream of *BUB1B* as causative of PCS (MVA) syndrome. PNAS. 2014. PMID 24344301.
7. Matsuura S, et al. Monoallelic *BUB1B* mutations and defective mitotic-spindle checkpoint in seven families with premature chromatid separation (PCS) syndrome. Am J Med Genet A. 2006. PMID 16411201.
8. Rio Frio T, et al. Homozygous *BUB1B* mutation and susceptibility to gastrointestinal neoplasia. N Engl J Med. 2010. PMID 21190457.
9. Kato M, et al. PCS/MVA syndrome caused by an Alu insertion in the *BUB1B* gene. Hum Genome Var. 2017. PMID 28611924.
10. Qu J, et al. Genetic analysis of two fetuses with mosaic variegated aneuploidy syndrome caused by compound heterozygous variants in *BUB1B* and its upstream regulatory elements. Chin J Med Genet. 2025. PMID 40555658.
11. Gardner EJ, et al. The Mobile Element Locator Tool (MELT). Genome Res. 2017. PMID 28855259.
12. Chuang NT, et al. Mutagenesis of human genomes by endogenous mobile elements on a population scale. Genome Res. 2021. PMID 34772701.
13. Chu C, et al. Comprehensive identification of transposable element insertions using multiple sequencing technologies. Nat Commun. 2021. PMID 34158502.
14. Wijngaard R, et al. Mobile element insertions in rare diseases: a comparative benchmark and reanalysis of 60,000 exome sequences. Eur J Hum Genet. 2024. PMID 37853102.
15. GeneDx. scramble. https://github.com/GeneDx/scramble.
16. Patterson M, et al. WhatsHap: weighted haplotype assembly for future-generation sequencing reads. J Comput Biol. 2015. PMID 25658651.
17. Hofmeister RJ, et al. Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank. Nat Genet. 2023. PMID 37386248.
18. Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021. PMID 33590861.
19. Exomiser documentation. https://exomiser.readthedocs.io/en/latest/.
20. Cooperstein IB, et al. An optimized variant prioritization process for rare disease diagnostics: recommendations for Exomiser and Genomiser. Genome Med. 2025. PMID 41121346.
21. Terao C, et al. Chromosomal alterations among age-related haematopoietic clones in Japan. Nature. 2020. PMID 32581364.
22. Loh PR, et al. Monogenic and polygenic inheritance become instruments for clonal selection. Nature. 2020. PMID 32581363.
23. Zekavat SM, et al. Hematopoietic mosaic chromosomal alterations and risk for infection among 767,891 individuals without hematologic cancer. medRxiv. 2020. PMID 33236019.
