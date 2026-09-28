# 02. Evidence matrix: MVA and PCS genetics

Human genetics of mosaic variegated aneuploidy (MVA) and premature chromatid separation (PCS), 1991-2026.

Abbreviations: SAC, spindle assembly checkpoint. MCC, mitotic checkpoint complex. PCD, premature centromere division (older term for PCS). CRC, colorectal cancer. MDS, myelodysplastic syndrome. AR, autosomal recessive.

## Summary

- OMIM numbers seven MVA-type phenotypes across seven genes (*BUB1B*, *CEP57*, *TRIP13*, *CENATAC*, *SLF2*, *SMC5*, *MAD1L1*). Three more genes have human evidence but no MVA number: *BUB1* (biallelic, microcephaly), *CDC20* (de novo dominant-negative, premature ageing) and *MAD2L1BP* (biallelic, mosaic aneuploidy with granulosa-cell tumours).
- ClinGen rates only *BUB1B*-MVA1 and *CEP57*-MVA2, both Definitive. No other MVA gene has a ClinGen curation, so no other gene should be called "ClinGen validated".
- A panel of *BUB1B*, *CEP57* and *TRIP13* alone misses seven genes with human evidence.
- Several papers show that a *BUB1B* second allele can be noncoding: an upstream intergenic substitution (Ochiai 2014), an upstream regulatory variant (Qu 2025), an Alu insertion (Kato 2017) or an intronic transposon (Kawakami 2024). A pipeline that keeps only coding consequences will miss these ([document 03](03_variant_spectrum.md), sections 5a and 5h).
- Cancer spectra differ by gene. *BUB1B* and *TRIP13* carry embryonal tumours (Wilms, rhabdomyosarcoma). *CEP57* has no reported cancer. *SLF2*/*SMC5* (Atelis syndrome) was reclassified in 2026 as an inherited bone-marrow-failure syndrome with early MDS (Shibata 2026). A Track 2 claim that "MVA causes childhood cancer" does not hold across the gene family.
- *MAD2L1BP* (p31comet), *MAD1L1* (MAD1, MVA7) and *MAD2L1* (MAD2) are three different genes. *MAD2L1* itself has no MVA evidence.
- *TRIM37* (Mulibrey nanism) belongs in the phenotype differential: Hakonen 2026 found mosaic aneuploidy in patient fibroblasts and amniocytes.

## A. Gene matrix

| Gene | Phenotype (OMIM) | Mechanism | Variant classes | Tumours | Reported cases | ClinGen | Key PMIDs |
|---|---|---|---|---|---|---|---|
| ***BUB1B*** | MVA1, 257300 (gene 602860) | MCC subunit that inhibits APC/C-Cdc20 and recruits PP2A-B56 to kinetochores. Hypomorph: typically one null allele plus one unstable missense | Nonsense, frameshift, missense near the kinase domain, splice-site, intronic de novo splice site, intergenic substitution 44 kb upstream, Alu insertion | Wilms, embryonal rhabdomyosarcoma; GI adenocarcinoma in adult homozygotes. Monoallelic variants reported in early-onset CRC and prostate cancer | 31 cases compiled by Pavone 2022 (not all genotyped); 5 families (Hanks 2004); 7 Japanese PCS families (Matsuura 2006) | Definitive (Hereditary Cancer GCEP, 2019-11-22) | 15475955; 16411201; 21190457; 24344301; 28611924; 35804254; 39014450 |
| ***CEP57*** | MVA2, 614114 | Centrosomal protein that organises pericentriolar material and binds microtubules. Loss of function; *Cep57* haploinsufficiency impairs tumour suppression in mice | Recurrent c.915_925dup11 (several ethnicities), nonsense, splice-site | None reported in any confirmed case | 12-14 cases in 7-11 families, depending on the review | Definitive (Hereditary Cancer GCEP) | 21552266; 24259107; 30035751; 32861809; 34500087; 35434947; 36635612 |
| ***TRIP13*** | MVA3, 617598 | AAA+ ATPase that, with p31comet, disassembles Mad2 complexes to silence the SAC. Loss of function; SAC impairment as severe as *BUB1B* | Biallelic | Wilms tumour in all 6 original patients; Wilms plus orbital rhabdomyosarcoma in a 2026 case | 1 family (6 individuals); 1 further case | Not curated | 28553959; 42595739 |
| ***CENATAC*** (*CCDC84*) | MVA4, 620153 (gene 620142) | Minor-spliceosome component. Truncating mutations retain AT-AN minor introns in about 100 genes, causing segregation errors. A splicing defect, not a checkpoint defect | Biallelic truncating | Not reported | Exome cohort; family count not in the abstract | Not curated | 34009673 |
| ***MAD1L1*** | MVA7, 620189 | SAC component forming the MAD1-MAD2 template. Biallelic loss of full-length protein | Biallelic | "A dozen of neoplasias" in the one patient | 1 patient, aged 36, 30-40% aneuploid blood cells, systemic inflammation | Not curated | 36322655 |
| ***SLF2*** | Atelis 1, 620184 (alias MVA5) | Recruits SMC5/6 to replication stress. Loss causes cohesion failure and segmented chromosomes ("mosaic variegated hyperploidy") | Biallelic | No solid tumours; bone-marrow failure and early MDS from premature stem-cell ageing | 11 patients across *SLF2* and *SMC5* (Grange 2022) | Not curated | 36333305; 42562919 |
| ***SMC5*** | Atelis 2, 620185 (alias MVA6) | Core ATPase of SMC5/6; same mechanism as *SLF2* | Biallelic; one attenuated case with p.Glu741Ter and p.Asn1022Ser plus near-tetraploidy | As *SLF2* | Shared cohort with *SLF2*; 1 attenuated family (Yang 2025) | Not curated | 36333305; 41374403; 42562919 |
| *BUB1* | MCPH30, 620183 | SAC kinase that recruits BUBR1, Aurora B, SGO1 and TOP2A. Biallelic loss leaves the SAC intact but disrupts centromeric recruitment | Biallelic | None in 2 patients | 2 patients (Carvalhal 2022) | Limited, for colorectal cancer only | 35044816; 23747338 |
| *CDC20* | No OMIM phenotype | APC/C activator. De novo p.Arg286Ser acts as a dominant negative with reduced BUBR1 binding | De novo heterozygous missense | Not reported | 1 patient, premature-ageing phenotype | Not curated | 33094908 |
| *MAD2L1BP* | No phenotype MIM (gene 618136) | p31comet. Nonsense alleles impair SAC silencing and homology-directed repair | Homozygous p.Arg253Ter, p.Arg227Ter | Juvenile granulosa-cell tumours | 2 families, 3 patients | Not curated | 37796616 |

Genes with no human MVA/PCS evidence:
- ***BUB3***: only monoallelic CRC risk jointly with *BUB1* (de Voer 2013). A larger series found such variants rare and without aneuploidy (Mur 2018).
- ***KNL1***: causes microcephaly in mice through progenitor apoptosis ([document 04](04_bubr1_biology.md)); no human MVA case.
- ***CENPE***: tested and excluded in two MVA siblings (Flejter 1998).
- ***TTK*** and ***MAD2L1***: no papers linking them to MVA or PCS.

## B. Paper matrix

83 papers (84 PMIDs; one case was reported twice). T1 = relevant to variant identification, T2 = relevant to mechanism or drugs.

| PMID | First author | Year | Journal | Gene | Study type | N patients/families | Variants as written | Key finding | Relevance |
|---|---|---|---|---|---|---|---|---|---|
| 1809239 | Warburton D | 1991 | Annales de genetique | none identified (pre-*BUB1B*) | Case report | 1 | see paper | First description proposing MVA with microcephaly as "a new human mitotic mutant" | T1 |
| 7985982 | Bitoun P | 1994 | Annales de genetique | none identified (pre-*BUB1B*) | Case report | 1 | see paper | Optic atrophy, microcephaly, mental retardation and MVA as a human mitotic mutation | T1 |
| 9450856 | Flejter WL | 1998 | Am J Med Genet | CENP-E (tested, negative) | Case report + candidate-gene test | 2 siblings | CENP-E not mutated; centromeres functional | Negative result for CENP-E; supports a recessive mitotic mutant of unknown gene at the time | T1 |
| 9677059 | Kajii T | 1998 | Am J Med Genet | none identified (pre-*BUB1B*) | Case report | 1 | see paper | Homozygosity for "total PCS trait" with MVA and congenital abnormalities | T1 |
| 10429359 | Kawame H | 1999 | J Hum Genet | none identified | Case report + review | 1 + review | see paper | Microcephaly, Dandy-Walker malformation and Wilms tumour caused by MVA with PCD | T1 |
| 9916837 | Limwongse C | 1999 | Am J Med Genet | none identified | Case report | 1 | see paper | MVA with embryonal rhabdomyosarcoma | T1 |
| 10877982 | Matsuura S | 2000 | Am J Hum Genet | (localises defect to SAC, pre-*BUB1B*) | Mechanism/cell biology | 1 family | see paper | Total-PCS/MVA chromosomal instability syndrome shown to be defective in the mitotic-spindle checkpoint | T1/T2 |
| 11528703 | Kawame H | 2001 | Ryoikibetsu Shokogun Shirizu | n/a | Japanese-language review | review | n/a | Mosaic variegated aneuploidy, review | T1 |
| 11746029 | Kajii T | 2001 | Am J Med Genet | none identified | Case series | 5 infants | see paper | "Cancer-prone syndrome" of MVA and total PCS, 5 infants | T1 |
| 11932988 | Jacquemont S | 2002 | Am J Med Genet | none identified | Case series/review | n/a | see paper | High risk of malignancy in MVA syndrome | T1 |
| 12116237 | Lane AH | 2002 | Am J Med Genet | none identified | Case report | 1 | see paper | MVA with growth hormone deficiency and congenital heart defects | T1 |
| 12567425 | Hunter A | 2003 | Am J Med Genet A | none identified | Case series/review | n/a | see paper | High risk of malignancy in MVA syndrome (same title as Jacquemont 2002; possibly a republication) | T1 |
| 12567426 | Bower KL | 2003 | Am J Med Genet A | none identified | Case report | 1 | see paper | New case of "apple-peel" intestinal atresia and ocular anomalies with MVA | T1 |
| 14608201 | Furukawa T | 2003 | J Pediatr Hematol Oncol | none identified | Case report | 1 | see paper | Cystic partially differentiated nephroblastoma, embryonal rhabdomyosarcoma and multiple anomalies with MVA/PCD | T1 |
| 14755404 | Chen CP | 2004 | Prenat Diagn | none identified | Case report | 1 | see paper | Prenatal diagnosis of PCD-related MVA | T1 |
| 15098245 | Kajii T | 2004 | Am J Med Genet A | none identified | Case report | 1 | see paper | Prenatal diagnosis of a heterozygous PCS-trait carrier | T1 |
| 15098246 | Kajii T | 2004 | Am J Med Genet A | n/a | Nomenclature note | n/a | n/a | Premature chromatid separation (PCS) versus premature centromere division (PCD): the two terms used in this literature are argued to describe the same phenomenon | T1 |
| 15108198 | Ikeuchi T | 2004 | Am J Med Genet A | pre-*BUB1B* | Cytogenetic study | PCS-trait carriers and controls | see paper | Induction of PCS in individuals with the PCS trait and in normal controls | T1/T2 |
| 15475955 | Hanks S | 2004 | Nat Genet | *BUB1B* | Discovery genetics | 5 families | truncating and missense | Biallelic *BUB1B* mutations cause constitutional aneuploidy and cancer predisposition (MVA1 discovery) | T1/T2 |
| 15767571 | Rao CV | 2005 | Proc Natl Acad Sci USA | *BUB1B* (mouse) | Mouse model | mouse | n/a | BubR1+/-;ApcMin/+ compound mutant mice: colonic tumorigenesis linked to premature sister-chromatid separation and genomic instability | T2 |
| 16059936 | Callier P | 2005 | Am J Med Genet A | *BUB1B*-era case | Case report | 1 | see paper | Microcephaly is not mandatory for the diagnosis of MVA syndrome | T1 |
| 16174301 | Corona-Rivera A | 2005 | BMC Med Genet | PCD | Cell biology | not stated in title | see paper | Cell cycle and centromere FISH studies in premature centromere division | T1/T2 |
| 16182441 | Hanks S | 2006 | Cancer Lett | *BUB1B* | Genetics | childhood cancer cohort | see paper | Comparative genomic hybridisation and *BUB1B* mutation analysis in childhood cancers associated with MVA | T1/T2 |
| 16411201 | Matsuura S | 2006 | Am J Med Genet A | *BUB1B* | Genetics | 7 families | monoallelic mutations, 1833delT recurrent in 4 families, plus splice-site, nonsense, missense | Monoallelic *BUB1B* mutations and defective mitotic-spindle checkpoint in 7 PCS families; no coding second allele found | T1/T2 |
| 17632782 | Micale MA | 2007 | Am J Med Genet A | *BUB1B*-era case | Case report | 1 | see paper | MVA without microcephaly: implications for cytogenetic diagnosis | T1 |
| 17981217 | Petković I | 2007 | Cancer Genet Cytogenet | PCS-era case | Case report | 2 (mother and son) | see paper | PCS in a woman with cervical carcinoma in situ and in her son with keratoacanthoma | T1 |
| 18548531 | García-Castillo H | 2008 | Am J Med Genet A | *BUB1B* (screened) | Genetics/clinical | 21 of 35 screened | see paper | Delineation of 3 clinical subtypes by *BUB1B* genotype; monoallelic group most severely affected in this series (Dandy-Walker 7/8, cataract 6/6, Wilms 7/8) | T1 |
| 18932004 | Bohers E | 2008 | Hum Genet | *BUB1B* | Cell biology | patient cells | see paper | Gradual reduction of BUBR1 protein levels results first in premature sister-chromatid separation, then in aneuploidy | T1/T2 |
| 19373780 / 19856400 | Sari N / Alter BP | 2009 / 2010 | Pediatr Blood Cancer | PCS + Fanconi anemia | Case report (same case, two reports) | 1 | see paper | Wilms tumour, AML and medulloblastoma in a child with cancer-prone total-PCS syndrome and co-occurring Fanconi anaemia | T1 |
| 19503101 | Izumi H | 2009 | Oncogene | *BUB1B* (mechanism) | Cell biology | cell lines | n/a | BubR1 localises to centrosomes and suppresses centrosome amplification via regulating Plk1 activity in interphase | T2 |
| 19690971 | Castronovo P | 2009 | Chromosome Res | PCS (differential) | Cytogenetic study | Cornelia de Lange syndrome cohort | n/a | PCS is not a useful diagnostic marker for Cornelia de Lange syndrome | T1 |
| 20516114 | Suijkerbuijk SJ | 2010 | Cancer Res | *BUB1B* | Cell biology, patient cells | patient-derived lines | see paper | Molecular causes for BUBR1 dysfunction in MVA1: truncating alleles give no transcript, missense alleles are unstable; ectopic BubR1 rescues checkpoint activity | T1/T2 |
| 21190457 | Rio Frio T | 2010 | N Engl J Med | *BUB1B* | Case report | 1 | homozygous intronic c.2386-11A>G creating a de novo splice site | Homozygous *BUB1B* mutation with adult-onset gastrointestinal neoplasia, no classical MVA1 syndrome | T1/T2 |
| 21389084 | Miyamoto T | 2011 | Hum Mol Genet | *BUB1B* | Cell biology | patient cells / model | n/a | Insufficiency of BUBR1 causes impaired ciliogenesis in vertebrates | T2 |
| 21449797 | de Voer RM | 2011 | N Engl J Med | *BUB1*/*BUB3* (correspondence) | Correspondence | n/a | n/a | Spindle-assembly checkpoint and gastrointestinal cancer (letter responding to Rio Frio 2010) | T1/T2 |
| 21552266 | Snape K | 2011 | Nat Genet | *CEP57* | Discovery genetics | initial cohort | biallelic loss-of-function | Mutations in *CEP57* cause MVA syndrome (MVA2 discovery) | T1/T2 |
| 23300461 | Wijshake T | 2012 | PLoS Genet | *BUB1B* (mouse) | Mouse model | mouse | n/a | Reduced life- and healthspan in mice carrying a mono-allelic BubR1 MVA mutation | T2 |
| 23602569 | Baker DJ | 2013 | Cell Rep | *BUB1B* (mouse) | Mouse model | mouse | n/a | p21 both attenuates and drives senescence and ageing in BubR1 progeroid mice | T2 |
| 23747338 | de Voer RM | 2013 | Gastroenterology | *BUB1*, *BUB3* | Genetics, cohort | 208 CRC patients, mutations in 2.9% | haploinsufficiency/heterozygous mutations | Germline *BUB1*/*BUB3* mutations are risk factors for early-onset/familial colorectal cancer and cause variegated aneuploidy with dysmorphic features in carriers | T1/T2 |
| 23789096 | Xu P | 2013 | Biol Open | *BUB1B* (mechanism) | Cell biology | n/a | n/a | BUBR1 recruits PP2A via the B56 family of targeting subunits to promote chromosome congression | T2 |
| 23916859 | Akasaka N | 2013 | Pediatr Neurol | *BUB1B*-era case | Case report | 1 | see paper | Refractory infantile spasms associated with MVA syndrome | T1 |
| 24259107 | Pinson L | 2014 | Am J Med Genet A | *CEP57* | Case report | 1 (5th proband reported) | c.915_925dup11, p.Leu309ProfsX9 | *CEP57* mutation in a girl with MVA; rhizomelic limb shortening and absence of tumour proposed as *CEP57*-suggestive features | T1 |
| 24344301 | Ochiai H | 2014 | Proc Natl Acad Sci USA | *BUB1B* | Functional genetics | families with monoallelic coding findings | intergenic substitution upstream of *BUB1B* TSS | TALEN-validated intergenic upstream mutation as the "missing" second allele in PCS/MVA families | T1/T2 |
| 24981203 | Schmid M | 2014 | Chromosome Res | *BUB1B* (mouse + human) | Mouse model + human note | mouse; human pregnancy-loss note | n/a | MVA in mouse BubR1-deficient embryos and a note on human pregnancy loss | T1/T2 |
| 25696020 | Cho CH | 2015 | Ann Clin Lab Sci | not stated | Case report + review | 1 fetus | see paper | Fetus with mosaic autosomal variegated aneuploidies, literature review | T1 |
| 26309087 | Miyamoto T | 2015 | Oncotarget | *BUB1B* | Cell biology | patient cells | n/a | Ciliopathy in PCS (MVA) syndrome, extending the ciliogenesis finding from 2011 | T2 |
| 27239782 | Hahn MM | 2016 | Genes Chromosomes Cancer | *BUB1B* | Genetics, cohort | 192 early-onset CRC patients | p.Glu390del, p.Cys945Tyr (heterozygous) | Germline *BUB1B* mutations are not enriched in early-onset CRC; no biallelic carriers found; the two heterozygous variants did not affect protein level or *BUB1*/*BUB3* binding | T1/T2 |
| 27663603 | Shimada K | 2017 | Pediatr Blood Cancer | PCS-era case | Case report | 1 | see paper | Wilms tumour accompanied by premature chromatid separation | T1 |
| 27931980 | Chaker F | 2017 | Ann Endocrinol (Paris) | not stated | Case report | 1 | see paper | Polycystic ovary syndrome proposed as a new phenotype in MVA syndrome | T1 |
| 28383136 | Yang Z | 2017 | Aging Cell | *BUB1B* (mouse) | Mouse model | mouse | n/a | Age-related decline in BubR1 impairs adult hippocampal neurogenesis | T2 |
| 28553959 | Yost S | 2017 | Nat Genet | *TRIP13* | Discovery genetics | 6 individuals | see paper | Biallelic *TRIP13* mutations predispose to Wilms tumour and chromosome missegregation (MVA3 discovery); contrasted with *CEP57* (no embryonal tumours) | T1/T2 |
| 28611924 | Kato M | 2017 | Hum Genome Var | *BUB1B* | Case report | 1 | Alu insertion in *BUB1B* | PCS/MVA syndrome caused by an Alu insertion in *BUB1B* | T1 |
| 29448935 | Mur P | 2018 | Mol Cancer | *BUB1*, *BUB3* | Genetics, cohort | 456 familial CRC + 88 polyposis | one splice-site (c.1965-1G>A) functionally validated, others VUS | Germline *BUB1*/*BUB3* mutations are infrequent in familial CRC/polyposis; carriers had neither variegated aneuploidy nor dysmorphic traits | T1/T2 |
| 29673003 | Yamaguchi T | 2018 | J Obstet Gynaecol Res | PCS/MVA-era case | Case report | 1 | see paper | Prenatal diagnosis of PCS/MVA syndrome | T1 |
| 30035751 | Aziz K | 2018 | J Clin Invest | *CEP57* (mouse) | Mouse model | mouse | n/a | MVA-syndrome mutation or haploinsufficiency in *Cep57* impairs tumour suppression | T1/T2 |
| 30147898 | Brightman DS | 2018 | Clin Case Rep | *CEP57* | Case report | 1 | see paper | MVA syndrome caused by a *CEP57* mutation, diagnosed by whole-exome sequencing | T1 |
| 30512160 | He T | 2018/2019 | Chin J Med Genet | not stated | Case report (Chinese) | 1 child | see paper | Clinical features and genetic analysis of a child with MVA syndrome | T1 |
| 30010053 | De la Torre-García O | 2019 | Eur J Med Genet | *CEP57* | Case report | 1 | homozygous c.915_925dupCAATGTTCAGC | *CEP57* patient with rhizomelic shortening of upper and lower limbs and narrow thorax; second Mexican family with the same duplication, suggesting a founder effect | T1 |
| 30668728 | Simmons AJ | 2019 | Hum Mol Genet | *BUB1B* (mouse) | Mouse model | mouse | n/a | Nearly complete deletion of BubR1 causes microcephaly through shortened mitosis and massive cell death | T2 |
| 30804344 | Watanabe K | 2019 | Nat Commun | *CEP57* (mechanism) | Cell biology | n/a | n/a | The Cep57-pericentrin module organises pericentriolar material expansion and centriole engagement | T2 |
| 31053147 | Laberko A | 2019 | Orphanet J Rare Dis | *BUB1B* | Case report | 1 | see paper | Haematopoietic stem cell transplantation in a patient with MVA1 (*BUB1B*) | T1 |
| 31081598 | Ochiai K | 2019 | Pediatr Blood Cancer | *BUB1B*-era case | Case report | 1 | see paper | Long-term remission of bilateral Wilms tumours in PCS/MVA syndrome via bilateral nephrectomy and peritoneal dialysis | T1 |
| 31184400 | Nishitani-Isa M | 2019 | Pediatr Int | PCS/MVA-era case | Case report | 1 | see paper | Rhabdomyosarcoma with PCS-MVA syndrome treated with reduced-intensity chemotherapy | T1 |
| 31738183 | Sieben CJ | 2020 | J Clin Invest | *BUB1B* (mouse) | Mouse model | mouse allelic series | n/a | BubR1 allelic effects drive phenotypic heterogeneity in MVA progeria syndrome (PMID 33136097 appears to index the same article) | T1/T2 |
| 31943948 | Pezzani L | 2020 | Mol Genet Genomic Med | *CEP57* | Case report | 1 | see paper | Double homozygosity in *CEP57* and *DYNC2H1* in one infant by trio-WES; a dual diagnosis in a consanguineous family | T1 |
| 31963583 | Akutsu SN | 2020 | Cells | multiple | Review | review | n/a | Applications of genome-editing technology in research on chromosome-aneuploidy disorders (review, title only) | T2 |
| 32861809 | Dery T | 2020 | Eur J Med Genet | *CEP57* | Case report | 2 (adult brothers) | homozygous c.915_925dup11 | Follow-up of two adult brothers expands the MVA2 phenotype into adulthood; proposes a FoSTeS mechanism for the recurrent duplication; states only 7 *CEP57* patients reported to date at time of writing | T1 |
| 32884756 | Lin SM | 2020 | Clin Case Rep | *BUB1B* | Case report | 1 | see paper | Prenatal diagnosis and long-term follow-up of a *BUB1B*-related MVA patient in Hong Kong | T1 |
| 33094908 | Fujita H | 2020 | Aging Cell | *CDC20* | Discovery genetics | 1 | de novo heterozygous c.856C>A, p.Arg286Ser | Dominant-negative *CDC20* mutation causes a premature-aging syndrome with random chromosome-number instability and SAC failure | T1/T2 |
| 33892776 | Bloom CR | 2021 | Cell Biosci | *BUB1B* (review) | Review | review | n/a | Physiological relevance of post-translational regulation of BubR1 (review, title only) | T2 |
| 34009673 | de Wolf B | 2021 | EMBO J | *CENATAC* | Discovery genetics | exome cohort (n not stated) | biallelic truncating | Chromosomal instability by mutations in the novel minor-spliceosome component *CENATAC* (MVA4 discovery) | T1/T2 |
| 34500087 | Santos-Simarro F | 2021 | Eur J Med Genet | *CEP57* | Case report + review | 2 new + review of 12 total | homozygous c.915_925dup11 | MVA2 phenotype review: growth retardation, facial dysmorphism, endocrine/cardiovascular/skeletal anomalies; explicitly, malignancies not reported in the reviewed cohort | T1 |
| 35044816 | Carvalhal S | 2022 | Sci Adv | *BUB1* | Discovery genetics | 2 patients | biallelic (specific HGVS not restated) | Biallelic *BUB1* mutations cause microcephaly, developmental delay and variable cohesion/segregation defects: *BUB1* loss and *BUB1* kinase-dead act through different downstream defects | T1/T2 |
| 35434947 | Feng B | 2022 | Mol Genet Genomic Med | *CEP57* | Case report | 1 | homozygous NM_014679.4:c.312T>G, p.Tyr104* | Novel *CEP57* nonsense variant; growth-hormone treatment reported ineffective in this patient | T1 |
| 35804254 | Pavone P | 2022 | Neurol Sci | *BUB1B* | Case report + systematic review | 1 new + review of 31 total | see section C | Reappraisal of MVA1 severity; tabulates phenotype prevalence across 31 published cases (Table 1, read in full) | T1 |
| 36322655 | Villarroya-Beltri C | 2022 | Sci Adv | *MAD1L1* | Discovery genetics | 1 | biallelic (specific HGVS not restated) | Biallelic *MAD1L1* mutations cause a new MVA variant with systemic inflammation and unprecedented tumour susceptibility (MVA7 discovery) | T1/T2 |
| 36333305 | Grange LJ | 2022 | Nat Commun | *SLF2*, *SMC5* | Discovery genetics | 11 patients | biallelic, both genes (specific HGVS not restated) | Pathogenic *SLF2*/*SMC5* variants cause segmented chromosomes and mosaic variegated hyperploidy, named Atelis syndrome | T1/T2 |
| 36588761 | Hübner CT | 2022 | Mol Syndromol | not stated | Case report | 1 family | see paper | MVA syndrome and Noonan syndrome co-occurring in the same family | T1 |
| 36635612 | Langeh N | 2023 | Clin Genet | *CEP57* | Case report | 2 siblings | homozygous splice-site c.382+2T>C, exon 3 skipping | First reported *CEP57* splice-site variant causing MVA2; states only 11 families and 5 pathogenic variants published to date | T1 |
| 37796616 | Abdel-Salam GMH | 2023 | JCI Insight | *MAD2L1BP* | Discovery genetics | 2 families / 3 patients | homozygous nonsense R253* and R227* | Biallelic *MAD2L1BP* (p31comet) mutation causes mosaic aneuploidy and juvenile granulosa-cell tumours | T1/T2 |
| 37900274 | Sterling NA | 2023 | Front Cell Dev Biol | *BUB1B* (mechanism) | Cell biology | n/a | n/a | p53-independent pathogenic mechanisms contribute to BubR1 microcephaly | T2 |
| 37981762 | Guo J | 2024 | HGG Adv | *CEP192* (candidate) | Case report | 1 | see paper | MVA syndrome with tetraploidy and male-infertility predisposition linked to mutant *CEP192*; *CEP192* is not yet an OMIM-numbered MVA gene | T1 |
| 38102195 | Kawakami R | 2024 | J Hum Genet | *BUB1B* | Case report + RNA-seq | 1 | intronic transposon insertion causing aberrant splicing | RNA sequencing and targeted long-read sequencing reveal an intronic transposon insertion in *BUB1B* causing aberrant splicing | T1 |
| 38664248 | Esposito M | 2024 | Int J Legal Med | *BUB1B* | Forensic case report | 1 | see paper | Sudden unexpected postnatal collapse and a *BUB1B* mutation; first forensic report | T1 |
| 38665733 | Taguchi S | 2024 | Cureus | PCS | Case report | 1 couple | see paper | Genetic counselling for an infertile couple with PCS syndrome | T1 |
| 38699879 | Sukla S | 2024 | Proteins | *CEP57* (structure) | Structural biology | n/a | n/a | Crystal structure of the human Cep57 C-terminal domain reveals a leucine zipper and a candidate microtubule-binding region | T2 |
| 38770464 | Kawamura S | 2024 | Cureus | PCS | Case report | 1 | see paper | Premature chromatid separation trait found during male-infertility work-up | T1 |
| 38857398 | Yeh HW | 2024 | Proc Natl Acad Sci USA | *CEP57* (mechanism) | Structural/cell biology | n/a | n/a | Cep57 regulates human centrosomes through multivalent interactions | T2 |
| 39014450 | Silva MP | 2024 | J Biomed Sci | *BUB1B* | Genetics, cohort | 462 early-onset/familial prostate cancer + 1,416 other hereditary-cancer patients | rare heterozygous variants in ~1.9% of prostate cases (specific HGVS not restated) | Monoallelic *BUB1B* germline variants contribute to prostate-cancer predisposition by triggering chromosomal instability and premature chromatid separation, and drive Taxol resistance | T1/T2 |
| 39169218 | Malumbres M | 2024 | Nat Rev Genet | multi-gene | Review | review | n/a | Mosaic variegated aneuploidy in development, ageing and cancer | T1/T2 |
| 39264246 | Nakano Y | 2024 | Clin Cancer Res | multi-gene (surveillance guideline) | Review/guideline | n/a | n/a | Update on recommendations for cancer screening and surveillance in children with genomic instability disorders | T1/T2 |
| 39392177 | Frattini A | 2025 | Am J Med Genet A | not stated | Case report + review | 1 + review | see paper | A distinctive type of MVA: case report and review of the literature | T1 |
| 39490209 | Ferreira A | 2024 | Stem Cell Res | *BUB1* | iPSC resource | 2 isogenic lines from 1 patient | compound heterozygous | Generation and characterisation of 2 isogenic iPSC lines from a young female with microcephaly and compound heterozygous *BUB1* mutation | T1/T2 |
| 40055864 | Pun R | 2025 | J Am Heart Assoc | *BUB1B* (mouse) | Mouse model | mouse | n/a | BubR1 controls heart development by promoting expression of cardiogenesis regulators | T2 |
| 40441073 | Milán M | 2025 | Curr Opin Cell Biol | multi-gene (review) | Review | review | n/a | Chromosomal instability in development and disease: beyond cancer evolution (general CIN review, not MVA-specific) | T2 |
| 40555658 | Qu J | 2025 | Chin J Med Genet | *BUB1B* | Case report | 2 fetuses | compound heterozygous *BUB1B* coding plus upstream regulatory-element variant | Genetic analysis of 2 fetuses with MVA caused by compound heterozygous *BUB1B* coding and upstream regulatory variants, with literature review | T1 |
| 41374403 | Yang Y | 2025 | Diagnostics (Basel) | *SMC5* | Case report | 1 | compound heterozygous c.2221G>T (p.Glu741Ter) and c.3065A>G (p.Asn1022Ser) (PMC12691382) | Recessive *SMC5* variants in a family with near-tetraploidy/MVA; attenuated phenotype without the severe Atelis-syndrome features | T1 |
| 41700350 | Viudes CP | 2026 | J Pediatr Endocrinol Metab | *CEP57* | Case report | 1 | see paper | Novel *CEP57* mutation in MVA2 | T1 |
| 41702694 | Hakonen AH | 2026 | J Med Genet | *TRIM37* (differential) | Case series | 2 index siblings + 10 further Mulibrey patients | see paper | MVA as a novel feature in Mulibrey nanism/*TRIM37* patients: 7-36% aneuploid metaphases in index cases, low-level abnormality in 1 of 10 further patients | T1 |
| 41820377 | González-Blanco A | 2026 | Nat Commun | multi-gene (mechanism) | Cell/mouse biology | n/a | n/a | Proteostasis failure and mitochondrial dysfunction contribute to chromosomal-instability-induced microcephaly (general CIN-microcephaly mechanism, not confirmed *BUB1B*-specific) | T2 |
| 42434306 | Wei TY | 2026 | Front Endocrinol | *BUB1B* | Case report | 2 novel heterozygous variants | see paper | Functional and clinical evidence for 2 novel heterozygous *BUB1B* variants and their value in precision genetic counselling for recurrent pregnancy loss | T1 |
| 42562919 | Shibata S | 2026 | Leukemia | *SLF2*, *SMC5* | Mechanism + clinical follow-up | patient-derived iPSC/HPCs; longitudinal cohort not sized in abstract | n/a | *SLF2*/*SMC5* dysfunction drives HSC ageing and predisposes to MDS, reframing Atelis syndrome as an inherited bone-marrow-failure syndrome | T1/T2 |
| 42595739 | Vempuluru VS | 2026 | Orbit (Amst) | *TRIP13* | Case report | 1 | c.998_999delCT, p.Ser333Leufs*30 | Sequential Wilms tumour and orbital rhabdomyosarcoma in a child with MVA3 (*TRIP13*) | T1 |

## C. *BUB1B* genotype and phenotype

Pavone 2022 Table 1 (PMC9616775) lists 31 published MVA cases by feature, not by genotype. Many pre-2004 cases were never tested for *BUB1B*, so these rates describe MVA broadly.

| Feature | n/31 (%) |
|---|---|
| Pre/postnatal growth retardation | 29 (93.5) |
| Microcephaly | 27 (87.0) |
| Developmental delay / intellectual disability | 24 (77.4) |
| Facial dysmorphism | 19 (61.2) |
| Epileptic seizures | 14 (45.1) |
| Any malignancy | 12 (38.7) |
| Other cerebral anomalies | 10 (32.2) |
| Dandy-Walker malformation | 8 (25.8) |
| Ocular anomalies | 6 (19.3) |
| Endocrine, genital or renal anomalies | 4 each (12.9) |
| Female | 14 (45.1) |

Cataract appears only as free text, in at least 3 rows, so it has no rate.

Published *BUB1B* genotypes:

| Source | N | Allele 1 | Allele 2 | Phenotype notes |
|---|---|---|---|---|
| Matsuura 2006 | 7 families | 1833delT (4 families); splice-site, nonsense or missense (1 each) | No coding second allele found | |
| Ochiai 2014 | Families with one coding allele | Coding allele | Intergenic substitution 44 kb upstream of the TSS, validated by TALEN editing | |
| Rio Frio 2010 | 1 adult | c.2386-11A>G | Same (homozygous) | Ampullary adenocarcinoma at 34, then colon and gastric adenocarcinomas; no MVA1 features |
| Kato 2017 | 1 | Alu insertion | Not in abstract | |
| Qu 2025 | 2 fetuses | Coding variant | Upstream regulatory variant | Prenatal |
| Silva 2024 | ~1.9% of 462 prostate-cancer patients | Rare heterozygous variant | None | Adult prostate cancer; Taxol resistance in cells |

A full genotype table would need the ~20 source papers behind Pavone's table.

## D. Other points for Track 1 and Track 2

- *CDC20* disease is de novo and dominant negative. A pipeline that assumes recessive inheritance for every MVA-family gene will miss it.
- Biallelic *BUB1* loss causes microcephaly (Carvalhal 2022). Monoallelic *BUB1*/*BUB3* variants are a weak, poorly replicated CRC risk (de Voer 2013 positive, Mur 2018 negative). Keep the two findings apart.
- *CENATAC* has no drug-repurposing literature and a different mechanism from every other MVA gene.

## Open questions

- Gene-level OMIM numbers for *CEP57*, *TRIP13*, *MAD1L1*, *SLF2*, *SMC5*, *BUB1* and *CDC20*.
- HGVS strings for the *TRIP13*, *CENATAC*, *MAD1L1*, *SLF2*/*SMC5* and *BUB1* discovery papers and for rows marked "see paper".
- Whether PMIDs 31738183 and 33136097 (Sieben) are one article, and whether Hunter 2003 republishes Jacquemont 2002.
- Whether *CEP192* (Guo 2024) has an OMIM phenotype.

## References

The paper matrix (section B) is the bibliography: each row gives first author, year, journal and PMID. Other sources:

- NCBI E-utilities (`esearch`, `esummary`, `efetch`; db=pubmed, pmc, omim). https://eutils.ncbi.nlm.nih.gov/entrez/eutils/.
- Pavone 2022 full text, PMC9616775.
- ClinGen Gene-Disease Validity, search.clinicalgenome.org/kb/genes/<SYMBOL>, for *BUB1B*, *CEP57*, *TRIP13*, *MAD1L1*, *SLF2*, *SMC5*, *BUB1*, *BUB3*, *CDC20*, *MAD2L1BP*, *TRIM37* and *CENPE*.
- Monarch Initiative API v3. https://api-v3.monarchinitiative.org.
