# 01. MVA overview and evaluation criteria

Mosaic variegated aneuploidy (MVA): definition, subtypes, HPO annotations by subtype, surveillance guidance, and how hackathon entries are scored. HPO release 2026-09-01 (annotation file 2026-09-02).

## Summary

- MVA is an autosomal recessive disorder of constitutional mosaic aneuploidy: different cells carry different trisomies or, less often, monosomies. Premature chromatid separation (PCS) is the cytogenetic hallmark of the *BUB1B* form (Nakano 2024).
- OMIM defines seven phenotypes: MVA1 (*BUB1B*), MVA2 (*CEP57*), MVA3 (*TRIP13*), MVA4 (*CENATAC*), Atelis syndrome 1 (*SLF2*, alias MVA5), Atelis syndrome 2 (*SMC5*, alias MVA6) and MVA7 (*MAD1L1*).
- Cancer risk depends on the gene. Nakano 2024: "A convincing increased cancer risk particularly during childhood has been described only for MVA1 and MVA3." MVA1 tumours are mainly Wilms tumour and rhabdomyosarcoma, with MDS, AML and ALL also reported. None of 15 published MVA2 patients had cancer.
- There is no treatment. AACR recommends renal ultrasound every 3 months from birth to age 7 for all MVA, and notes that these children can be very sensitive to radiation and chemotherapy.

## 1. Definition and diagnosis

Nakano 2024 defines MVA as "constitutional mosaic patterns of variegated aneuploidy, including trisomies and, to a lesser extent, monosomies". The aneuploid fraction varies by patient and tissue: 12% on karyotype in one *BUB1B* case (Pavone 2022), about 30-40% of blood cells in the *MAD1L1* case. There is no diagnostic cut-off.

Diagnosis rests on genetic testing plus karyotype of non-malignant metaphases (Nakano 2024, Table 2). PCS is typical of *BUB1B* MVA and was seen in carrier parents in one series (Garcia-Castillo 2008). Patients "generally present with severe microcephaly, growth deficiency, congenital heart defects, skeletal malformations, and intellectual disability".

## 2. Subtypes

| Subtype | Gene | OMIM | Childhood cancer | Distinguishing features (HPO, section 3) |
|---|---|---|---|---|
| MVA1 | *BUB1B* | 257300 | Wilms, rhabdomyosarcoma; also MDS/AML/ALL | IUGR 9/9, microcephaly 9/9, Dandy-Walker, cataract, PCS |
| MVA2 | *CEP57* | 614114 | None in 15 patients | Rhizomelia 3/6, hypothyroidism 3/4, short stature 5/5 |
| MVA3 | *TRIP13* | 617598 | Wilms 6/6 in the discovery series | Nephroblastoma 6/6, PCS 3/3 |
| MVA4 | *CENATAC* | 620153 | Not established | Microcephaly 2/2 (3 annotations only) |
| Atelis 1 (MVA5) | *SLF2* | 620184 | Bone-marrow failure, MDS (Shibata 2026) | Developmental delay 7/7, cafe-au-lait spots 3/7, anaemia |
| Atelis 2 (MVA6) | *SMC5* | 620185 | As Atelis 1 | Anaemia, thrombocytopenia |
| MVA7 | *MAD1L1* | 620189 | Multiple neoplasias in the one reported patient | Inflammation |

All are autosomal recessive. The MVA5 and MVA6 aliases come from OMIM search matches

Differential genes with human evidence but no MVA number: *BUB1* (MCPH30, OMIM 620183), *CDC20* (de novo heterozygous), *MAD2L1BP*, and *TRIM37* (Mulibrey nanism, OMIM 253250), where mosaic aneuploidy was recently reported. Fanconi anaemia has co-occurred with total PCS in one child. PCS does not help separate Cornelia de Lange syndrome ([document 02](02_evidence_matrix.md)).

## 3. HPO annotations by subtype

From `phenotype.hpoa` (2026-09-02), filtered on the seven OMIM IDs. Cells show n/N where the source gives counts, "yes" where the term has no frequency. The table shows 44 of 190 terms: those annotated in at least two subtypes, plus key discriminators.

| HPO ID | Term | MVA1 | MVA2 | MVA3 | MVA4 | Atelis1 | Atelis2 | MVA7 |
|---|---|---|---|---|---|---|---|---|
| HP:0000252 | Microcephaly | 9/9 | 2/4 | 3/6 | 2/2 |  | 2/4 | 1/1 |
| HP:0000347 | Micrognathia | yes | 2/2 |  |  |  | 1/4 | 1/1 |
| HP:0000286 | Epicanthus | yes | 5/6 |  |  |  | 1/4 |  |
| HP:0004322 | Short stature | yes | 5/5 | 1/1 |  |  |  |  |
| HP:0000369 | Low-set ears | yes | 2/2 | 1/6 |  |  | 1/4 |  |
| HP:0001631 | Atrial septal defect | 2/10 | 1/4 |  |  | 1/7 |  |  |
| HP:0001263 | Global developmental delay |  |  | 2/6 |  | 7/7 |  | 1/1 |
| HP:0001250 | Seizure | 2/10 |  | 1/6 |  |  | 1/4 |  |
| HP:0001252 | Hypotonia |  | 1/4 |  |  | 2/7 |  |  |
| HP:0000494 | Downslanted palpebral fissures |  | 1/1 |  |  | 1/7 |  |  |
| HP:0000490 | Deeply set eye |  | 2/5 | 1/6 |  |  |  | 1/1 |
| HP:0000325 | Triangular face |  | 1/1 |  |  |  | 1/4 |  |
| HP:0000276 | Long face |  | 1/1 |  |  | 2/7 | 1/4 |  |
| HP:0000821 | Hypothyroidism |  | 3/4 |  |  | 2/7 |  |  |
| HP:0001511 | Intrauterine growth retardation | 9/9 | 2/4 |  |  |  |  |  |
| HP:0002007 | Frontal bossing |  | 1/1 |  |  |  | 1/4 |  |
| HP:0001510 | Growth delay |  | 6/6 | 2/6 |  |  |  |  |
| HP:0002247 | Duodenal atresia | 1/10 | 1/4 |  |  |  |  |  |
| HP:0000639 | Nystagmus | yes |  | 1/6 |  |  |  | 1/1 |
| HP:0000358 | Posteriorly rotated ears | yes |  |  |  |  |  | 1/1 |
| HP:0000518 | Cataract | yes |  |  |  | 1/7 |  |  |
| HP:0000175 | Cleft palate | yes |  | 1/6 |  |  |  |  |
| HP:0002667 | Nephroblastoma | yes |  | 6/6 |  |  |  |  |
| HP:0003196 | Short nose | yes | 1/1 |  |  |  |  |  |
| HP:0012745 | Short palpebral fissure |  | 1/4 |  |  |  | 1/4 |  |
| HP:0000957 | Cafe-au-lait spot |  | 1/4 |  |  | 3/7 |  |  |
| HP:0000954 | Single transverse palmar crease |  | 1/1 |  |  |  | 1/4 |  |
| HP:0001629 | Ventricular septal defect |  | 2/5 |  |  | 1/7 |  |  |
| HP:0000414 | Bulbous nose |  | 1/1 |  |  |  | 1/4 |  |
| HP:0030084 | Clinodactyly |  | 2/4 |  |  |  | 2/4 |  |
| HP:0005280 | Depressed nasal bridge | yes | 2/4 |  |  |  |  |  |
| HP:0001518 | Small for gestational age | 9/9 | 4/5 |  |  |  |  |  |
| HP:0011800 | Midface retrusion | yes |  |  |  |  |  | 1/1 |
| HP:0006743 | Embryonal rhabdomyosarcoma | 2/10 |  |  |  |  |  | 1/1 |
| HP:0000343 | Long philtrum | yes |  |  |  | 1/7 |  |  |
| HP:0001903 | Anemia |  |  |  |  | 2/7 | 1/4 |  |
| HP:0000218 | High palate |  |  |  |  | 1/7 | 1/4 |  |
| HP:0000448 | Prominent nose |  |  |  |  | 1/7 | 1/4 |  |
| HP:0007018 | Attention deficit hyperactivity disorder |  |  |  |  | 2/7 | 1/4 |  |
| HP:0001873 | Thrombocytopenia |  |  |  |  | 1/7 | 1/4 |  |
| HP:0001642 | Pulmonic stenosis | 1/10 |  |  |  |  | 1/4 |  |
| HP:0200024 | Premature chromatid separation | yes |  | 3/3 |  |  |  |  |
| HP:0001305 | Dandy-Walker malformation | yes |  |  |  |  |  |  |
| HP:0008905 | Rhizomelia |  | 3/6 |  |  |  |  |  |

Uses:
- Phenotype-driven ranking tools (Exomiser; LIRICAL subject to its licence, [document 07](07_data_and_tools.md)) take the patient's HPO terms from the gated phenotype file. This table shows which terms separate the subtypes.
- In the methods report, the patient's profile can be compared against each subtype: nephroblastoma points to MVA1 or MVA3, rhizomelia or hypothyroidism to MVA2.
- The counts come from small curated series and are not population risks.

## 4. Clinical frequencies

Pavone 2022 (PMID 35804254) tabulates 31 published MVA cases, mostly but not all *BUB1B*: growth retardation 29/31, microcephaly 27/31, developmental delay or intellectual disability 24/31, facial dysmorphism 19/31, seizures 14/31, malignancy 12/31, other cerebral anomalies 10/31, Dandy-Walker malformation 8/31, ocular anomalies 6/31. These describe MVA broadly, not *BUB1B* specifically.

## 5. Evaluation criteria

### 5.1 Track 1 scoring

From `evaluation.py` ([document 06](06_hackathon_rules.md)):
- The row matching the answer key earns points by its rank after sorting by EPCR: rank 1 = 100, 2-3 = 50, 4-5 = 25, 6-10 = 10. A row with only one allele of the true pair earns half.
- A compound heterozygote must sit on one row (`_1` and `_2`) to count as a full match.
- F-max is computed per variant across EPCR thresholds, so any extra row at or above the top EPCR adds false positives.
- `finding_type` does not change the score. Secondary rows need an EPCR below the primary row.
- Matching is exact on chrom, pos, ref and alt. Use the `chr` prefix as in the template, and left-normalise REF/ALT against GRCh38 ([document 03](03_variant_spectrum.md), section 5f).

### 5.2 Standards for the Track 1 pipeline

The organisers judge the methods write-up separately and say it "counts just as much" as the score (thread #18). The Track 1 pipeline in this repository aims to:

1. Put the causal variant or pair on one row with the highest EPCR and nothing tied with it.
2. Produce the ranking automatically, end to end from the organiser files, so it can be rerun from this repository.
3. Handle the hard cases explicitly: estimate phase from read-backed phasing, population panels and published precedent, and report it as a probability; search for a noncoding second allele (Genomiser mode, the 44 kb upstream window); state what was done about mobile-element insertions (documents [03](03_variant_spectrum.md) and [07](07_data_and_tools.md)).
4. Compare the patient's HPO profile against all seven subtypes and the differential genes.
5. Keep secondary findings out of this repository while thread #20 is open, with EPCR below the primary row.
6. Record the generative-AI provider, tier and data setting, every tool licence, and the deletion plan ([document 07](07_data_and_tools.md), section 5.4).
7. Optionally, scan per-chromosome depth and allele balance for mosaic aneuploidy, with its sensitivity limits ([document 03](03_variant_spectrum.md), section 6).

### 5.3 Track 2 judging

Scientific Rigor 35, Potential Impact 25, Innovation 25, Scalability 15. Impact asks how much a proposal "could contribute to understanding or diagnosing MVA". The clinical gap is surveillance-only care, embryonal and haematological cancer risk in MVA1, and sensitivity to genotoxic therapy. Documents [04](04_bubr1_biology.md) and [05](05_drug_evidence.md) cover the evidence.

## Open questions

- OMIM alternative titles for the MVA5 and MVA6 aliases.
- The primary sources behind each HPO frequency.

## References

1. Nakano Y, et al. Update on recommendations for cancer screening and surveillance in children with genomic instability disorders. Clin Cancer Res. 2024. PMID 39264246; PMC11705613.
2. Pavone P, et al. Pathogenic correlation between mosaic variegated aneuploidy 1 (MVA1) and a novel *BUB1B* variant: a reappraisal of a severe syndrome. Neurol Sci. 2022. PMID 35804254.
3. Human Phenotype Ontology. phenotype.hpoa 2026-09-02; hp.json 2026-09-01. https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa.
4. Hackathon Space `evaluation.py` and discussion threads #18 and #20 (see [document 06](06_hackathon_rules.md)).
