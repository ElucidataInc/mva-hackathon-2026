# 05. Drug evidence for SAC, BubR1 and aneuploidy biology

Regulatory status, label warnings and published evidence for compounds that could be relevant to Track 2. This document grades evidence. It does not rank or recommend compounds and gives no dosing or clinical advice.

Scope (thread #5): market-approved drugs only; combinations of approved drugs allowed; investigational compounds excluded.

Sources: ChEMBL_37 (release 2026-05-01), openFDA drug labels (266,327 records) and Drugs@FDA (29,352 records), ClinicalTrials.gov API, PubMed.

## Summary

- Of 18 compound classes examined, these approved drugs are in scope: dasatinib, niacin, nicotinamide (as a multivitamin component only), sirolimus, everolimus, metformin, venetoclax and four PARP inhibitors.
- The only BubR1-specific pharmacology is NMN raising BubR1 protein in mice (PMID 24825348). NMN is not an approved drug. The approved NAD+ precursors (niacin, nicotinamide) have never been tested on BubR1.
- In humans, NR raises NAD+ metabolites but has not changed functional endpoints (PMIDs 31412242, 32320006).
- Dasatinib + quercetin clears senescent cells in mice, but the progeroid model was *Ercc1*, not BubR1 (PMID 25754370). Quercetin is a supplement. Human data are small open-label studies in unrelated diseases.
- mTOR inhibitors carry boxed warnings for malignancy and immunosuppression. PARP inhibitors carry MDS/AML warnings and are clastogenic. Both profiles sit badly with a cancer-predisposition syndrome.
- KIF18A, MPS1/TTK and Aurora kinase inhibitors are all investigational. MPS1 and Aurora inhibitors would also weaken the checkpoint that MVA1 already lacks.
- No trial of any of these compounds in MVA or any SAC disorder exists.

## Approved drugs

| Drug | ChEMBL (max phase, first approval) | FDA status | Indications | Paediatric labelling | Key warnings | Link to this biology | Against | Grade |
|---|---|---|---|---|---|---|---|---|
| Dasatinib | CHEMBL1421 (4, 2006) | NDA021986 (SPRYCEL); generics ANDA216547, ANDA216261 | Ph+ CML, Ph+ ALL | From 1 year | Myelosuppression; clastogenic in vitro, negative in rat micronucleus. No boxed warning | Senolytic with quercetin in aged, irradiated and *Ercc1* mice (25754370); reduced senescent cells in human fat and skin (31542391) | Not tested in BubR1 mice; human studies small and uncontrolled (IPF pilot n = 14, 30616998) | Plausible |
| Venetoclax | CHEMBL3137309 (4, 2016) | NDA208573 (VENCLEXTA) | CLL/SLL; AML in combination | Not established | Tumour lysis, neutropenia | BCL-2 family is the senolytic node | BCL-2 selective; lacks the BCL-xL activity behind navitoclax senolysis | Weak |
| Niacin | CHEMBL573 (4, 1982) | NIASPAN (NDA020381) discontinued; generic ER, e.g. ANDA203899 | Hyperlipidaemia; secondary MI prevention | Not established to 16 years | Myopathy with statins; hepatic effects | Approved NAD+ precursor | No BubR1 study; cardiovascular add-on trials reported null | Weak |
| Nicotinamide | CHEMBL1140 (4, 1982) | Only in parenteral multivitamins, e.g. NDA021265 (INFUVITE PEDIATRIC), ANDA210671 | Vitamin maintenance in parenteral nutrition | Paediatric products exist | Those of the multivitamin | Same salvage pathway as NMN | Inhibits sirtuins at high concentration. Approval as a multivitamin component may not count as "market-approved" | Weak |
| Sirolimus | CHEMBL413 (4, 1999) | NDA021110, NDA021083 (RAPAMUNE); topical NDA213478 (HYFTOR) | Renal transplant; LAM; HYFTOR for TSC facial angiofibroma | Transplant from 13 years; HYFTOR from 6 years | Boxed: infection, lymphoma and other malignancy. Lymphoma in an 86-week mouse study | General mTOR ageing rationale | No BubR1 study; malignancy profile | Weak |
| Everolimus | CHEMBL1908360 (4, 2009) | NDA022334, NDA203985 (AFINITOR, DISPERZ); NDA021560 (ZORTRESS) | Oncology; TSC angiomyolipoma, SEGA, seizures | SEGA from 1 year; TSC seizures from 2 years | ZORTRESS boxed: malignancy, infection, graft thrombosis, nephrotoxicity. Myelosuppression, pneumonitis | As sirolimus; precedent of paediatric labelling in a rare genetic disease | No BubR1 study; malignancy profile | Weak |
| Metformin | CHEMBL1431 (4, 1995) | Generics, e.g. ANDA215008 (IR), ANDA090692 (ER) | Type 2 diabetes | IR from 10 years; ER not established | Boxed: lactic acidosis. Genotoxicity negative | General ageing agent only | No SAC or aneuploidy mechanism | Weak |
| PARP inhibitors (olaparib, niraparib, rucaparib, talazoparib) | CHEMBL521686, CHEMBL1094636, CHEMBL1173055, CHEMBL3137320 (all 4) | NDA208558, NDA214876, NDA209115, NDA211651 | HRD/BRCA cancers | Not established | MDS/AML warning, including fatal cases; clastogenic | DNA-repair synthetic lethality | Warnings are the main finding; no MVA study read | Weak |

Everolimus is the clearest example of an mTOR inhibitor labelled for children with a rare genetic disorder (TSC).

## Not approved

| Compound | ChEMBL (max phase) | Status | Note |
|---|---|---|---|
| Quercetin | CHEMBL50 (3) | Supplement | Only relevant as the partner of dasatinib |
| Fisetin | CHEMBL31574 (2) | Supplement | Several trials withdrawn (NCT05505747, NCT04770064) or suspended (NCT05276895) |
| Navitoclax | CHEMBL443684 (3) | Investigational | Thrombocytopenia is dose-limiting; NCT04480086 terminated, NCT06007911 withdrawn |
| Nicotinamide riboside | CHEMBL438497 (3) | Supplement | Raised muscle NAD+ metabolome without changing mitochondrial function (31412242) or insulin sensitivity (32320006) |
| NMN | CHEMBL610238 (none) | Supplement; US status contested | Raised BubR1 in mouse testis after 7 days. The lifespan gain in the same paper came from SIRT2 overexpression, not NMN |
| Resveratrol | CHEMBL165 (3) | Supplement | Sirtuin-activation claims contested; no BubR1 study |
| Spermidine | CHEMBL19612 (3) | Supplement | No BubR1 or SAC link |
| KIF18A inhibitors (sovilnesib, VLS-1488, MEN2501, GenSci122, ATX-295) | Sovilnesib CHEMBL5084301 (1) | Investigational | Aneuploid cancer cells depend on KIF18A (33505028); says nothing about non-malignant mosaic aneuploidy. Trials NCT04293094, NCT06084416, NCT05902988, NCT07226427, NCT06772415, NCT06799065 |
| MPS1/TTK inhibitors (empesertib, BAY-1217389, luvixasertib) | Phase 1-2 | Investigational | Would weaken a checkpoint that is already weak |
| Aurora kinase inhibitors (alisertib, barasertib, tozasertib, danusertib) | Phase 2-3 | Investigational | Aurora B drives error correction; inhibition increases missegregation (not separately cited) |

## Genetic evidence behind the pharmacology

These mouse and human genetics results are what any drug argument rests on:
- BubR1 hypomorphic mice develop progressive aneuploidy, short lifespan, cataract, fat loss and poor wound healing (PMID 15208629).
- Clearing p16-positive cells in these mice delays age-related disease in fat, muscle and eye (PMID 22048312).
- High BubR1 expression reduces tumours, corrects checkpoint and attachment defects and extends lifespan (PMID 23242215).
- Biallelic *BUB1B* mutations cause MVA with microcephaly and childhood cancer (PMID 15475955). A second allele can be an intergenic substitution 44 kb upstream (PMID 24344301).
- Malumbres 2024 reviews MVA causes as SAC and centrosome genes, with tumour risk depending on the gene (PMID 39169218).

## Open questions

1. A ChEMBL identifier for VLS-1488 (trial NCT05902988).
2. The US regulatory status of NMN.
3. Whether nicotinamide inhibits SIRT2 at achievable concentrations, which would undercut the SIRT2 rationale.
4. Whether PARP inhibitors affect aneuploid cells specifically.
5. EMA status and paediatric labelling for each approved drug.

## References

1. Baker DJ, et al. BubR1 insufficiency causes early onset of aging-associated phenotypes and infertility in mice. Nat Genet. 2004. PMID 15208629.
2. Baker DJ, et al. Clearance of p16Ink4a-positive senescent cells delays ageing-associated disorders. Nature. 2011. PMID 22048312.
3. Baker DJ, et al. Increased expression of BubR1 protects against aneuploidy and cancer and extends healthy lifespan. Nat Cell Biol. 2013. PMID 23242215.
4. Zhu Y, et al. The Achilles' heel of senescent cells: from transcriptome to senolytic drugs. Aging Cell. 2015. PMID 25754370.
5. North BJ, et al. SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. EMBO J. 2014. PMID 24825348.
6. Hanks S, et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in *BUB1B*. Nat Genet. 2004. PMID 15475955.
7. Ochiai H, et al. TALEN-mediated single-base-pair editing identification of an intergenic mutation upstream of *BUB1B* as causative of PCS (MVA) syndrome. PNAS. 2014. PMID 24344301.
8. Cohen-Sharir Y, et al. Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition. Nature. 2021. PMID 33505028.
9. Malumbres M, Villarroya-Beltri C. Mosaic variegated aneuploidy in development, ageing and cancer. Nat Rev Genet. 2024. PMID 39169218.
10. Justice JN, et al. Senolytics in idiopathic pulmonary fibrosis: results from a first-in-human, open-label, pilot study. EBioMedicine. 2019. PMID 30616998.
11. Hickson LJ, et al. Senolytics decrease senescent cells in humans: preliminary report from a clinical trial of dasatinib plus quercetin in individuals with diabetic kidney disease. EBioMedicine. 2019. PMID 31542391.
12. Elhassan YS, et al. Nicotinamide riboside augments the aged human skeletal muscle NAD+ metabolome and induces transcriptomic and anti-inflammatory signatures. Cell Rep. 2019. PMID 31412242.
13. Remie CME, et al. Nicotinamide riboside supplementation alters body composition and skeletal muscle acetylcarnitine concentrations in healthy obese humans. Am J Clin Nutr. 2020. PMID 32320006.
14. Wirth M, et al. Correction: Effects of spermidine supplementation on cognition and biomarkers in older adults with subjective cognitive decline (SmartAge). Alzheimers Res Ther. 2022. PMID 35690844.
15. ChEMBL_37 (release 2026-05-01), records for every compound above. https://www.ebi.ac.uk/chembl/.
16. openFDA Drug Labeling. SPL set IDs: 4764f37b-c9e6-4ede-bcc2-8a03b7c521df (SPRYCEL), 58ce42c6-10bd-4bb4-bf89-eae8ed82c339 (dasatinib ANDA216547), 5908cd1a-fc5a-462f-99ed-1d8983e253c9 (sirolimus), edb3ea90-5adc-48ec-99f5-ab963e302f18 (HYFTOR), 2150f73a-179b-4afc-b8ce-67c85cc72f04 (AFINITOR), e082a024-7850-400b-a5c2-2a140612562a (ZORTRESS), d5228e81-9c86-4781-b984-20feaa2ae779 and 7cc02a26-5c22-445b-ad8f-3e7570c143d3 (metformin), 8eff9205-84ac-e69d-e053-2995a90aad4b (niacin ER), 741ff3e3-dc1a-45a6-84e5-2481b27131aa (Lynparza), b7f675e2-159c-490c-b6f4-3f16d9492b7d (ZEJULA), 0295d202-1cfe-7659-e063-6294a90a476e (Rubraca), 13839d4f-6acf-4ffb-a128-79df59319273 (Talzenna), b118a40d-6b56-cee3-10f6-ded821a97018 (Venclexta).
17. Drugs@FDA applications: NDA021986, NDA021110, NDA021083, NDA213478, NDA022334, NDA203985, NDA021560, NDA020357, NDA020381, NDA021265, ANDA210671, NDA208558, NDA214876, NDA208447, NDA209115, NDA211651, NDA208573.
18. ClinicalTrials.gov records: NCT02874989, NCT02848131, NCT04313634, NCT04685590, NCT04063124, NCT04733534, NCT06018467, NCT04946383, NCT05724329, NCT04785300, NCT05838560, NCT05422885, NCT05653258, NCT05506488, NCT06940297, NCT07025226, NCT07270120, NCT07000734, NCT06355037, NCT04994561, NCT04771611, NCT05595499, NCT06431932, NCT07195318, NCT06796374, NCT06819254, NCT07279714, NCT05276895, NCT05025956, NCT05505747, NCT04770064, NCT07676435, NCT05902988, NCT04293094, NCT06084416, NCT07226427, NCT06772415, NCT06799065, NCT04480086, NCT00888108, NCT00878449, NCT05054465, NCT06007911.
19. Open Targets Platform: *BUB1B* ENSG00000156970, *TRIP13* ENSG00000071539, SIRT2 ENSG00000068903, *CEP57* ENSG00000166037.
