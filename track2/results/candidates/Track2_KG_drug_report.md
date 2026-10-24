# Track 2: drugs acting on the MVA marker genes and Track 1 genes (Polly KG)

Prepared 23 Oct 2026. Source: Polly KG `1668691496_pollybkg` (PollyBKG v2), read-only queries. Companion table: `Track2_drug_table_16genes.tsv` (27 rows).

## 1. Scope

**16 genes.** The 11 marker genes (CDC20, BUB1, TRIP13, MAD2L1BP, MAD1L1, SMC5, SLF2, CEP57, CENATAC, BUB1B, TRIM37; confirmed from the team gene file and the 8 Oct meeting) plus the final Track 1 genes in the team drug sheet (ACLY, DIDO1, ARL9, LZTR1, GNRHR).

**Partner genes** checked for indirect drugs: AURKB, PLK1, TTK (MPS1), PLK4, MTOR, NAMPT, SIRT2, BUB3, MAD2L1, SMC6, NSMCE2.

**Questions asked for every gene:**
1. Which drugs act on it, directly or through a partner?
2. What is the mechanism and the target, and how far along is the drug (phase, approval)?
3. Is the drug's target allele similar to the proband's allele?
4. Does the drug restore, compensate, exploit or worsen the proband's loss of function?

Queries used gene names and variant types only, never patient coordinates.

## 2. How the KG stores drug information

| What | KG structure | Key fields used |
|---|---|---|
| Compound acts on protein | `(gene)-[:encodes]->(protein)<-[:binds_target or :inhibits_target]-(molecule)` (ChEMBL) | molecule: `name`, `max_phase`, `first_approval`, `black_box_warning`; edge: `action_type`, `pchembl_value` |
| Clinical drug for a target | `(gene)-[:gene_participates_in]->(has_drug_target)` (Open Targets clinical precedence) | `drug_id`, `drug_from_source`, `clinical_stage`, `direction_on_target` (LoF = inhibits, GoF = activates), `disease_from_source`, `clinical_report_id` |
| Synthetic lethality / therapeutic target | `(gene)-[:gene_participates_in]->(synthetic_lethal_target_in_disease or therapeutic_target_for_disease)-(disease)` | `evidence_text`, `pmids`, `claim_polarity` |
| Drug resistance / predictive biomarker | `confers_drug_resistance_in_disease`, `predictive_biomarker_for` | disease, PMIDs |

**Not in the KG for these genes:** pharmacogenomic links (variant to drug response). None of the 34 association types attached to the 16 genes is a drug-response type, so the "similar allele" question was answered by allele-type reasoning (Section 4.4).

## 3. Headline findings

1. **No drug acts directly on 14 of the 16 genes in a useful direction.** Only ACLY and GNRHR have clinical drugs (Open Targets clinical-precedence edges). BUB1 has one clinical binder (AT-9283, an inhibitor). The other marker genes have no clinical drug edge.
2. **Almost every drug in the checkpoint neighbourhood is an inhibitor.** Aurora B, PLK1, MPS1/TTK, PLK4, BUB1 and NAMPT inhibitors weaken the checkpoint or lower NAD+. For a patient whose BUBR1 is already low, these are **contraindications**, not candidates.
3. **The team's BUB1B candidates stand.** The KG adds no new direct BUBR1-raising drug. Niacin (NAD+/SIRT2 axis) and gentamicin (read-through of the UGA stop allele) remain the only restorative ideas; mTOR inhibitors and senolytics act on consequences.
4. **BUB1B cancer-context paper (cited under avoid):** BUBR1-depleted MCF-7 cells died more with cisplatin plus the research compound CXI-benzo-84 (PMID 30468712). Not a lead: research compound, cancer cells only.
5. **Correction to the team sheet:** ACLY does have KG drug evidence. Bempedoic acid has 22 Open Targets clinical records (approval stage, 12 registered trials plus FDA and EMA label records, cardiovascular and lipid indications). The verdict stays "exclude": the drug inhibits ACLY, the same direction as the proband's predicted damaging missense, and ACLY is not an MVA gene.
6. **Allele-type matches:** BUB1B p.Leu737Ter (UGA) and LZTR1 p.Tyr748Ter (UAG) are both stop-gains, so read-through agents apply in principle (UGA is the more responsive stop). GNRHR p.Arg139His belongs to a receptor class where pharmacoperones rescue many misfolded missense alleles; whether R139H is a folding defect needs checking. Both LZTR1 and GNRHR are secondary findings.

## 4. Results by step

### 4.1 Direct compounds (ChEMBL edges)

| Gene | Molecules in KG | Clinical (max_phase ≥ 1) | Named / clinical compounds |
|---|---|---|---|
| BUB1 | 709 | 1 | AT-9283 (CHEMBL495727), binds BUB1 at pChEMBL 7.03, phase 2, multi-kinase inhibitor |
| GNRHR | 704 | 12 | Agonists (gonadorelin, leuprolide) and antagonists (degarelix, relugolix, elagolix, linzagolix, cetrorelix, ganirelix, abarelix with black box warning); sufugolix (TAK-013) and acyline at phase 2 |
| ACLY | 2 | 0 | Assay compounds only (bempedoic acid sits on the Open Targets edge instead) |
| TRIP13 | 1 | 0 | One assay compound |
| BUB1B, CDC20, CENATAC, CEP57, MAD1L1, MAD2L1BP, SLF2, SMC5, TRIM37, DIDO1, ARL9, LZTR1 | 0 | 0 | None |

### 4.2 Clinical drugs (Open Targets clinical precedence, `has_drug_target`)

| Gene | Records | Drugs (stage, direction) |
|---|---|---|
| ACLY | 22 | Bempedoic acid / ETC-1002: approved (FDA 2020, EMA), inhibitor (LoF); dyslipidaemia, atherosclerosis, ASCVD, familial hypercholesterolaemia; 12 registered trials plus label records |
| GNRHR | 817 | Agonists (GoF): leuprolide (149 trials), goserelin, triptorelin, buserelin, histrelin, nafarelin, gonadorelin. Antagonists (LoF): degarelix (80), relugolix (37), elagolix, linzagolix, cetrorelix, ganirelix, abarelix |
| 11 marker genes, DIDO1, ARL9, LZTR1 | 0 | None |

**Partner genes (all inhibitors, direction LoF):**

| Partner | Clinical drugs in KG | Highest stage | Reading for the proband |
|---|---|---|---|
| AURKB | 18 agents incl. barasertib (AZD1152), alisertib, danusertib, chiauranib, AT-9283 | Phase 3 | Contraindication: weakens checkpoint signalling |
| PLK1 | Volasertib, onvansertib, BI 2536, rigosertib and others | Phase 3 | Contraindication: anti-mitotic |
| TTK / MPS1 | BAY1161909, BAY1217389 | Phase 1 | Contraindication: overrides the checkpoint |
| PLK4 | CFI-400945, ocifisertib, RG1530 | Phase 2 | Not relevant (TRIM37-amplified tumours) |
| MTOR | Sirolimus (approved, boxed warning), everolimus (approved; boxed warning only on the transplant label Zortress, not on Afinitor) plus 20+ clinical inhibitors | Approval | Already in sheet as rank 4 (acts on a consequence) |
| NAMPT | Daporinad, teglarinad, GMX1777 | Phase 2 | Contraindication: lowers NAD+, opposite of the niacin rationale |
| SIRT2 | None clinical (103 assay molecules) | None | No activator available |
| BUB3, MAD2L1, SMC6, NSMCE2 | None clinical | None | No drug |

### 4.3 Synthetic lethality and therapeutic-target evidence

| Gene | Synthetic lethal records | Therapeutic-target records | Main contexts and drug angle |
|---|---|---|---|
| BUB1B | 11 | 33 | Dependency in TP53-altered paediatric AML (PMID 37002311) and MAPK-active glioblastoma (36191188, 28855212). **New:** BUBR1-depleted MCF-7 cells more sensitive to cisplatin plus anti-microtubule agent (30468712) |
| BUB1 | 4+ | 40+ | ATR or PARP inhibitor synergy in breast cancer (38396175); mesothelioma, sarcoma, lipedema |
| CDC20 | 6 | 124 | High CDC20 in aneuploid AML; apcin, proTAME, compound E1 (preclinical) |
| SLF2 | 4 | 1 | SUMOylation inhibitor TAK-981 in SLF2-deficient lymphoma (37485814) |
| SMC5 | 4 | 1 | APOBEC3A activity, SETX deficiency (38886582, 41533569); no drug |
| TRIM37 | 1 | 8 | PLK4 inhibition in TRIM37-amplified tumours (35864973) |
| TRIP13 | 1 | 33 | Selinexor suppresses TRIP13 (Wilms tumour, 38589567) |
| ACLY | 7 | 100+ | GOT plus ACLY inhibition in KRAS cancers (36269753); lipid lowering |
| LZTR1 | 3 | 5 | Two-hit loss in glioblastoma (23917401); KRAS-LZTR1 interface (39194017) |
| MAD1L1 | 0 | 1 | Expression in cancer (38143740) |
| CEP57, CENATAC, MAD2L1BP, DIDO1, ARL9, GNRHR | 0 | 0 to 99 | GNRHR therapeutic records are hormone indications |

These records describe cancers where the gene is **over-expressed or needed**, not germline loss. For a germline loss-of-function patient they point to drugs to avoid or to cancer-only strategies.

### 4.4 Allele similarity (drug classes that target the same type of allele)

| Proband allele | Allele type | Drug class that targets this allele type | Fit |
|---|---|---|---|
| BUB1B p.Leu737Ter | Stop-gain, UGA | Aminoglycoside read-through (gentamicin); ataluren no longer authorised in the EU | Best allele-level match; UGA is the most read-through-prone stop. Limit: truncating MVA alleles leave little transcript |
| BUB1B p.Asn1002Lys | Kinase-domain missense, unstable protein | Protein stabilisation via NAD+/SIRT2 (niacin); no chaperone in KG | Partial: targets low BUBR1 abundance |
| LZTR1 p.Tyr748Ter | Stop-gain, UAG (c.2244C>G turns TAC into TAG) | Read-through agents | Possible, less efficient than UGA; secondary finding, low priority |
| GNRHR p.Arg139His | Missense in the conserved DRS motif | Pharmacoperones (small-molecule GnRHR chaperones) | Class match; R139H may impair signalling rather than folding, so rescue needs confirmation; carrier finding |
| ACLY p.Arg244Trp | Missense (predicted damaging) | None (no activator or chaperone) | No match; the approved drug inhibits ACLY |
| DIDO1 5′ UTR splice deletion | Splice | Bespoke antisense oligonucleotide only | No off-the-shelf drug |
| ARL9 p.Thr240Argfs*16 | Frameshift | None | Not targetable |
| SMC5 upstream deletion | Non-coding structural | None | Not targetable |

## 5. Changes against the team sheet

| Item | Team sheet | This KG pass |
|---|---|---|
| ACLY drug evidence | "absent from Polly KG v2 drug edges" | Present: 22 Open Targets clinical records for bempedoic acid (approved; 12 registered trials plus FDA/EMA label records); verdict unchanged |
| Contraindications | Only BUB1/CDC20/TRIP13 inhibitors marked wrong direction | Added Aurora B, PLK1, MPS1/TTK and NAMPT inhibitors as explicit contraindications |
| BUB1B cancer context | Dependency rows only | Added cisplatin plus anti-microtubule combination in BUBR1-low cells (cancer only) |
| Allele-type rows | BUB1B UGA only | Added LZTR1 UAG read-through and GNRHR pharmacoperone class |
| SIRT2 / NAMPT | Mentioned | Confirmed: no clinical SIRT2 activator; NAMPT drugs all inhibitors |
| Pharmacogenomics | Not checked | KG holds no variant-to-drug-response records for the 16 genes |

## 6. Shortlist for the deep dive

| Rank | Idea | Why it stays | Main open question |
|---|---|---|---|
| 1 | Niacin / NAD+ precursors | Only approved route that aims to raise BUBR1 (mouse evidence with NMN and SIRT2) | Any effect on human BUBR1 levels; paediatric dosing |
| 2 | Gentamicin read-through (UGA allele) | Allele-specific, approved, paediatric pilot in another nonsense disorder | Is enough BUB1B transcript left from the stop-gain allele? |
| 3 | mTORC1 inhibition (sirolimus/everolimus) | Acts on a measured consequence in BubR1 mice | Boxed warning on sirolimus and on the everolimus transplant label (Zortress, not Afinitor); benefit unproven |
| 4 | Senolytics (dasatinib with quercetin) | Senescent-cell clearance helped BubR1 progeroid mice (genetic, not drug) | No drug study in BubR1 models |
| Avoid | Aurora B, PLK1, MPS1, BUB1, NAMPT inhibitors; anti-mitotic chemotherapy | Further weaken the checkpoint or lower NAD+ | List as contraindications in the report |

All entries are research hypotheses for functional testing, not clinical advice.

## 7. Query log (read-only Cypher)

| # | Purpose | Pattern |
|---|---|---|
| Q1 | Schema of compound edges | `MATCH (g:gene {id:'699'})-[:encodes]->(p)<-[r]-(m) RETURN type(r), labels(m), keys(m), keys(r)` |
| Q2 | Molecule properties | keys of molecules on MTOR (`binds_target`, `inhibits_target`) |
| Q3 | Molecule relationship types | relationships of approved MTOR molecules |
| Q4 | Direct compounds, 16 genes | `encodes`→protein←molecule; counts, `max_phase`, `first_approval`, `black_box_warning`, `action_type`, `pchembl_value` |
| Q5 | Partner compounds | same for AURKB, PLK1, TTK, PLK4, MTOR, NAMPT, SIRT2, BUB3, MAD2L1, SMC6, NSMCE2 |
| Q6 | Association types on the 16 genes | `(g)-[:gene_participates_in]->(a) RETURN labels(a)[0], count(*)` |
| Q7 | Synthetic lethal and therapeutic-target records | association nodes with `evidence_text`, `pmids`, linked disease |
| Q8 | BUB1B synthetic lethal detail | all 11 records with evidence text |
| Q9 | Drug resistance / predictive biomarker / has_drug_target | per gene, linked diseases and PMIDs |
| Q10 | Open Targets clinical drugs | `has_drug_target` properties: `drug_id`, `drug_from_source`, `clinical_stage`, `direction_on_target`, `disease_from_source`, `clinical_report_id`, for the 16 genes and partners |

Gene IDs (Entrez): BUB1 699, BUB1B 701, CDC20 991, TRIP13 9319, MAD2L1BP 9587, MAD1L1 8379, SMC5 23137, SLF2 55719, CEP57 9702, CENATAC 338657, TRIM37 4591, ACLY 47, DIDO1 11083, ARL9 132946, LZTR1 8216, GNRHR 2798; partners MAD2L1 4085, BUB3 9184, TTK 7272, AURKB 9212, PLK1 5347, SIRT2 22933, NAMPT 10135, SMC6 79677, NSMCE2 286053, PLK4 10733, MTOR 2475.

## 8. Limits

- Open Targets clinical-precedence edges cover drugs with a declared target; off-target effects (for example sunitinib on TTK) appear only as ChEMBL binding edges.
- Synthetic-lethal and therapeutic-target records are text-mined claims; each was read as evidence text, not as a validated finding.
- The pharmacoperone and read-through notes for GNRHR and LZTR1 come from allele-type reasoning and need literature confirmation.
- No query used patient coordinates.
