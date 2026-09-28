# 04. BubR1 biology and intervention points

What BubR1 does, what the mouse models show, what goes wrong in MVA1 cells, and which nodes a drug could act on.

Abbreviations: SAC, spindle assembly checkpoint. MCC, mitotic checkpoint complex. APC/C, anaphase-promoting complex. K-MT, kinetochore-microtubule. KARD, kinetochore attachment regulatory domain. CIN, chromosomal instability. MEF, mouse embryonic fibroblast.

## Summary

- BubR1 has two mitotic jobs. It is an MCC subunit that inhibits APC/C-Cdc20, and its phosphorylated KARD recruits PP2A-B56 to stabilise K-MT attachments. Its kinase domain does no essential catalysis; it holds the protein's fold (refs 1, 3, 12).
- MVA1 is a dosage defect. A truncating allele makes no transcript and a missense allele makes unstable protein, so total BubR1 is low. Patient cells have a weak checkpoint and misaligned chromosomes, and adding back BubR1 rescues them (ref 4).
- Raising BubR1 has the best mouse evidence: 10-fold overexpression extends lifespan, cuts tumours and prevents age-related aneuploidy (ref 5).
- SIRT2 stabilises BubR1. SIRT2 overexpression extended median lifespan of BubR1 hypomorphs by 58%, but only in males. NMN raised BubR1 protein in mouse testis after 7 days (ref 6).
- Clearing senescent cells in BubR1 hypomorphs delayed sarcopenia, cataract and fat loss but did not extend survival (ref 7).
- *BUB1B* has no approved or clinical-stage drug in any modality (Open Targets). Any drug has to act indirectly.
- Drugs that weaken the SAC or destabilise K-MT attachments should make MVA1 worse: MPS1/TTK, TRIP13 or PLK1 inhibitors, and anything that releases APC/C-Cdc20.

## 1. Molecular function

### 1.1 The checkpoint

The SAC holds anaphase until every chromosome is attached to both spindle poles. Unattached kinetochores catalyse the MCC, which binds and inhibits APC/C, the ubiquitin ligase that triggers segregation (refs 9, 12).

The fission-yeast MCC structure (Mad2, Mad3 and Cdc20; ref 12) shows that:
- MCC blocks the degron-recognition sites on Cdc20 and displaces it, so the D-box receptor with Apc10 cannot form.
- Closed Mad2 positions the Mad3 KEN box to bind Cdc20.
- Mad3 and p31comet compete for the same Mad2 surface, which is how p31comet switches the checkpoint off.

Mad3 is the yeast orthologue of BubR1. The KEN-box logic applies to human BUBR1 (ref 3), but the human MCC structure is inferred rather than solved.

### 1.2 Attachment through KARD and PP2A-B56

Suijkerbuijk 2012 (ref 1) showed that:
- Deleting KARD, or blocking its phosphorylation, stops kinetochore microtubules from forming. Inhibiting Aurora B reverses this, so the cause is excess Aurora B phosphorylation at the outer kinetochore.
- PLK1 phosphorylates KARD, which lets BUBR1 bind PP2A-B56alpha.
- Removing BUBR1 or inhibiting PLK1 lowers PP2A-B56alpha at kinetochores.

Low BubR1 therefore costs MVA1 cells both the checkpoint and the phosphatase at once, which explains why they show a weak checkpoint and misaligned chromosomes together (refs 1, 4).

### 1.3 The kinase domain

Human BUBR1 is a pseudokinase (ref 3). Its catalytic activity is not needed for accurate segregation, but the residues that would contact ATP are needed to keep the protein folded. MVA1 missense mutations that speed up turnover cluster in or near this domain (ref 4). Two consequences for drug design:
- A BubR1 "kinase activator" makes no sense.
- An ATP-site ligand would more likely destabilise the protein than help it.

## 2. Mouse models

| Study | Model | Result | Limits |
|---|---|---|---|
| Baker 2004 (ref 13) | BubR1 hypomorph (H/H) | Progressive aneuploidy, short life, dwarfism, lordokyphosis, cataract, fat loss, poor wound healing, infertility. Wild-type mice lose BubR1 in several tissues with age | A progeroid phenotype, without microcephaly or childhood cancer |
| Baker 2013 (ref 5) | BubR1 overexpression (line T23, ~10x in MEFs) | Longer median and maximum lifespan (max p = 0.05 and 0.0056). DMBA tumours in 33% vs 100% of wild type. Kras lung tumours 8 vs 17. No age-related aneuploidy in lung, muscle, kidney, eye and heart | No effect in bone marrow or gut; osteoporosis and cataract not delayed; median gain not stated |
| Baker 2011 (ref 7) | Senescent-cell clearance (INK-ATTAC + AP20187) in H/H mice | Lifelong clearance delayed lordokyphosis and cataract, preserved muscle, exercise capacity and fat. Clearance from 5 months slowed progression | Survival "not substantially extended". p16-independent cardiac and arterial defects untouched. AP20187 is a tool compound |
| Zhu 2015 (ref 14) | Dasatinib + quercetin (D+Q) in aged, irradiated and *Ercc1*(-/delta) mice | Fewer senescent cells, better cardiac and vascular function, longer healthspan in *Ercc1* mice | The progeroid model is *Ercc1*, a DNA-repair model, not BubR1. Dasatinib and quercetin target different senescent cell types |
| North 2014 (ref 6) | SIRT2 overexpression in H/H; NMN in wild type | SIRT2 keeps BubR1 K668 deacetylated against CBP. Median lifespan +58% overall (n = 31 vs 33, p = 0.038): +123% in males (p = 0.004), no change in females (p = 0.79). NMN 500 mg/kg/day i.p. for 7 days restored testis BubR1 in old mice | Effect is male-only. The NMN data cover one tissue, one week and protein level only, in wild-type mice |

A full BubR1 knockout is embryonic lethal between E7.5 and E13.5, with over 50% aneuploid cells and massive apoptosis (ref 15). MVA1 is a hypomorphic state. Age-related BubR1 decline has been shown only in mice.

## 3. MVA1 in human cells

### 3.1 Patient cells

In MVA patient cell lines (Suijkerbuijk 2010, ref 4), BubR1 is low, the checkpoint is weak and chromosomes misalign. Ectopic BubR1 restores checkpoint activity, so the BubR1 defect causes the errors. Some mutations disable a specific function; others lower abundance, through absent transcript (truncating alleles) or fast turnover (missense alleles near the kinase domain). The most tractable axis is BubR1 abundance and stability.

### 3.2 Aneuploidy and PCS

- PCS accompanies MVA in *BUB1B* patients (refs 15, 16, 17).
- In García-Castillo 2008 (21 of 35 MVA patients screened for *BUB1B*), the monoallelic group was the most severe: Dandy-Walker 7/8, cataract 6/6, Wilms 7/8, PCS in 8/8 probands and 7/7 carrier parents. The biallelic group was moderate. Patients without a *BUB1B* mutation were mild, without cancer, and rarely showed PCS (1/7). Many "monoallelic" patients may carry an undetected noncoding second allele (Ochiai 2014), which would make this consistent with a dosage model.
- Two *BUB1B*-negative siblings had about half of their lymphocytes aneuploid: trisomies 42%, monosomies 28% (ref 17).
- One MVA1 case showed 12% mosaicism on karyotype (ref 16). The aneuploid fraction varies widely between patients and tissues.

### 3.3 Microcephaly

No patient brain data exist for *BUB1B*. The model below comes from other SAC genes and aneuploidy studies.
- Deleting *Knl1* in embryonic mouse brain causes missegregation, DNA damage and p53-driven apoptosis and microglial clearance of neural progenitors. Removing p53 rescued head size but killed every pup before weaning (ref 19).
- In *Drosophila*, aneuploid neural stem cells lengthen G1, exit the cycle and differentiate early, shrinking the brain without apoptosis (ref 20).

Both routes deplete progenitors. Blocking apoptosis does not fix either, and ref 19 shows it can be lethal.

### 3.4 Cancer

Biallelic *BUB1B* and *TRIP13* patient cells have severe SAC impairment and embryonal tumours; *CEP57* and unexplained MVA cells have minimal SAC impairment and no embryonal tumours (ref 21). Checkpoint strength, not aneuploidy alone, tracks tumour risk. *BUB1B* MVA tumours include rhabdomyosarcoma and Wilms (refs 17, 22). A homozygous intronic c.2386-11A>G carrier developed ampullary adenocarcinoma at 34, then colon and gastric adenocarcinomas (ref 2).

## 4. Vulnerabilities of aneuploid cells

| Finding | System | Ref |
|---|---|---|
| Aneuploid strains share cell-cycle defects, higher glucose uptake and sensitivity to protein-synthesis and folding stress, caused by the extra genes | Yeast disomes | 8 |
| AICAR, 17-AAG (HSP90 inhibitor) and chloroquine selectively slow aneuploid cells; AICAR triggers p53 apoptosis in trisomic MEFs | Trisomic MEFs, cancer lines | 10 |
| Aneuploid cells overload lysosomes, activating TFEB | Human cells | 11 |
| Across ~1,000 cancer lines, aneuploid cells are more sensitive to genetic SAC loss but less sensitive to short-term SAC inhibitors, and depend on KIF18A | Cancer lines | 9 |
| KIF18A inhibitor ATX020 arrests CIN-high ovarian cancer cells | Ovarian cancer lines | 23 |

Of these compounds only chloroquine is approved. AICAR, 17-AAG and KIF18A inhibitors are not.

All of these studies use cancer or engineered cells that tolerate aneuploidy. MVA1 tissue is aneuploid but not transformed, and killing aneuploid cells could deepen the progenitor loss behind microcephaly (section 3.3). Selective killing makes sense only for tumour prevention or treatment.

## 5. Intervention points

The MVA1 defect is low BubR1, a weak SAC and unstable K-MT attachment.

| Node | Direction | Evidence | PMIDs | Strength |
|---|---|---|---|---|
| BubR1 abundance | Increase | Rescue in patient cells; overexpression protects mice | 20516114; 23242215 | Strong, but no drug does it |
| Missense BubR1 turnover | Decrease | Missense protein is degraded fast and is the only protein made | 20516114 | Strong mechanism, no drug |
| SIRT2 / BubR1 K668 | Increase SIRT2 | Stabilises BubR1; extends lifespan in male hypomorphs | 24825348 | Moderate |
| NAD+ (NMN, NR) | Increase | Raises BubR1 in mouse testis | 24825348 | Weak |
| CBP at K668 | Decrease | Opposes SIRT2 | 24825348 | Weak; CBP inhibition is broad |
| PP2A-B56 at kinetochores | Preserve | Needed for stable attachment | 23079597 | Moderate; no selective agent |
| Aurora B | Decrease modestly | Reverses the KARD attachment defect | 23079597 | High risk; Aurora B also drives error correction |
| p16+ senescent cells | Decrease | Healthspan gains in hypomorphs | 22048312 | Moderate; no survival gain |
| Senolytics (D+Q) | Decrease | Healthspan in *Ercc1* mice | 25754370 | Weak; not a BubR1 model |
| HSP90 (17-AAG) | Tumour only | Aneuploidy-selective | 21315436 | Harmful to progenitors |
| AMPK (AICAR) | Tumour only | p53 apoptosis in trisomic MEFs | 21315436 | AICAR not approved |
| Autophagy-lysosome (chloroquine) | Unclear | Aneuploidy-selective | 26404941; 21315436 | Approved, but benefit in MVA1 unknown |
| KIF18A | Decrease, tumour only | Aneuploid cancer cells depend on it | 33505028; 41369352 | No approved agent |
| p53 clearance of damaged progenitors | Do not suppress | Rescued head size, killed pups | 31197172 | Strong negative result |

### 5.1 Interventions expected to cause harm

| Node | Intervention | Why | PMIDs |
|---|---|---|---|
| MPS1/TTK | Inhibit | Weakens a checkpoint that is already weak | 33505028; 20516114 |
| TRIP13 | Inhibit | *TRIP13* loss causes MVA3 with Wilms tumour | 28553959 |
| PLK1 | Inhibit | Removes PP2A-B56alpha from kinetochores, the axis already broken | 23079597 |
| p31comet | Activate | Competes with BUBR1 for Mad2 and dismantles the MCC | 22437499 |
| CDC20 | Release from APC/C inhibition | Equivalent to switching the checkpoint off | 22437499 |
| *BUB1B* ATP pocket | Bind | Likely destabilises an already scarce protein | 22698286; 20516114 |
| Aurora B | Strong inhibition | Blocks error correction; dose-dependent | 23079597 |
| p53 in brain | Suppress | Lethal in *Knl1* mice | 31197172 |

### 5.2 *BUB1B* tractability

Open Targets (ENSG00000156970):
- Small molecule: druggable family and high-quality ligand flags only; no approved, clinical or ligand-bound structure, and no high- or medium-quality pocket.
- Antibody: none (intracellular target).
- PROTAC: literature and ubiquitination evidence only.

MVA1 needs BubR1 stabilised, not degraded. A repurposed approved drug would have to act through abundance regulators (NAD+/SIRT2) or downstream effects (senescent cells, tumour surveillance).

## Tensions

- Cohen-Sharir 2021: aneuploid cells are more sensitive to genetic SAC loss but less sensitive to short-term SAC inhibitors. "SAC inhibitors kill aneuploid cells" needs that caveat.
- Baker 2011 improved healthspan without survival gain, because the cardiac cause of death is p16-independent.
- Killing aneuploid cells (oncology) and preserving progenitors (neurodevelopment) pull in opposite directions. A systemic intervention that does one will tend to undo the other.

## Open questions

- Whether BubR1 declines with age in human tissue, as it does in mice.
- Whether any study of brain tissue from *BUB1B* patients exists.

## References

1. Suijkerbuijk SJE, et al. Integration of kinase and phosphatase activities by BUBR1 ensures formation of stable kinetochore-microtubule attachments. Dev Cell. 2012. PMID 23079597.
2. Rio Frio T, et al. Homozygous *BUB1B* mutation and susceptibility to gastrointestinal neoplasia. N Engl J Med. 2010. PMID 21190457.
3. Suijkerbuijk SJE, et al. The vertebrate mitotic checkpoint protein BUBR1 is an unusual pseudokinase. Dev Cell. 2012. PMID 22698286.
4. Suijkerbuijk SJE, et al. Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy. Cancer Res. 2010. PMID 20516114; PMC2887387.
5. Baker DJ, et al. Increased expression of BubR1 protects against aneuploidy and cancer and extends healthy lifespan. Nat Cell Biol. 2013. PMID 23242215; PMC3707109.
6. North BJ, et al. SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. EMBO J. 2014. PMID 24825348; PMC4194088.
7. Baker DJ, et al. Clearance of p16Ink4a-positive senescent cells delays ageing-associated disorders. Nature. 2011. PMID 22048312; PMC3468323.
8. Torres EM, et al. Effects of aneuploidy on cellular physiology and cell division in haploid yeast. Science. 2007. PMID 17702937.
9. Cohen-Sharir Y, et al. Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition. Nature. 2021. PMID 33505028.
10. Tang YC, et al. Identification of aneuploidy-selective antiproliferation compounds. Cell. 2011. PMID 21315436.
11. Santaguida S, et al. Aneuploidy-induced cellular stresses limit autophagic degradation. Genes Dev. 2015. PMID 26404941.
12. Chao WCH, et al. Structure of the mitotic checkpoint complex. Nature. 2012. PMID 22437499.
13. Baker DJ, et al. BubR1 insufficiency causes early onset of aging-associated phenotypes and infertility in mice. Nat Genet. 2004. PMID 15208629.
14. Zhu Y, et al. The Achilles' heel of senescent cells: from transcriptome to senolytic drugs. Aging Cell. 2015. PMID 25754370.
15. Schmid M, et al. Mosaic variegated aneuploidy in mouse BubR1 deficient embryos and pregnancy loss in human. Chromosome Res. 2014. PMID 24981203.
16. Pavone P, et al. Pathogenic correlation between mosaic variegated aneuploidy 1 (MVA1) and a novel *BUB1B* variant. Neurol Sci. 2022. PMID 35804254.
17. García-Castillo H, et al. Clinical and genetic heterogeneity in patients with mosaic variegated aneuploidy: delineation of clinical subtypes. Am J Med Genet A. 2008. PMID 18548531.
18. Open Targets Platform GraphQL API, target ENSG00000156970. https://api.platform.opentargets.org/api/v4/graphql. 
19. Shi L, et al. Robust elimination of genome-damaged cells safeguards against brain somatic aneuploidy following *Knl1* deletion. Nat Commun. 2019. PMID 31197172.
20. Gogendeau D, et al. Aneuploidy causes premature differentiation of neural and intestinal stem cells. Nat Commun. 2015. PMID 26573328.
21. Yost S, et al. Biallelic *TRIP13* mutations predispose to Wilms tumor and chromosome missegregation. Nat Genet. 2017. PMID 28553959.
22. Hanks S, et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in *BUB1B*. Nat Genet. 2004. PMID 15475955.
23. Nair J, et al. The KIF18A inhibitor ATX020 induces mitotic arrest and DNA damage in chromosomally instable high-grade serous ovarian cancer cells. Cells. 2025. PMID 41369352.
