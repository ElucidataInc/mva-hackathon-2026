# Open Targets check of Top10_helper variants with evidence count >= 3

Sample WGS_EX2312012 (child with suspected mosaic variegated aneuploidy, MVA). Source of the variant list: the Top10_helper sheet. Results come from `opentargets_check.tsv`, produced by `check_opentargets.py` (second, fixed run). Research use only.

## 1. What we searched

- **Selection:** every row of Top10_helper with an evidence count of 3 or more. That is 18 rows: evidence 6 (2 rows), 5 (1), 4 (1) and 3 (14).
- **Database and access:** Open Targets Platform GraphQL API (`api.platform.opentargets.org/api/v4`). The Open Targets release number was not recorded by the script.
- **Variant-level lookup (14 rows):** rows whose label has full alleles were queried as GRCh38 variant IDs in the form `chr_pos_ref_alt`. From each hit we read rsID, most severe consequence, population allele frequencies and variant-effect predictions (VEP, LOFTEE, SIFT, GERP, FoldX, AlphaMissense).
- **Gene-level lookup (all rows):** each gene symbol was mapped to its Ensembl ID, then the top 300 disease associations by score were pulled. The script flagged associations whose names contain MVA, aneuploidy, chromosome, microcephaly, premature chromatid separation, cancer predisposition, Wilms, mitotic or segregation.
- **Rows with no variant-level lookup (4):** the BUB1B/CHST14 inversion, the BUB1B duplication and the SMC5 deletion are structural variants, and Open Targets has no structural-variant pages. The ADAMTS1 p.Thr220ArgfsTer76 allele is truncated in the sheet (`GTCCTCGC..(15)>G`). These rows got the gene-level check only.
- **First run:** the variant queries failed because of a field-name error in my script, so the "not found" results from that run were invalid and are not used. The results below are from the fixed run.

## 2. How to read the results

- **Found / not found:** "not found" means Open Targets returned no record for that variant ID. It does not mean benign or rare. Open Targets holds only a subset of variants, and very rare ones (AF about 1e-6) are often absent. Indels may also be written differently in its index.
- **AFs below** are Open Targets' population values; `nfe` is non-Finnish European, `remaining` is the gnomAD remaining group.
- **Gene scores** (0 to 1) are gene-disease association scores. They say whether the gene is linked to a disease, not whether this variant is pathogenic.
- **Low scores:** association hits below 0.1 for terms like "chromosome ... deletion syndrome" are keyword noise and are ignored.

## 3. Summary

- **Confirmed:** BUB1B p.Leu737Ter is a rare, high-confidence loss-of-function stop in a gene strongly linked to MVA (scores 0.81 and 0.79).
- **Common, so excluded:** FANCD2 splice donor (AF 0.35) and MAD1L1 p.Arg59Cys (AF 0.0057).
- **Gene supports the variant class, variant not found:** SLF2 and SMC5 are linked to Atelis syndrome (0.71 and 0.69), which fits the MVA-like picture.
- **No MVA link at all:** ACLY, ARL9, CACNB1, SERPINA1, GNRHR, LZTR1, IRX4 and ADAMTS1 (a weak microcephaly score of 0.03 only).
- **Sheet flag problem:** the sheet's "rare" flag was set for FANCD2 because no AF was reported, but Open Targets shows 35%. Blank AF means unknown, not rare.

## 4. Findings per variant

### Row 2. BUB1B p.Leu737Ter (chr15:40209701 T>G), evidence 6

This is the strongest result. Open Targets confirms a rare, high-confidence stop-gained in BUB1B, and BUB1B has the clearest MVA association of any gene in the list.

- Variant ID searched: `15_40209701_T_G`, **found**, rs759242053.
- Consequence: stop_gained; VEP 0.90; LOFTEE high-confidence LoF; GERP 0.969.
- AF: nfe 9.83e-5, remaining 6.4e-5; 0 in sas, mid, fin, eas, asj and amr. This matches the sheet's 9.98e-5.
- Gene BUB1B (ENSG00000156970), 589 associations: mosaic variegated aneuploidy syndrome 0.81, MVA syndrome 1 0.79, premature chromatid separation trait 0.63, microcephaly 0.13.
- Other top hits: neurodegenerative disease 0.54, colorectal cancer 0.50.

### Row 3. BUB1B p.Asn1002Lys (chr15:40220612 T>G), evidence 6

The variant is not in Open Targets, so there is no new variant-level information. Its AlphaMissense score of 0.92 comes from the sheet, not from Open Targets.

- Variant ID searched: `15_40220612_T_G`, **not found**. This is consistent with the sheet's AF of about 9e-7.
- Gene-level results are the same as row 2 (BUB1B, MVA 0.81 and 0.79).

### Row 4. BUB1B, CHST14 inversion (chr15:38852774), evidence 5

Only BUB1B has an MVA link. CHST14 and EIF2AK4 have none, so the inversion gets no gene-level support beyond lying 735 kb from BUB1B.

- No variant-level lookup (structural variant).
- BUB1B: as row 2.
- CHST14 (ENSG00000169105), 428 associations: top hits are Ehlers-Danlos syndrome, musculocontractural type 0.84, cardiovascular abnormality 0.47 and Ehlers-Danlos syndrome 0.42. The only microcephaly-type terms are about 0.03 to 0.04.
- EIF2AK4 (ENSG00000128829), 288 associations: pulmonary veno-occlusive disease 2 0.77 and pulmonary arterial hypertension 0.55 to 0.56. No MVA-related terms.

### Row 5. BUB1B duplication (chr15:40202405), evidence 4

The duplication sits inside the BUB1B gene span, but Open Targets cannot assess structural variants.

- No variant-level lookup.
- Gene BUB1B: as row 2.

### Row 6. MAD1L1 p.Arg59Cys (chr7:2225526 G>A), evidence 3

The gene is a real MVA gene, but this variant is common and predicted tolerated, so it does not support an MVA explanation.

- Variant ID searched: `7_2225526_G_A`, **found**, rs121908982, missense_variant.
- AF: nfe 5.68e-3, remaining 3.18e-3, fin 2.67e-3, amr 8.16e-4, mid 1.65e-4.
- Predictions: VEP 0.68; SIFT tolerated (0.06); GERP 1.96; FoldX 1.35.
- Gene MAD1L1 (ENSG00000002822), 226 associations: MVA syndrome 7 with inflammation and tumour predisposition 0.53, MVA syndrome 1 0.43.
- Other top hits: schizophrenia 0.50, neurodegenerative disease 0.45.

### Row 7. CACNB1 p.Ter599TrpextTer12 (chr17:39175193 T>C), evidence 3

A stop-loss in a gene with no MVA link. The variant is not in Open Targets.

- Variant ID searched: `17_39175193_T_C`, **not found**.
- Gene CACNB1 (ENSG00000067191), 374 associations: epilepsy 0.61, neuropathic pain 0.60, seizure 0.60. The only chromosome-type term is maternal uniparental disomy of chromosome 21 at 0.04.

### Row 8. ADAMTS1 p.Glu225LeufsTer13 (chr21:26844275 GCCCCTTC>G), evidence 3

Not found, and the gene has no real MVA link. The sheet also flags this variant and row 9 as likely on the same copy (cis).

- Variant ID searched: `21_26844275_GCCCCTTC_G`, **not found** (indel representation may differ).
- Gene ADAMTS1 (ENSG00000154734), 494 associations: microcephaly 20, primary, autosomal recessive 0.03; Wilms tumour 0.02. Both are weak.
- Top hits: autoimmune thyroid disease 0.43, neurodegenerative disease 0.33.

### Row 9. ADAMTS1 p.Thr220ArgfsTer76 (chr21:26844283), evidence 3

No variant-level lookup, because the allele is truncated in the sheet. Gene-level results are the same as row 8.

### Row 10. ARL9 p.Thr240ArgfsTer16 (chr4:56523792 CCTGA>C), evidence 3

- Variant ID searched: `4_56523792_CCTGA_C`, **not found**.
- Gene ARL9 (ENSG00000196503), 25 associations, none MVA-related. Top hit: colon adenocarcinoma 0.07.

### Row 11. SLF2 intron variant (chr10:100939590 A>G), evidence 3

The variant is not in Open Targets and is intronic. The gene is the notable part: SLF2 is linked to Atelis syndrome.

- Variant ID searched: `10_100939590_A_G`, **not found**.
- Gene SLF2 (ENSG00000119906), 57 associations: Atelis syndrome 1 0.71, abnormality of the skeletal system 0.47. Microcephaly 0.03.

### Row 12. FANCD2 splice donor variant (chr3:10046723 AG>A), evidence 3

Open Targets shows this variant is common, so it is not a real loss-of-function candidate. The sheet's "rare" flag came from a missing AF.

- Variant ID searched: `3_10046723_AG_A`, **found**, rs750338758, splice_donor_variant.
- AF: **nfe 0.352**, fin 0.390, mid 0.347, remaining 0.320, asj 0.291, eas 0.244, amr 0.226, sas 0.141.
- Predictions: VEP 0.93; LOFTEE high-confidence (despite the high AF); GERP 1.83.
- Gene FANCD2 (ENSG00000144554), 595 associations: Fanconi anaemia complementation group D2 0.83, Fanconi anaemia 0.81, acute myeloid leukaemia 0.64. No MVA terms.

### Row 13. IRX4 p.Arg165Cys (chr5:1879747 G>A), evidence 3

- Variant ID searched: `5_1879747_G_A`, **not found**.
- Gene IRX4 (ENSG00000113430), 480 associations: prostate carcinoma 0.57, prostate cancer 0.53, ventricular septal defect 1 0.43. Microcephaly-cardiomyopathy syndrome 0.06 and microcephaly with increased sister chromatid exchange 2 0.04 are weak.

### Row 14. SMC5 deletion (chr9:70286200), evidence 3

Gene support is good, but Open Targets cannot assess this depth-based deletion.

- No variant-level lookup.
- Gene SMC5 (ENSG00000198887), 311 associations: Atelis syndrome 2 0.69, abnormality of the skeletal system 0.30. Microcephaly 0.04; hereditary thrombocytopenia and haematological cancer predisposition syndrome 0.03.

### Row 15. ACLY p.Arg244Trp (chr17:41907459 G>A), evidence 3

- Variant ID searched: `17_41907459_G_A`, **not found**.
- Gene ACLY (ENSG00000131473), 650 associations, none MVA-related. Top hits: familial hypercholesterolaemia 0.53 and atherosclerosis 0.47.

### Row 16. GNRHR p.Arg139His (chr4:67753920 C>T), evidence 3

This is a known hypogonadism variant and has no link to MVA.

- Variant ID searched: `4_67753920_C_T`, **found**, rs104893842, missense_variant.
- AF: nfe 1.55e-4, fin 5.16e-4, amr 1.83e-4, remaining 4.80e-5, eas 2.23e-5, sas 1.10e-5, mid and asj 0.
- Predictions: VEP 0.68; SIFT deleterious (0); GERP 2.47; AlphaMissense ambiguous (0.506); FoldX 4.46.
- Gene GNRHR (ENSG00000109163), 646 associations: hypogonadotropic hypogonadism 0.83, prostate cancer 0.61.
- Co-annotated gene UBA6-DT (ENSG00000248049): 0 associations.

### Row 17. LZTR1 p.Tyr748Ter (chr22:20996720 C>G), evidence 3

An ultra-rare stop in a real disease gene, but the disease is Noonan syndrome and schwannomatosis, not MVA.

- Variant ID searched: `22_20996720_C_G`, **found**, rs1682503990, stop_gained.
- AF: nfe 8.47e-7; 0 in every other group.
- Predictions: VEP 0.90; LOFTEE high-confidence; GERP -0.579 (low conservation).
- Gene LZTR1 (ENSG00000099949), 1,222 associations: Noonan syndrome 0.84, schwannomatosis 0.82, Noonan syndrome 10 0.79. No MVA terms.
- Co-annotated genes: THAP7 (42 associations, none relevant; top hit pain 0.10) and THAP7-AS1 (4 associations, none relevant).

### Rows 18 and 19. SERPINA1 p.(Glu347Glyfs*12) and p.(Thr346Serfs*7) (chr14:94379487 CCT>C and chr14:94379492 G>GAC), evidence 3

- Variant IDs searched: `14_94379487_CCT_C` and `14_94379492_G_GAC`, both **not found**.
- Gene SERPINA1 (ENSG00000197249), 1,408 associations, none MVA-related. Top hits: alpha-1 antitrypsin deficiency 0.77, chronic obstructive pulmonary disease 0.69, liver cirrhosis 0.57.

## 5. Limits

- Open Targets gives gene-disease links and population data. It does not classify these variants as pathogenic for this child.
- Variants "not found" cannot be ranked up or down on this evidence alone.
- Structural variants (rows 4, 5 and 14) were checked at gene level only.
- ClinVar status and AlphaMissense for variants Open Targets did not return come from the sheet, not from this search.
