# Track 2: drug repurposing (team RareELleles)

Drug-repurposing candidates for the proband's Track 1 result, a BUB1B compound-heterozygous pair: p.(Leu737Ter), a stop-gain, and p.(Asn1002Lys), a missense change in the kinase domain. Both allele classes lower BUBR1, so the therapeutic direction is to raise functional BUBR1 or treat a downstream consequence, and never to weaken the spindle assembly checkpoint further.

All candidates are computational hypotheses for testing in patient-derived cells, not treatment recommendations.

## How the candidates were found

1. **Genes:** the Track 1 top-10 rows and the 11 MVA marker genes.
2. **Mechanism:** what each variant does (Ensembl VEP, UniProt domains, Exomiser, Pfam and InterPro, GO and KEGG).
3. **Search:** drugs acting on the gene or its pathway (Polly biomedical knowledge graph v2, ChEMBL, Open Targets, PubMed).
4. **Four gates:** approved (ChEMBL max_phase 4); right direction (raise BUBR1 or treat a consequence, never inhibit the checkpoint); evidence (BubR1-specific data first); safety (boxed warnings, label warnings, paediatric labelling).
5. **Rank:** 1 niacin, 2 gentamicin, 3 dasatinib, 4 sirolimus or everolimus. Each has a first experiment in patient-derived cells, and a drug is dropped if missegregation does not fall.

A network-proximity screen was run alongside as a control. No drug survived multiple-testing correction at either network setting, and the drugs closest to significance were mostly chemotherapy, the class this patient should avoid.

## Folder layout

| Path | Contents |
|---|---|
| `results/candidates/RareELleles_track2_drug_ranking.html` | Ranked candidates, the four gates, excluded drug classes, and the first experiment for each candidate. Every claim carries its PMID or database record. |
| `results/candidates/SAC_genes_drug_candidates.xlsx` | Evidence table for the spindle assembly checkpoint genes and candidate drugs. |
| `results/candidates/Track2_drug_table_16genes.tsv`, `Track2_KG_drug_report.md` | Knowledge-graph drug table for 16 genes and its report. |
| `results/candidates/OpenTargets_findings.md`, `opentargets_check.tsv`, `ACLY_OpenTargets_records_PollyKG.tsv` | Open Targets checks of the submitted variants and genes. |
| `results/candidates/gene_module_network.pdf` | The gene module network: 76 genes in 13 modules. |
| `results/network_screen/` | The network-proximity screen: ranked output at STRING confidence 700 (strict) and 400 (looser), the per-drug comparison of both runs, the drugs below FDR 0.25 in both, and the screen's reports. |
| `scripts/checkopentargets.py` | Open Targets GraphQL checks for the shortlisted variants (Python 3, standard library only). |
| `scripts/annotate_candidates.py` | Annotates shortlisted variants from gnomAD, Ensembl VEP, SpliceAI, Open Targets and Europe PMC. Run with `--help` for options. |
| `methods/RareELleles_methods_description.xlsx` | The hackathon methods description form, Track 1 and Track 2. |

## Network screen settings

From the screen's run parameters: module of 22 seed genes; STRING confidence cutoffs 700 and 400; drugs mapped through ChEMBL and DGIdb (5,460 gene-drug rows, 2,963 drugs); for each drug, the average shortest-path distance from its targets to the nearest module gene, compared with 1,000 degree-matched random target sets (z, empirical p, Benjamini-Hochberg FDR across all drugs). With 1,000 permutations the smallest possible p is 1/1001, so the best reachable FDR is 0.174 at STRING 700.

## Data

Only public drug and target databases and literature were used. No individual-level sequencing data is included in this folder.
