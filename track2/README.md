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
| `scripts/network_screen/` | The network-proximity screen pipeline: `01_network.py` to `04_rank.py`, shared helpers (`common.py`), the seed list (`seeds.csv`) and the manual rescue targets (`rescue_targets.csv`). |
| `scripts/checkopentargets.py` | Open Targets GraphQL checks for the shortlisted variants (Python 3, standard library only). |
| `scripts/annotate_candidates.py` | Annotates shortlisted variants from gnomAD, Ensembl VEP, SpliceAI, Open Targets and Europe PMC. Run with `--help` for options. |
| `methods/RareELleles_methods_description.xlsx` | The hackathon methods description form, Track 1 and Track 2. |

## Network screen

Four steps, run in order from `scripts/network_screen/` (each script's header lists its options):

1. `01_network.py`: builds a STRING network (2 hops) and an OmniPath directed network around the 44 seed genes in `seeds.csv`, and reports how each seed connects to the module.
2. `02_drugs.py`: maps drugs to the seeds, the top network neighbours and the rescue targets through UniProt, ChEMBL and DGIdb (5,460 gene-drug rows, 2,963 drugs).
3. `03_ot_evidence.py`: adds Open Targets target evidence (expression, safety, mouse knockouts), FAERS adverse events and ClinicalTrials.gov trial outcomes.
4. `04_rank.py`: for each drug, the average shortest-path distance from its targets to the nearest gene of the 22-gene module (seed tiers A and B) on the human STRING network, compared with 1,000 degree-matched random target sets (z, empirical p, Benjamini-Hochberg FDR across all drugs), then direction-of-effect and safety tiers.

Recorded settings (the `params` sheet of each ranked workbook): module tiers A and B, STRING cutoff 700 (main run) and 400 (sensitivity run), 1,000 permutations, z cutoff -1.64. The two runs share steps 1 to 3; only step 4 was rerun at cutoff 400. With 1,000 permutations the smallest possible p is 1/1001, so the best reachable FDR is 0.174 at STRING 700.

Needs Python 3 with requests, pandas, numpy, networkx, matplotlib and openpyxl.

```bash
cd scripts/network_screen
python 01_network.py --seeds seeds.csv --out out --score 0.7 --add-nodes 40 --score2 0.9
python 02_drugs.py --nodes out/nodes.csv --seeds seeds.csv --rescue rescue_targets.csv --top-neighbours 25
python 03_ot_evidence.py --drugs out/drug_target_long.csv --out out
python 04_rank.py --links 9606.protein.links.v12.0.txt.gz --info 9606.protein.info.v12.0.txt.gz --string-cutoff 700 --out out
python 04_rank.py --links 9606.protein.links.v12.0.txt.gz --info 9606.protein.info.v12.0.txt.gz --string-cutoff 400 --out out_400
```

Step 4 needs the STRING v12.0 human flat files (`9606.protein.links.v12.0.txt.gz`, `9606.protein.info.v12.0.txt.gz`) from string-db.org. The seed list's free-text notes are left out of `seeds.csv`; no script reads them.

## Data

Only public drug and target databases and literature were used. No individual-level sequencing data is included in this folder.
