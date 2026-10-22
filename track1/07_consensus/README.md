# List comparison

Compares the six ranked lists, flags the evidence for every candidate and supports the manual review that sets the submitted rows. The comparison was built as a workbook from the six list files; the workbook holds individual-level results and is not part of this repository.

## Inputs

| List | Approach | Pipeline step |
|---|---|---|
| a | Exomiser 15.1.0 on the DeepVariant calls | 01_alignment_qc_calling |
| b | Exomiser 15.1.0 on the provided calls | 01_alignment_qc_calling |
| c | Tiered rule-based ranking of small variants, SVs, phase, Exomiser and AlphaGenome | 05_vcf_qc_annotation_ranking |
| d | Rare variants ranked by AlphaGenome Atlas (AVI) score | 06_alphagenome_atlas |
| e | AlphaGenome scoring of rare variants within 50 kb of the 11 MVA marker genes | 04_alphagenome |
| f | ACMG point-based scoring of marker-gene variants and nearby SVs | 08_acmg_scoring |

## Steps

1. **Top 30 per list.** The first 30 distinct small variants (SNVs and indels) of each list.
   - Exomiser ranks genes, so every contributing variant of a ranked gene is a row and two variants can share a rank.
   - List d uses the first 30 rows of the rare-only file (population AF below 1% or absent). Its rank numbers come from the full list, so they are not consecutive.
2. **Match by variant.** Variants are matched on chromosome, position, REF and ALT. List f gives positions without alleles, so its variants are matched on position.
3. **Count agreement.** For each variant, record which lists have it in their top 30, and build a pairwise table of shared variants.
4. **Evidence flags.** Each variant gets six yes/no flags. Their sum (0 to 6) is used only to sort the table for review; it is not a probability and not the EPCR.
   - MVA gene: in one of the 11 marker genes.
   - Two or more independent methods: lists a and b count as one method (both Exomiser).
   - Protein-altering, or an SV/CNV.
   - Rare: max population AF below 0.1%, or no AF reported.
   - Pathogenic evidence: ClinVar pathogenic/likely pathogenic, or AlphaMissense likely pathogenic.
   - Fits recessive inheritance: homozygous, or two or more hits in the same gene.
5. **Warnings.** Flag ClinVar benign calls, AF of 1% or more, same-gene pairs in cis and annotation errors.
6. **Phase.** For same-gene pairs, look up the WhatsHap phase set and genotype in both phased call sets: same phase set with opposite orientation means trans, same orientation means cis, different phase sets means unphased.
7. **Structural variants.** The SV and CNV rows of lists c and f are compared separately. Two events match when they are on the same chromosome, have the same SV type and their intervals overlap. Read-depth exon CNVs from list c have no SV caller call, so they are reported but cannot match caller-based SVs.
8. **Manual review.** Top candidates are checked in the BAM (both call sets, allele counts by strand at MAPQ 0 and 20) and against the scoring rules. The reviewed table sets the final rows, their order and their EPCR values.

## Result

- All six lists ranked BUB1B c.2210T>G p.(Leu737Ter) and c.3006T>G p.(Asn1002Lys) first and second.
- The six top-30 lists contain 148 distinct variants: 126 (85%) in one list, 18 in two, 2 in three, none in four or five, and 2 (the BUB1B pair) in all six.
- Pairwise shared variants in the top 30:

| | a | b | c | d | e | f |
|---|---|---|---|---|---|---|
| a | 30 | 16 | 4 | 2 | 2 | 2 |
| b | 16 | 30 | 6 | 2 | 2 | 2 |
| c | 4 | 6 | 30 | 2 | 4 | 2 |
| d | 2 | 2 | 2 | 30 | 3 | 2 |
| e | 2 | 2 | 4 | 3 | 30 | 3 |
| f | 2 | 2 | 2 | 2 | 3 | 30 |

- SVs: list c has 17 distinct SV/CNV events (6 in MVA genes) and list f has 22 SVs near the panel genes. One event is in both: an inversion at chr15:38,852,774-39,425,831 near BUB1B.
- Low overlap below the top two is expected, because the approaches score different things (phenotype match, rule-based tiers, regulatory prediction, ACMG points).
