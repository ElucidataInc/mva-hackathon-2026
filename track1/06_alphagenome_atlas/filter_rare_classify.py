#!/usr/bin/env python3
"""Filter ranked AlphaGenome-Atlas variants to RARE variants and classify them
as coding / synonymous / non-coding candidates, using the VEP CSQ and ClinVar
annotations in the VCF INFO field.

Input
-----
  --ranked   <prefix>.ranked_input_variants.tsv  (from rank_atlas_variants.py)
  --vcf      the VCF used as --variants (only its header is read, to get the
             CSQ field names). Not needed if the ranked file came from an
             annotated TSV that already has Consequence / MAX_AF columns.

Output
------
  <ranked>.rare.tsv        rare variants only (rank order kept), or every
                           variant with --keep-all (use the is_rare column).
                           ALL columns of the ranked file are kept unchanged and
                           in their original order; the columns below are added
                           after them.

Added columns
-------------
  csq_<FIELD>              fields of the representative transcript for your
                           ALT: PICK=1 > MANE Select > CANONICAL > most severe
                           (Consequence, IMPACT, SYMBOL, Gene, Feature, BIOTYPE,
                           EXON, INTRON, HGVSc, HGVSp, Existing_variation, MAX_AF,
                           MAX_AF_POPS, gnomADe_AF, gnomADg_AF, CLIN_SIG,
                           SIFT, PolyPhen, am_class, am_pathogenicity,
                           SpliceAI*, CADD*, LoF* ... whatever your VEP run has)
  csq_n_transcripts        transcripts annotated for this allele
  pop_af                   population AF used for rarity (MAX_AF, else the max
                           gnomAD/1000G AF in CSQ); empty = not in those databases
  pop_af_source            column it came from
  is_rare                  pop_af < --max-af, or absent from the databases
  candidate_class          coding | synonymous | non-coding
  candidate_subclass       LoF | missense/inframe | canonical_splice |
                           synonymous_splice_region | synonymous | splice_region |
                           5'UTR | 3'UTR | intronic | upstream | downstream |
                           non-coding_transcript | regulatory | intergenic | other
  any_transcript_coding    True if ANY transcript of the allele is protein-altering
                           (e.g. intronic in MANE but missense in an isoform)
  worst_consequence        most severe consequence across all transcripts
  clinvar_stars            0-4 from CLNREVSTAT
  zygosity                 het | hom_alt | hemi/hom (from GT, single sample)

Usage
-----
  python filter_rare_classify.py \\
      --ranked results/v3/atlas_scores_all_hc.ranked_input_variants.tsv \\
      --vcf    results/v1/mva_candidate_variants.vcf.gz \\
      --max-af 0.01
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys

import numpy as np
import pandas as pd

# Ensembl VEP consequence terms, most -> least severe
SEVERITY = [
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
    "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "feature_elongation",
    "feature_truncation", "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "splice_donor_5th_base_variant", "splice_region_variant",
    "splice_donor_region_variant", "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant",
    "start_retained_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
    "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant",
    "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant", "coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification",
    "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification",
    "regulatory_region_variant", "intergenic_variant", "sequence_variant",
]
RANK = {t: i for i, t in enumerate(SEVERITY)}

LOF = {"transcript_ablation", "stop_gained", "frameshift_variant", "stop_lost", "start_lost"}
CANON_SPLICE = {"splice_acceptor_variant", "splice_donor_variant"}
MISSENSE = {"missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant",
            "incomplete_terminal_codon_variant", "feature_elongation", "feature_truncation",
            "transcript_amplification"}
PROTEIN_ALTERING = LOF | CANON_SPLICE | MISSENSE
SYNONYMOUS = {"synonymous_variant", "stop_retained_variant", "start_retained_variant"}
SPLICE_REGION = {"splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant",
                 "splice_polypyrimidine_tract_variant"}

CLNREVSTAT_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
}

AF_PREFERENCE = ["MAX_AF", "gnomADe_AF", "gnomADg_AF", "gnomAD_AF", "AF"]


# ----------------------------------------------------------------------------
def csq_fields_from_vcf(path: str) -> list[str] | None:
    opener = gzip.open if path.endswith((".gz", ".bgz")) else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            if line.startswith("##INFO=<ID=CSQ") or line.startswith("##INFO=<ID=ANN"):
                m = re.search(r"Format:\s*([^\"]+)\"", line)
                if m:
                    return [f.strip() for f in m.group(1).split("|")]
    return None


def parse_info(info) -> dict:
    out = {}
    if not isinstance(info, str) or info in ("", "."):
        return out
    for kv in info.split(";"):
        k, sep, v = kv.partition("=")
        out[k] = v if sep else "True"
    return out


def vep_allele(ref: str, alt: str) -> str:
    """VEP writes indel alleles without the shared anchor base ('-' if empty)."""
    ref, alt = str(ref), str(alt)
    if len(ref) == 1 and len(alt) == 1:
        return alt
    if ref and alt and ref[0] == alt[0]:
        return alt[1:] or "-"
    return alt


def num(v):
    """Largest numeric value in a VEP field ('&'-separated, '.'/'' = missing)."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    vals = []
    for p in re.split(r"[&,]", str(v)):
        try:
            vals.append(float(p))
        except ValueError:
            pass
    return max(vals) if vals else np.nan


def worst(cons_list) -> str:
    terms = [t for c in cons_list for t in str(c).split("&") if t]
    return min(terms, key=lambda t: RANK.get(t, len(SEVERITY))) if terms else ""


def classify(consequence: str) -> tuple[str, str]:
    terms = set(str(consequence).split("&")) if consequence else set()
    if terms & LOF:
        return "coding", "LoF"
    if terms & CANON_SPLICE:
        return "coding", "canonical_splice"
    if terms & MISSENSE:
        return "coding", "missense/inframe"
    if terms & SYNONYMOUS:
        return "synonymous", "synonymous_splice_region" if terms & SPLICE_REGION else "synonymous"
    if terms & SPLICE_REGION:
        return "non-coding", "splice_region"
    for t, sub in (("5_prime_UTR_variant", "5'UTR"), ("3_prime_UTR_variant", "3'UTR"),
                   ("mature_miRNA_variant", "non-coding_transcript"),
                   ("non_coding_transcript_exon_variant", "non-coding_transcript"),
                   ("intron_variant", "intronic"), ("upstream_gene_variant", "upstream"),
                   ("downstream_gene_variant", "downstream"), ("TF_binding_site_variant", "regulatory"),
                   ("regulatory_region_variant", "regulatory"), ("intergenic_variant", "intergenic"),
                   ("non_coding_transcript_variant", "non-coding_transcript")):
        if t in terms:
            return "non-coding", sub
    return ("non-coding", "other") if terms else ("unknown", "no_annotation")


def pick_block(blocks: list[dict]) -> dict:
    for test in (lambda b: b.get("PICK") == "1",
                 lambda b: bool(b.get("MANE_SELECT")) or b.get("MANE") == "MANE_Select",
                 lambda b: b.get("CANONICAL") == "YES"):
        hit = [b for b in blocks if test(b)]
        if hit:
            return min(hit, key=lambda b: RANK.get(worst([b.get("Consequence", "")]), 99))
    return min(blocks, key=lambda b: RANK.get(worst([b.get("Consequence", "")]), 99))


def zygosity(gt) -> str:
    if not isinstance(gt, str) or gt in ("", "."):
        return ""
    al = re.split(r"[/|]", gt)
    if any(a == "." for a in al):
        return "partial/missing"
    al = [a for a in al if a != "0"]
    if len(re.split(r"[/|]", gt)) == 1:
        return "hemi/hom" if al else "ref"
    if not al:
        return "ref"
    return "hom_alt" if len(al) == 2 and al[0] == al[1] else "het"


# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ranked", required=True, help="<prefix>.ranked_input_variants.tsv")
    ap.add_argument("--vcf", help="VCF whose header defines the CSQ field order")
    ap.add_argument("--csq-format", help="CSQ field names 'Allele|Consequence|...' if no --vcf")
    ap.add_argument("--max-af", type=float, default=0.01,
                    help="rare = population AF below this (default 0.01 = 1%%; use 0.001 for dominant / "
                         "ultra-rare disease)")
    ap.add_argument("--exclude-af-missing", action="store_true",
                    help="treat variants absent from gnomAD/1000G as NOT rare (default: rare)")
    ap.add_argument("--keep-all", action="store_true", help="write all variants with flags, not only rare")
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    out = a.out or re.sub(r"\.tsv$", "", a.ranked) + (".annotated.tsv" if a.keep_all else ".rare.tsv")

    df = pd.read_csv(a.ranked, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
    info_col = "INFO" if "INFO" in df else None
    infos = [parse_info(x) for x in df[info_col]] if info_col else [{} for _ in range(len(df))]
    csq_key = next((k for k in ("CSQ", "ANN") if any(k in d for d in infos)), None)
    if csq_key is None and "INFO_CSQ" in df:
        csq_key = "CSQ"
        for d, v in zip(infos, df["INFO_CSQ"]):
            d["CSQ"] = v

    fields = None
    if csq_key:
        fields = (a.csq_format.split("|") if a.csq_format else csq_fields_from_vcf(a.vcf) if a.vcf else None)
        if fields is None:
            sys.exit("INFO has CSQ but its field names are unknown: pass --vcf <your.vcf.gz> "
                     "(header is read) or --csq-format 'Allele|Consequence|...'")

    # ---- per-variant CSQ for YOUR allele -------------------------------------------
    rows = []
    for i, (d, ref, alt) in enumerate(zip(infos, df["REF"], df["ALT"])):
        rec = {}
        if csq_key and d.get(csq_key):
            blocks = [dict(zip(fields, b.split("|"))) for b in d[csq_key].split(",")]
            al = vep_allele(ref, alt)
            mine = [b for b in blocks if b.get("Allele") == al] or blocks
            p = pick_block(mine)
            rec = {f"csq_{k}": v for k, v in p.items() if k != "Allele"}
            rec["csq_n_transcripts"] = len(mine)
            rec["worst_consequence"] = worst([b.get("Consequence", "") for b in mine])
            rec["any_transcript_coding"] = any(
                set(b.get("Consequence", "").split("&")) & PROTEIN_ALTERING for b in mine)
            cons = p.get("Consequence", "")
        else:  # annotated TSV input: use its own columns
            cons = df.at[i, "Consequence"] if "Consequence" in df else ""
            rec["worst_consequence"] = worst([cons])
            rec["any_transcript_coding"] = bool(set(cons.split("&")) & PROTEIN_ALTERING)
        rec["candidate_class"], rec["candidate_subclass"] = classify(cons)

        # population AF: CSQ first, then same-named columns of an annotated TSV
        af, src = np.nan, ""
        for name in AF_PREFERENCE:
            for col, val in ((f"csq_{name}", rec.get(f"csq_{name}")), (name, df.at[i, name] if name in df else None)):
                if name == "AF" and col == "AF":
                    continue  # top-level AF is the call-set AF, never population AF
                v = num(val)
                if not np.isnan(v):
                    af, src = v, col
                    break
            if src:
                break
        if not src:  # any other gnomAD AF field
            cand = [(k, num(v)) for k, v in rec.items() if re.match(r"csq_gnomAD[eg]?_.*AF$", k)]
            cand = [(k, v) for k, v in cand if not np.isnan(v)]
            if cand:
                src, af = max(cand, key=lambda kv: kv[1])
        rec["pop_af"], rec["pop_af_source"] = af, src

        # ClinVar stars
        rs = d.get("CLNREVSTAT") or (df.at[i, "ClinVar_CLNREVSTAT"] if "ClinVar_CLNREVSTAT" in df else "")
        rec["clinvar_stars"] = CLNREVSTAT_STARS.get(str(rs), 0 if rs else "")
        rec["zygosity"] = zygosity(df.at[i, "GT"]) if "GT" in df else ""
        rows.append(rec)

    ann = pd.DataFrame(rows, index=df.index)
    absent = ann["pop_af"].isna()
    ann["is_rare"] = (ann["pop_af"] < a.max_af) | (absent & (not a.exclude_af_missing))
    res = pd.concat([df, ann[[c for c in ann.columns if c not in df.columns]]], axis=1)

    # Column order: the ranked file exactly as it was (all columns, same order,
    # values untouched), then the new classification columns, then csq_* fields.
    key = ["is_rare", "pop_af", "pop_af_source", "candidate_class", "candidate_subclass",
           "any_transcript_coding", "worst_consequence", "clinvar_stars", "zygosity"]
    new = [c for c in res.columns if c not in df.columns]
    key = [c for c in key if c in new]
    res = res[list(df.columns) + key + [c for c in new if c not in key]]

    kept = res if a.keep_all else res[res["is_rare"]]
    kept.to_csv(out, sep="\t", index=False)

    n = len(res)
    print(f"{n} variants; rare (pop AF < {a.max_af}{'' if a.exclude_af_missing else ' or absent'}): "
          f"{int(res['is_rare'].sum())}  (absent from AF databases: {int(absent.sum())})")
    print("AF source:", ann["pop_af_source"].replace("", "none").value_counts().to_dict())
    rare = res[res["is_rare"]]
    print("rare by class:", rare["candidate_class"].value_counts().to_dict())
    print("rare by subclass:", rare["candidate_subclass"].value_counts().to_dict())
    print(f"-> {out}")


if __name__ == "__main__":
    sys.exit(main())