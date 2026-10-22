#!/usr/bin/env python3
"""Rank your input variants using the outputs of alphagenome_atlas_gene_scores.py.

Keeps ONLY your input variants (your observed ALT alleles) from
  <prefix>.variants_atlas_summary.tsv   (one row per variant x target gene)
  <prefix>.position_alleles.tsv         (all 3 alleles per position)
collapses them to one row per variant, ranks them, and writes
  <prefix>.ranked_input_variants.tsv

Ranking
-------
1. Variants with an Atlas score for their ALT (atlas_status == scored) come
   first, ordered by AVI PHRED (genome-wide rank of the AVI score; higher =
   more impactful). Ties are broken by the strongest per-track effect.
2. If AVI is missing (run with --no-avi), the strongest per-track effect is
   used instead.
3. Unscored input variants (indels, outside gene windows, REF/ALT swapped, ...)
   are listed after, ordered by the strongest AVI among the 3 possible SNVs at
   their position (avi_max_phred_any_allele), as a weaker positional hint.

Added columns
-------------
  rank                      1 = highest priority
  avi_tier                  AVI PHRED band: top_0.1% (>=30) | top_1% (20-30) |
                            top_10% (10-20) | below_top_10% (<10) | no_avi
  target_genes              all target genes whose window contains the variant
  closest_target_gene, closest_dist_bp, closest_position_vs_gene
  top_scorer                per-track scorer with the strongest effect for the
                            variant (any target gene)
  top_scorer_value          that signed value (calibrated if available)
  top_scorer_track          track (tissue/cell type/assay) giving that value
  top_scorer_gene           target gene it refers to (gene-level scorers)
  n_scorers_strong          number of scorers whose |value| >= --strong
  allele_rank_at_position   1-3: rank of YOUR allele among the 3 possible SNVs
                            at the position by AVI (1 = the most damaging)
  other_alleles_avi_phred   AVI PHRED of the other two alleles, e.g. "C:12.1;T:25.3"

Usage
-----
  python rank_atlas_variants.py --prefix results/v3/atlas_scores_all_hc
  python rank_atlas_variants.py --prefix results/v3/atlas_scores_all_hc --top 50
"""

from __future__ import annotations

import argparse
import re
import sys

import numpy as np
import pandas as pd

SCORER_SUFFIX = "__max_abs_score"


def avi_tier(phred: float) -> str:
    if pd.isna(phred):
        return "no_avi"
    if phred >= 30:
        return "top_0.1%"
    if phred >= 20:
        return "top_1%"
    if phred >= 10:
        return "top_10%"
    return "below_top_10%"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", required=True, help="--out-prefix used for alphagenome_atlas_gene_scores.py")
    ap.add_argument("--out", help="output path (default <prefix>.ranked_input_variants.tsv)")
    ap.add_argument("--strong", type=float, default=0.99,
                    help="|per-scorer value| counted as strong in n_scorers_strong (default 0.99, "
                         "i.e. a calibrated quantile in the top 1%%; adjust if your scores are raw)")
    ap.add_argument("--top", type=int, default=20, help="print the top N to the screen")
    a = ap.parse_args(argv)
    out = a.out or f"{a.prefix}.ranked_input_variants.tsv"

    s = pd.read_csv(f"{a.prefix}.variants_atlas_summary.tsv", sep="\t", low_memory=False)
    try:
        pa = pd.read_csv(f"{a.prefix}.position_alleles.tsv", sep="\t", low_memory=False)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        pa = pd.DataFrame()

    scorer_cols = [c for c in s.columns if c.endswith(SCORER_SUFFIX) and not c.endswith("_any_allele")]
    scorers = [c[: -len(SCORER_SUFFIX)] for c in scorer_cols]

    # ---- strongest per-track effect per summary row (variant x target gene) ----
    if scorer_cols:
        vals = s[scorer_cols].astype(float)
        absv = vals.abs()
        has = absv.notna().any(axis=1)
        best_col = pd.Series(pd.NA, index=s.index, dtype="object")
        best_col[has] = absv[has].idxmax(axis=1)
        s["_best_scorer"] = best_col.str.replace(SCORER_SUFFIX, "", regex=False)
        s["_best_value"] = [vals.at[i, c] if isinstance(c, str) else np.nan for i, c in zip(s.index, best_col)]
        s["_best_track"] = [s.at[i, f"{sc}__top_track"] if isinstance(sc, str) and f"{sc}__top_track" in s else None
                            for i, sc in zip(s.index, s["_best_scorer"])]
        s["_n_strong"] = (absv >= a.strong).sum(axis=1)
    else:
        s["_best_scorer"], s["_best_value"], s["_best_track"], s["_n_strong"] = None, np.nan, None, 0

    # ---- collapse to one row per input variant ----------------------------------
    s["_abs_best"] = s["_best_value"].abs()
    s["_absdist"] = pd.to_numeric(s.get("dist_to_gene_bp"), errors="coerce").abs()
    # representative row = the target gene with the strongest per-track effect,
    # falling back to the closest gene
    rep = (s.sort_values(["variant_key", "_abs_best", "_absdist"], ascending=[True, False, True],
                         na_position="last")
             .drop_duplicates("variant_key").set_index("variant_key"))
    genes = s.groupby("variant_key")["target_gene"].agg(
        lambda x: ",".join(sorted({str(g) for g in x if pd.notna(g)})))
    closest = (s.sort_values(["variant_key", "_absdist"], na_position="last")
                 .drop_duplicates("variant_key").set_index("variant_key"))
    n_strong = s.groupby("variant_key")["_n_strong"].max()
    # per-scorer: strongest value across target genes
    per_scorer = {}
    for sc, col in zip(scorers, scorer_cols):
        g = s[["variant_key", col]].dropna()
        if len(g):
            idx = g[col].abs().groupby(g["variant_key"]).idxmax()
            per_scorer[col] = g.loc[idx].set_index("variant_key")[col]

    v = rep.copy()
    v["target_genes"] = genes
    v["closest_target_gene"] = closest.get("target_gene")
    v["closest_dist_bp"] = closest.get("dist_to_gene_bp")
    v["closest_position_vs_gene"] = closest.get("position_vs_gene")
    v["top_scorer"] = v["_best_scorer"]
    v["top_scorer_value"] = v["_best_value"]
    v["top_scorer_track"] = v["_best_track"]
    v["top_scorer_gene"] = v.get("target_gene")
    v["n_scorers_strong"] = n_strong
    for col, ser in per_scorer.items():
        v[col] = ser
    v = v.reset_index()

    # ---- compare your allele with the other 2 at the position --------------------
    if len(pa) and "avi_phred" in pa.columns:
        pa = pa.copy()
        pa["_r"] = pa.groupby("pos_key")["avi_phred"].rank(ascending=False, method="min")
        mine = pa[pa["atlas_alt"].notna()][["pos_key", "atlas_alt", "_r"]].rename(
            columns={"atlas_alt": "ALT", "_r": "allele_rank_at_position"})
        v = v.merge(mine, on=["pos_key", "ALT"], how="left")
        by_pos = {}
        for pk, alt, ph in pa[["pos_key", "atlas_alt", "avi_phred"]].dropna().itertuples(index=False):
            by_pos.setdefault(pk, []).append((alt, ph))
        v["other_alleles_avi_phred"] = [
            ";".join(f"{al}:{ph:.1f}" for al, ph in sorted(by_pos.get(pk, [])) if al != alt)
            for pk, alt in zip(v["pos_key"], v["ALT"])]
    else:
        v["allele_rank_at_position"] = np.nan
        v["other_alleles_avi_phred"] = ""

    # ---- rank ---------------------------------------------------------------------
    if "avi_phred" not in v:
        v["avi_phred"] = np.nan
    if "avi_max_phred_any_allele" not in v:
        v["avi_max_phred_any_allele"] = np.nan
    v["avi_tier"] = v["avi_phred"].map(avi_tier)
    scored = v["atlas_status"].eq("scored") if "atlas_status" in v else v["avi_phred"].notna()
    v["_grp"] = np.where(scored, 0, 1)
    v["_k1"] = np.where(scored, v["avi_phred"].fillna(-1), v["avi_max_phred_any_allele"].fillna(-1))
    v["_k2"] = v["top_scorer_value"].abs().fillna(-1)
    v = v.sort_values(["_grp", "_k1", "_k2"], ascending=[True, False, False]).reset_index(drop=True)
    v.insert(0, "rank", np.arange(1, len(v) + 1))

    # ---- column order -------------------------------------------------------------
    front = ["rank", "variant_key", "ID", "CHROM", "POS", "REF", "ALT", "SYMBOL", "Consequence", "IMPACT",
             "atlas_status", "avi_tier", "avi_phred", "avi_top_percent", "avi_raw",
             "avi_top_feature", "avi_top_feature_value", "allele_rank_at_position", "other_alleles_avi_phred",
             "avi_max_phred_any_allele", "avi_top_allele_any_allele",
             "top_scorer", "top_scorer_value", "top_scorer_track", "top_scorer_gene", "n_scorers_strong",
             "target_genes", "closest_target_gene", "closest_dist_bp", "closest_position_vs_gene",
             "ref_matches_atlas", "ClinVar_CLNSIG", "am_class", "am_pathogenicity", "gnomADg_AF", "MAX_AF"]
    front = [c for c in front if c in v.columns]
    drop = {"_grp", "_k1", "_k2", "_best_scorer", "_best_value", "_best_track", "_n_strong", "_abs_best",
            "_absdist", "target_gene", "target_gene_id", "dist_to_gene_bp", "position_vs_gene"}
    rest = [c for c in v.columns if c not in front and c not in drop
            and not re.search(r"(__max_abs_score_any_allele|__top_allele_any_allele)$", c)]
    v = v[front + rest]
    v.to_csv(out, sep="\t", index=False)

    # ---- screen summary -----------------------------------------------------------
    n = len(v)
    print(f"{n} input variants -> {out}")
    if "atlas_status" in v:
        print("status:", v["atlas_status"].value_counts().to_dict())
    print("AVI tiers (scored variants):", v.loc[v["atlas_status"].eq("scored") if "atlas_status" in v
                                                  else v["avi_phred"].notna(), "avi_tier"].value_counts().to_dict())
    show = [c for c in ["rank", "variant_key", "SYMBOL", "Consequence", "avi_phred", "avi_tier", "avi_top_feature",
                        "allele_rank_at_position", "top_scorer", "top_scorer_value", "top_scorer_track",
                        "closest_target_gene", "closest_dist_bp"] if c in v.columns]
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 40):
        print(v[show].head(a.top).to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())