#!/usr/bin/env python3
"""10: add readable splice-site + expression columns to the FINAL table, from our raw AlphaGenome
score files (scores/<fid>.tsv.gz) + VEP fields. No API call. Writes a NEW file."""
import os, re, pandas as pd
AG = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/alphagenome/v2_team11/alphagenome_model"))
F = pd.read_csv(f"{AG}/FINAL_alphagenome_variant_table.tsv", sep="\t", low_memory=False)
OUT = f"{AG}/FINAL_alphagenome_variant_table_splice_expression.tsv"
pct = lambda l2: round((2 ** l2 - 1) * 100, 1)
rna_scorers = set(); rows = []
for r in F.itertuples():
    o = {"fid": r.fid}; g = r.marker_gene
    o["vep_splice_region"] = "yes" if "splice" in str(r.consequence) else "no"
    m = re.search(r"c\.[*-]?\d+([+-]\d+)", str(r.hgvsc)); o["intron_offset_bp"] = int(m.group(1)) if m else None
    d = pd.read_csv(f"{AG}/scores/{r.fid}.tsv.gz", sep="\t", low_memory=False)
    d["aq"] = d.quantile_score.abs(); sc = d.variant_scorer.astype(str)
    rna_scorers |= set(sc[d.output_type == "RNA_SEQ"].str.split("(").str[0])
    def own(ot):
        x = d[d.output_type == ot]; return x[x.gene_name == g] if (x.gene_name == g).any() else x
    site, use, jun = 0, 0, 0
    ss = d[d.output_type == "SPLICE_SITES"]
    if len(ss):
        t = ss.loc[ss.aq.idxmax()]; site = abs(t.raw_score)
        o["ag_splice_site_change"] = f"{t.track_name} ({t.track_strand}) {t.raw_score:+.3f}"
    su = own("SPLICE_SITE_USAGE")
    if len(su):
        t = su.loc[su.aq.idxmax()]; use = abs(t.raw_score); o["ag_splice_usage_change"] = f"{t.raw_score:+.3f} ({t.biosample_name})"
    sj = own("SPLICE_JUNCTIONS")
    if len(sj):
        t = sj.loc[sj.aq.idxmax()]; jun = abs(t.raw_score); o["ag_junction_change"] = f"{t.raw_score:+.3f} ({t.biosample_name})"
    o["ag_splice_verdict"] = ("strong" if site >= 0.5 or jun >= 2 else
                              "moderate" if site >= 0.2 or use >= 0.1 or jun >= 1 else
                              "weak" if site >= 0.1 or use >= 0.05 or jun >= 0.5 else "none")
    lfc = d[(d.output_type == "RNA_SEQ") & sc.str.startswith("GeneMaskLFCScorer") & (d.gene_name == g)]
    if len(lfc):
        top = lfc.loc[lfc.aq.idxmax()]
        o["ag_expr_direction"] = "up" if top.raw_score > 0 else "down"
        o["ag_expr_change_pct_top"] = pct(top.raw_score)
        o["ag_expr_median_log2fc"] = round(lfc.raw_score.median(), 4)
        o["ag_expr_n_tracks"] = len(lfc); o["ag_expr_n_tracks_q95"] = int((lfc.aq > 0.95).sum())
        t3 = lfc.sort_values("aq", ascending=False).head(3)
        o["ag_expr_top3"] = "; ".join(f"{b} {pct(x):+.1f}%" for b, x in zip(t3.biosample_name, t3.raw_score))
        ch = abs(pct(top.raw_score)) if top.aq >= 0.95 else 0
        o["ag_expr_verdict"] = "strong" if ch >= 30 else "moderate" if ch >= 10 else "small" if ch >= 2 else "none"
    else:
        o["ag_expr_verdict"] = f"no {g} RNA rows (gene not in scored window/annotation)"
    act = d[(d.output_type == "RNA_SEQ") & sc.str.startswith("GeneMask") & ~sc.str.startswith("GeneMaskLFC") & (d.gene_name == g)]
    if len(act):
        o["ag_gene_baseline_activity_median"] = round(act.raw_score.median(), 3)
        o["ag_gene_most_active_in"] = act.loc[act.raw_score.idxmax(), "biosample_name"]
    rows.append(o)
N = F.merge(pd.DataFrame(rows), on="fid", how="left")
N.to_csv(OUT, sep="\t", index=False)
print("RNA scorers seen:", sorted(rna_scorers))
print(f"rows: {len(N)} -> {OUT} (server only)")
print("\n== splice verdict ==\n" + N.ag_splice_verdict.value_counts().to_string())
print("\n== expression verdict ==\n" + N.ag_expr_verdict.value_counts().to_string())
print("\n== tier A + any non-'none' (no coordinates) ==")
x = N[(N.final_tier.str.startswith("A")) | (N.ag_splice_verdict != "none") | (~N.ag_expr_verdict.isin(["none"]))]
print(x[["marker_gene", "consequence", "quality_ok", "vep_splice_region", "intron_offset_bp", "ag_splice_site_change",
         "ag_splice_verdict", "ag_expr_direction", "ag_expr_change_pct_top", "ag_expr_n_tracks_q95", "ag_expr_verdict"]].to_string(index=False))
