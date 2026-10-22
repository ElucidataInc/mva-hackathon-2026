#!/usr/bin/env python3
"""07: Score rare marker-region variants with the AlphaGenome API, then summarise per variant.
NOTE: sends variant coordinates to the AlphaGenome API -> run only with team approval (logged).
Full per-variant scores are server-only (delete at hackathon close)."""
import os, time, argparse, datetime, hashlib
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/alphagenome/v2_team11"))
IN   = f"{BASE}/marker_region_vep.tsv"
OUT  = f"{BASE}/alphagenome_model"
SC   = f"{OUT}/scores"
LOG  = f"{OUT}/logs/07_score.log"
SUMM = f"{OUT}/rare_alphagenome_summary.tsv"

ap = argparse.ArgumentParser()
ap.add_argument("--approved-by", required=True, help="who approved sending variants to the API, and date")
ap.add_argument("--ranks", default="1,2,3,4,5", help="rank_key values to score")
ap.add_argument("--max-variants", type=int, default=0, help="limit (pilot run)")
ap.add_argument("--summarise-only", action="store_true")
a = ap.parse_args()
os.makedirs(SC, exist_ok=True); os.makedirs(f"{OUT}/logs", exist_ok=True)
log = open(LOG, "a")
def say(m): print(m, flush=True); log.write(m + "\n"); log.flush()

# ---- rare variants + rank_key (same rule as rare_shortlist) ----
v = pd.read_csv(IN, sep="\t", dtype=str, low_memory=False)
v = v[v["rare"] == "yes"].copy()
def rank(r):
    if r["impact"] == "HIGH": return 1
    if r["impact"] == "MODERATE": return 2
    if str(r["clinvar_exact_clnsig"]).startswith(("Pathogenic", "Likely_pathogenic")) or r["am_class"] == "likely_pathogenic": return 3
    if r["MAX_AF"] == ".": return 4
    return 5
v["rank_key"] = v.apply(rank, axis=1)
v = v[v["rank_key"].isin({int(x) for x in a.ranks.split(",")})]
v = v[~v["alt"].isin(["*", "."])].sort_values(["rank_key", "marker_gene"])
if a.max_variants: v = v.head(a.max_variants)
v["key"] = v.chrom + ":" + v.pos + ":" + v.ref + ":" + v.alt
v["fid"] = v.key.map(lambda k: hashlib.md5(k.encode()).hexdigest()[:12])
say(f"\n== {datetime.datetime.now():%Y-%m-%d %H:%M} approved_by='{a.approved_by}' ranks={a.ranks} n_variants={len(v)}")

# ---- scoring (resumable: skips variants already scored) ----
if not a.summarise_only:
    model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
    scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
    for i, r in enumerate(v.itertuples(), 1):
        fn = f"{SC}/{r.fid}.tsv.gz"
        if os.path.exists(fn): say(f"[{i}/{len(v)}] {r.marker_gene} rank{r.rank_key}: already scored"); continue
        var = genome.Variant(chromosome=r.chrom, position=int(r.pos), reference_bases=r.ref, alternate_bases=r.alt)
        iv = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        for attempt in range(4):
            try:
                s = model.score_variant(interval=iv, variant=var, variant_scorers=scorers); break
            except Exception as e:
                err = str(e)[:150]; time.sleep(20 * (attempt + 1))
        else:
            say(f"[{i}/{len(v)}] {r.marker_gene} rank{r.rank_key}: FAILED ({err})"); continue
        df = variant_scorers.tidy_scores([s])
        df.to_csv(fn, sep="\t", index=False, compression="gzip")
        say(f"[{i}/{len(v)}] {r.marker_gene} rank{r.rank_key}: ok ({len(df)} rows)")

# ---- per-variant summary: EFFECT scorers only ----
def summarise(df, gene):
    df = df.copy(); df["aq"] = df.quantile_score.abs(); sc = df.variant_scorer.astype(str)
    out = {}
    def top(sub, p):
        if len(sub) == 0:
            out.update({f"{p}_q": None, f"{p}_raw": None, f"{p}_where": ""}); return
        t = sub.loc[sub.aq.idxmax()]
        out.update({f"{p}_q": round(t.quantile_score, 4), f"{p}_raw": round(t.raw_score, 4),
                    f"{p}_where": f"{t.biosample_name}|{t['Assay title']}"})
    rna = df[(df.output_type == "RNA_SEQ") & sc.str.startswith("GeneMaskLFCScorer")]
    own = rna[rna.gene_name == gene]
    top(own, "expr_own"); out["expr_own_ntracks_q95"] = int((own.aq > 0.95).sum())
    oth = rna[rna.gene_name != gene]
    top(oth, "expr_other"); out["expr_other_gene"] = oth.loc[oth.aq.idxmax(), "gene_name"] if len(oth) else ""
    for ot, p in [("SPLICE_SITES", "splice_site"), ("SPLICE_SITE_USAGE", "splice_usage"), ("SPLICE_JUNCTIONS", "splice_junction")]:
        sub = df[df.output_type == ot]; g = sub[sub.gene_name == gene]
        top(g if len(g) else sub, p)
    for ot in ["DNASE", "ATAC", "CAGE", "PROCAP", "CHIP_TF", "CHIP_HISTONE"]:
        top(df[(df.output_type == ot) & sc.str.contains("DIFF_LOG2_SUM")], ot.lower())
    return out

rows = []
for r in v.itertuples():
    fn = f"{SC}/{r.fid}.tsv.gz"
    if not os.path.exists(fn): continue
    d = pd.read_csv(fn, sep="\t", low_memory=False)
    rows.append({"key": r.key, "fid": r.fid, **summarise(d, r.marker_gene)})
ag = pd.DataFrame(rows)
res = v.merge(ag, on=["key", "fid"], how="left")

# heuristic evidence flags (thresholds are adjustable; documented in log)
f = lambda c: pd.to_numeric(res[c], errors="coerce").fillna(0)
res["ag_expression"] = (f("expr_own_q").abs() >= 0.95) & (f("expr_own_raw").abs() >= 0.1)
res["ag_splicing"]   = ((f("splice_site_q") >= 0.95) & (f("splice_site_raw") >= 0.2)) | \
                       ((f("splice_usage_q").abs() >= 0.95) & (f("splice_usage_raw").abs() >= 0.1)) | \
                       ((f("splice_junction_q").abs() >= 0.99) & (f("splice_junction_raw").abs() >= 0.5))
reg = [c for c in ["dnase", "atac", "cage", "procap", "chip_tf", "chip_histone"]]
res["ag_regulatory"] = pd.concat([(f(f"{c}_q").abs() >= 0.99) & (f(f"{c}_raw").abs() >= 0.5) for c in reg], axis=1).any(axis=1)
res["ag_scored"] = res["expr_own_q"].notna() | res["splice_site_q"].notna()
res["ag_class"] = res.apply(lambda r: ";".join([n for n in ["splicing", "expression", "regulatory"] if r[f"ag_{n}"]]) or ("none" if r["ag_scored"] else "not_scored"), axis=1)
res["ag_max_q"] = pd.concat([f("expr_own_q").abs(), f("splice_site_q").abs(), f("splice_usage_q").abs()], axis=1).max(axis=1)
# final order: coding/ClinVar first (rank 1-3), then AlphaGenome-supported, then rest
res["final_tier"] = res.apply(lambda r: "A_coding_or_clinvar" if r.rank_key <= 3 else ("B_noncoding_AG_supported" if r.ag_class not in ("none", "not_scored") else "C_noncoding_no_AG_signal"), axis=1)
res = res.sort_values(["final_tier", "rank_key", "ag_max_q"], ascending=[True, True, False])
res.to_csv(SUMM, sep="\t", index=False)

say("thresholds: expression |q|>=0.95 & |log2FC|>=0.1; splicing site q>=0.95 & raw>=0.2 OR usage |q|>=0.95 & |raw|>=0.1 OR junction |q|>=0.99 & |raw|>=0.5; regulatory |q|>=0.99 & |log2|>=0.5")
say(f"scored: {int(res.ag_scored.sum())}/{len(res)}  -> {SUMM} (server only)")
say("\n== per gene (counts only) ==")
say(res.groupby("marker_gene").agg(rare=("key", "size"), scored=("ag_scored", "sum"),
    expression=("ag_expression", "sum"), splicing=("ag_splicing", "sum"), regulatory=("ag_regulatory", "sum")).to_string())
say("\n== final tiers ==\n" + res.final_tier.value_counts().to_string())
