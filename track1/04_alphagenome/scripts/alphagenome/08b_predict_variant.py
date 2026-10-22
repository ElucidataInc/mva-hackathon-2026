#!/usr/bin/env python3
"""08b: AlphaGenome predict_variant() for selected rare variants -> REF/ALT tracks, alt-ref summary, plots.
Sends variant coordinates to the AlphaGenome API (same approval as 07; logged). Outputs are server-only."""
import os, argparse, datetime, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/alphagenome/v2_team11/alphagenome_model"))
SUMM = f"{BASE}/rare_alphagenome_summary.tsv"; OUT = f"{BASE}/predict"; LOG = f"{BASE}/logs/08b_predict.log"
TERMS_DEFAULT = "CL:0011020,NTR:0000427,CL:0000047,UBERON:0000955,EFO:0003042,CL:2000001,CL:0002551,EFO:0002067"
MODS = {"RNA_SEQ": "rna_seq", "SPLICE_SITES": "splice_sites", "SPLICE_SITE_USAGE": "splice_site_usage",
        "SPLICE_JUNCTIONS": "splice_junctions", "DNASE": "dnase", "ATAC": "atac", "CAGE": "cage", "CHIP_HISTONE": "chip_histone"}
PLOT = ["rna_seq", "splice_site_usage", "dnase", "chip_histone"]

ap = argparse.ArgumentParser()
ap.add_argument("--approved-by", required=True)
ap.add_argument("--tiers", default="A_coding_or_clinvar,B_noncoding_AG_supported")
ap.add_argument("--max-variants", type=int, default=0)
ap.add_argument("--terms", default=TERMS_DEFAULT)
ap.add_argument("--crop-bp", type=int, default=20000)
a = ap.parse_args()
os.makedirs(OUT, exist_ok=True); log = open(LOG, "a")
def say(m): print(m, flush=True); log.write(m + "\n"); log.flush()

s = pd.read_csv(SUMM, sep="\t", dtype=str, low_memory=False)
s = s[s.final_tier.isin(a.tiers.split(","))]
if a.max_variants: s = s.head(a.max_variants)
terms = a.terms.split(",")
say(f"\n== {datetime.datetime.now():%Y-%m-%d %H:%M} approved_by='{a.approved_by}' n_variants={len(s)} terms={terms}")
model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])

allsum = []
for i, r in enumerate(s.itertuples(), 1):
    d = f"{OUT}/{r.fid}"; os.makedirs(d, exist_ok=True)
    if os.path.exists(f"{d}/diff_summary.tsv"):
        say(f"[{i}/{len(s)}] {r.marker_gene} {r.final_tier[:1]} {r.consequence}: already done")
        allsum.append(pd.read_csv(f"{d}/diff_summary.tsv", sep="\t").assign(fid=r.fid)); continue
    var = genome.Variant(chromosome=r.chrom, position=int(r.pos), reference_bases=r.ref, alternate_bases=r.alt)
    iv = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    try:
        out = model.predict_variant(interval=iv, variant=var, ontology_terms=terms,
                                    requested_outputs=[getattr(dna_client.OutputType, k) for k in MODS])
    except Exception as e:
        say(f"[{i}/{len(s)}] {r.marker_gene}: FAILED ({str(e)[:150]})"); continue
    rows = []
    for k, attr in MODS.items():
        R, A = getattr(out.reference, attr, None), getattr(out.alternate, attr, None)
        if R is None or A is None or getattr(R, "values", None) is None or R.values.size == 0: continue
        md = R.metadata.copy(); md.to_csv(f"{d}/metadata_{attr}.tsv", sep="\t", index=False)
        rv, av = np.asarray(R.values), np.asarray(A.values)
        if attr == "splice_junctions" or rv.ndim != 2 or not hasattr(R, "resolution"):
            np.savez_compressed(f"{d}/tracks_{attr}.npz", ref=rv, alt=av)
            diff = (av - rv).reshape(len(rv), -1) if rv.ndim > 1 else (av - rv)[:, None]
            for t in range(diff.shape[1] if diff.ndim == 2 else 0):
                j = int(np.abs(diff[:, t]).argmax())
                rows.append({"modality": attr, "track": t, "max_abs_diff": float(abs(diff[j, t])), "signed_diff": float(diff[j, t]),
                             "offset_bp": None, "ref_sum": float(rv.reshape(len(rv), -1)[:, t].sum()), "alt_sum": float(av.reshape(len(av), -1)[:, t].sum())})
            continue
        res = int(R.resolution); c = (int(r.pos) - R.interval.start) // res; w = a.crop_bp // res
        lo, hi = max(0, c - w), min(rv.shape[0], c + w + 1)
        rc, ac = rv[lo:hi], av[lo:hi]
        np.savez_compressed(f"{d}/tracks_{attr}.npz", ref=rc, alt=ac, resolution=res,
                            start=R.interval.start + lo * res, chrom=R.interval.chromosome)
        diff = ac - rc
        for t in range(diff.shape[1]):
            j = int(np.abs(diff[:, t]).argmax())
            rows.append({"modality": attr, "track": t, "max_abs_diff": float(abs(diff[j, t])), "signed_diff": float(diff[j, t]),
                         "offset_bp": int((lo + j - c) * res), "ref_sum": float(rc[:, t].sum()), "alt_sum": float(ac[:, t].sum())})
        if attr in PLOT:
            try:
                fig = plot_components.plot(
                    [plot_components.OverlaidTracks(tdata={"REF": R, "ALT": A}, colors={"REF": "dimgrey", "ALT": "red"})],
                    interval=var.reference_interval.resize(2**14),
                    annotations=[plot_components.VariantAnnotation([var], alpha=0.8)])
                fig.savefig(f"{d}/plot_{attr}.png", dpi=110, bbox_inches="tight"); plt.close(fig)
            except Exception as e:
                say(f"    plot {attr} skipped ({str(e)[:100]})")
    ds = pd.DataFrame(rows)
    if len(ds):
        for attr in ds.modality.unique():
            mdf = f"{d}/metadata_{attr}.tsv"
            if os.path.exists(mdf):
                m = pd.read_csv(mdf, sep="\t"); cols = [c for c in ["name", "biosample_name", "Assay title", "strand"] if c in m.columns]
                lab = m[cols].astype(str).agg("|".join, axis=1) if cols else pd.Series(range(len(m))).astype(str)
                ds.loc[ds.modality == attr, "track_label"] = ds.loc[ds.modality == attr, "track"].map(lambda t: lab.iloc[t] if t < len(lab) else "")
    ds.to_csv(f"{d}/diff_summary.tsv", sep="\t", index=False)
    json.dump({"marker_gene": r.marker_gene, "tier": r.final_tier, "consequence": r.consequence, "fid": r.fid, "terms": terms}, open(f"{d}/info.json", "w"))
    allsum.append(ds.assign(fid=r.fid))
    say(f"[{i}/{len(s)}] {r.marker_gene} {r.final_tier[:1]} {r.consequence}: ok ({len(ds)} tracks summarised)")

if allsum:
    A = pd.concat(allsum, ignore_index=True)
    top = A.loc[A.groupby(["fid", "modality"]).max_abs_diff.idxmax()]
    top = top.merge(s[["fid", "marker_gene", "final_tier", "consequence"]], on="fid", how="left")
    top.to_csv(f"{OUT}/predict_top_per_variant_modality.tsv", sep="\t", index=False)
    say("\n== strongest alt-ref change per variant x modality (no coordinates) ==")
    say(top[["marker_gene", "final_tier", "consequence", "modality", "max_abs_diff", "signed_diff", "offset_bp", "track_label"]]
        .sort_values(["final_tier", "marker_gene", "modality"]).to_string(index=False))
