#!/usr/bin/env python3
"""08c: re-summarise saved predict_variant tracks LOCALLY (+/-500 bp around variant): log2(ALT/REF) of window sums.
Reads saved .npz only - no API call. Avoids bfloat16 rounding noise from far peaks."""
import glob, json, os, numpy as np, pandas as pd
OUT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/alphagenome/v2_team11/alphagenome_model/predict"))
HALF = 500
rows = []
for d in sorted(glob.glob(f"{OUT}/*/")):
    if not os.path.exists(d + "info.json"): continue
    info = json.load(open(d + "info.json"))
    for f in sorted(glob.glob(d + "tracks_*.npz")):
        attr = os.path.basename(f)[7:-4]; z = np.load(f, allow_pickle=True)
        if "resolution" not in z.files: continue          # junctions: not positional
        ref, alt, res = z["ref"].astype(float), z["alt"].astype(float), int(z["resolution"])
        c = ref.shape[0] // 2; w = max(1, HALF // res); lo, hi = max(0, c - w), c + w + 1
        md = pd.read_csv(d + f"metadata_{attr}.tsv", sep="\t")
        eps = 1e-3 if attr.startswith("splice") else 1.0
        for t in range(ref.shape[1]):
            rs, as_ = ref[lo:hi, t].sum(), alt[lo:hi, t].sum()
            rows.append({"marker_gene": info["marker_gene"], "consequence": info["consequence"], "tier": info["tier"][:1],
                         "modality": attr, "biosample": md["biosample_name"].iloc[t] if "biosample_name" in md else "",
                         "track": md["name"].iloc[t] if "name" in md else t,
                         "local_ref": round(rs, 3), "local_alt": round(as_, 3),
                         "local_absdiff": round(abs(as_ - rs), 4), "local_log2": round(float(np.log2((as_ + eps) / (rs + eps))), 4)})
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/predict_local_summary.tsv", sep="\t", index=False)
df["a"] = df.local_log2.abs()
top = df.loc[df.groupby(["marker_gene", "consequence", "modality"]).a.idxmax()]
pd.set_option("display.width", 200)
print(f"== strongest LOCAL change (+/-{HALF} bp) per variant x modality ==")
print(top[["marker_gene", "tier", "consequence", "modality", "biosample", "local_ref", "local_alt", "local_log2", "local_absdiff"]].to_string(index=False))
