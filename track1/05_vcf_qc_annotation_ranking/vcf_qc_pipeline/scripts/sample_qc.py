#!/usr/bin/env python3
"""
Sample-level QC from the profiler's outputs (gt_matrix.npz, per_sample.tsv, profile.json):
  * KING-robust kinship + IBS0 -> relatedness / duplicate detection
  * within-cohort PCA (ancestry structure; reference projection comes from somalier if run)
  * sex check: chrX non-PAR heterozygosity and chrY call rate vs declared sex
  * contamination screen: autosomal het-rate outliers and skewed het allele balance
  * runs of homozygosity (bcftools roh output if supplied, otherwise a fallback scan)
  * outlier flags on per-sample metrics
"""

import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcf_common import read_sex_file, write_json, is_autosome


def mad_flags(vals, k=4.0, min_rel=0.15, min_abs=0.0):
    """Robust outlier flags: |x - median| > k*MAD and also > min_rel*|median| and > min_abs
    (the extra floors stop tiny absolute differences being flagged in homogeneous cohorts)."""
    v = np.array([np.nan if x is None else x for x in vals], dtype=float)
    if len(v) < 4 or not np.isfinite(v).any():
        return [False] * len(v), np.nan, np.nan
    med = np.nanmedian(v)
    mad = np.nanmedian(np.abs(v - med)) * 1.4826
    if not np.isfinite(mad) or mad == 0:
        return [False] * len(v), med, mad
    thr = max(k * mad, min_rel * abs(med), min_abs)
    return [bool(abs(x - med) > thr) if np.isfinite(x) else False for x in v], med, mad


def king_robust(G):
    """G: sites x samples int8 (0/1/2, -1 missing). Returns kinship, ibs0, n_shared matrices."""
    C = (G >= 0).astype(np.float32)
    H = (G == 1).astype(np.float32)
    A0 = (G == 0).astype(np.float32)
    A2 = (G == 2).astype(np.float32)
    nhh = H.T @ H
    nopp = A0.T @ A2 + A2.T @ A0
    hc = H.T @ C                       # het in i at sites where j is called
    denom = hc + hc.T
    with np.errstate(divide="ignore", invalid="ignore"):
        kin = (nhh - 2 * nopp) / denom
        nboth = C.T @ C
        ibs0 = nopp / nboth
    return kin, ibs0, nboth


def relation(k, ibs0):
    if k > 0.354:
        return "duplicate / MZ twin"
    if k > 0.177:
        return "1st degree (parent-offspring)" if ibs0 < 0.005 else "1st degree (full sibling)"
    if k > 0.0884:
        return "2nd degree"
    if k > 0.0442:
        return "3rd degree"
    return "unrelated"


def pca(G, n_pc=10):
    X = G.astype(float)
    X[G < 0] = np.nan
    p = np.nanmean(X, axis=1) / 2
    keep = (p > 0.01) & (p < 0.99)
    X, p = X[keep], p[keep]
    if X.shape[0] < 50 or X.shape[1] < 3:
        return None
    mu = 2 * p
    X = np.where(np.isnan(X), mu[:, None], X)
    sd = np.sqrt(2 * p * (1 - p))
    Z = (X - mu[:, None]) / sd[:, None]
    U, S, Vt = np.linalg.svd(Z.T, full_matrices=False)
    k = min(n_pc, len(S))
    var = S ** 2 / np.sum(S ** 2)
    return {"pcs": (U[:, :k] * S[:k]).tolist(), "var_explained": var[:k].tolist(), "n_sites": int(X.shape[0])}


def roh_fallback(G, chrom, pos, samples, min_sites=25, min_len=1_000_000):
    """Scan the down-sampled dosage matrix for runs of homozygous calls (0 or 2)."""
    res = {s: [] for s in samples}
    chroms = np.unique(chrom)
    for ch in chroms:
        if not is_autosome(str(ch)):
            continue
        m = chrom == ch
        g, p = G[m], pos[m]
        for j, s in enumerate(samples):
            col = g[:, j]
            start = None
            n = 0
            last = None
            for k in range(len(col)):
                v = col[k]
                if v < 0:
                    continue
                if v == 1:
                    if start is not None and n >= min_sites and p[last] - p[start] >= min_len:
                        res[s].append((str(ch), int(p[start]), int(p[last]), n))
                    start, n = None, 0
                else:
                    if start is None:
                        start = k
                    n += 1
                    last = k
            if start is not None and n >= min_sites and p[last] - p[start] >= min_len:
                res[s].append((str(ch), int(p[start]), int(p[last]), n))
    return res


def read_bcftools_roh(path):
    res = {}
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        for line in fh:
            if line.startswith("RG"):
                p = line.rstrip().split("\t")
                # RG sample chrom start end length nmarkers quality
                s, ch, st, en, ln, nm = p[1], p[2], int(p[3]), int(p[4]), int(p[5]), int(p[6])
                res.setdefault(s, []).append((ch, st, en, nm))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sex-file", help="sample<TAB>sex or PED")
    ap.add_argument("--roh", help="bcftools roh output (-O r)")
    ap.add_argument("--genome-size-mb", type=float, default=2875.0, help="autosomal size for F_ROH")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    prof = json.load(open(os.path.join(a.profile_dir, "profile.json")))
    samples = prof["samples"]
    rows = list(csv.DictReader(open(os.path.join(a.profile_dir, "per_sample.tsv")), delimiter="\t")) if samples else []
    fl = lambda r, k: float(r[k]) if r.get(k) not in (None, "") else None
    out = {"n_samples": len(samples)}
    if not rows:
        write_json(out, os.path.join(a.out, "sample_qc.json"))
        return

    # ---------------------------------------------------------------- outliers
    flags = {r["sample"]: [] for r in rows}
    for metric, label in (("missing_rate", "missingness"), ("ts_tv", "Ts/Tv"), ("het_homalt_ratio", "Het/HomAlt"),
                          ("mean_dp", "mean DP"), ("mean_gq", "mean GQ"), ("n_nonref", "variant count"),
                          ("singletons", "singletons"), ("het_rate", "het rate")):
        f, med, mad = mad_flags([fl(r, metric) for r in rows], min_abs=0.02 if metric == "missing_rate" else 0.0)
        for r, bad in zip(rows, f):
            if bad:
                flags[r["sample"]].append(f"{label} outlier")
    for r in rows:
        mr = fl(r, "missing_rate")
        if mr is not None and mr > 0.05:
            flags[r["sample"]].append(f"missingness {mr:.1%} > 5%")

    # ---------------------------------------------------------------- sex check
    declared = read_sex_file(a.sex_file)
    auto_het = {r["sample"]: fl(r, "het_rate") for r in rows}
    sex = []
    yrates = [fl(r, "chrY_call_rate") or 0 for r in rows]
    ymax = max(yrates) if yrates else 0
    for r in rows:
        s = r["sample"]
        xh = fl(r, "chrX_het_rate")
        xc = fl(r, "chrX_nonPAR_called") or 0
        yr = fl(r, "chrY_call_rate")
        ynr = fl(r, "chrY_nonref") or 0
        ah = auto_het.get(s)
        ratio = xh / ah if (xh is not None and ah) else None
        inferred = "unknown"
        if xc >= 20 and ratio is not None:
            if ratio < 0.2:
                inferred = "M"
            elif ratio > 0.4:
                inferred = "F"
            else:
                inferred = "ambiguous"
            if yr is not None and ymax > 0:
                yrel = yr / ymax
                if inferred == "M" and yrel < 0.2 and ynr == 0:
                    inferred = "ambiguous (X0?)"
                elif inferred == "F" and yrel > 0.6 and ynr > 0:
                    inferred = "ambiguous (XXY?)"
        d = declared.get(s, "NA")
        status = "n/a" if d in ("NA", "U") or inferred.startswith(("unknown", "ambig")) else (
            "OK" if d == inferred else "MISMATCH")
        if status == "MISMATCH":
            flags[s].append("sex mismatch")
        sex.append({"sample": s, "chrX_nonPAR_called": int(xc), "chrX_het_rate": xh, "auto_het_rate": ah,
                    "x_auto_het_ratio": ratio, "chrY_call_rate": yr, "chrY_nonref": int(ynr),
                    "declared": d, "inferred": inferred, "status": status})
    out["sex_check"] = sex

    # ---------------------------------------------------------------- contamination
    ab = prof.get("per_sample_ab_hist", {})
    het_flags, med, mad = mad_flags([fl(r, "het_rate") for r in rows], 3.0, min_rel=0.1)
    cont = []
    for r, hf in zip(rows, het_flags):
        s = r["sample"]
        h = ab.get(s) or []
        n = sum(h)
        centers = [(i + 0.5) / 20 for i in range(20)]
        mean_ab = sum(c * x for c, x in zip(centers, h)) / n if n else None
        central = sum(h[8:12]) / n if n else None
        skew = fl(r, "het_ab_outlier_frac")
        hr = fl(r, "het_rate")
        flag = []
        if hf and hr is not None and hr > med:
            flag.append("elevated het rate")
        if skew is not None and skew > 0.25:
            flag.append("skewed allele balance")
        if mean_ab is not None and abs(mean_ab - 0.5) > 0.08:
            flag.append("het AB off-centre")
        if flag:
            flags[s].append("possible contamination: " + ", ".join(flag))
        cont.append({"sample": s, "het_rate": hr, "n_het_with_AD": n, "mean_het_AB": mean_ab,
                     "frac_AB_0.4_0.6": central, "frac_AB_outside_0.25_0.75": skew,
                     "flag": "; ".join(flag) or "OK"})
    out["contamination"] = cont

    # ---------------------------------------------------------------- matrix-based
    npz = os.path.join(a.profile_dir, "gt_matrix.npz")
    G = chrom = pos = None
    if os.path.exists(npz):
        d = np.load(npz, allow_pickle=False)
        G, chrom, pos = d["dosage"], d["chrom"], d["pos"]
        out["gt_matrix_sites"] = int(G.shape[0])
        out["gt_matrix_eligible"] = int(d["n_eligible"])

    if G is not None and G.shape[1] >= 2 and G.shape[0] >= 100:
        kin, ibs0, nboth = king_robust(G)
        pairs = []
        n = len(samples)
        for i in range(n):
            for j in range(i + 1, n):
                k = float(kin[i, j]) if np.isfinite(kin[i, j]) else None
                if k is None:
                    continue
                rel = relation(k, float(ibs0[i, j]))
                if n <= 60 or rel != "unrelated":
                    pairs.append({"sample1": samples[i], "sample2": samples[j], "kinship": k,
                                  "ibs0": float(ibs0[i, j]), "n_sites": int(nboth[i, j]), "relationship": rel})
                if rel.startswith("duplicate"):
                    flags[samples[i]].append(f"duplicate of {samples[j]}")
                    flags[samples[j]].append(f"duplicate of {samples[i]}")
        pairs.sort(key=lambda p: -p["kinship"])
        out["kinship_pairs"] = pairs
        out["kinship_matrix"] = np.nan_to_num(kin).round(4).tolist() if n <= 100 else None
        pc = pca(G)
        if pc:
            out["pca"] = pc

    roh = read_bcftools_roh(a.roh)
    roh_src = "bcftools roh"
    if roh is None and G is not None:
        roh = roh_fallback(G, chrom, pos, samples)
        roh_src = f"fallback scan of {G.shape[0]:,} down-sampled sites (low resolution)"
    if roh is not None:
        summ = []
        for s in samples:
            segs = roh.get(s, [])
            segs = [x for x in segs if x[2] - x[1] >= 1_000_000]
            tot = sum(x[2] - x[1] for x in segs) / 1e6
            longest = max([(x[2] - x[1]) / 1e6 for x in segs], default=0)
            per_chr = {}
            for x in segs:
                per_chr[x[0]] = per_chr.get(x[0], 0) + (x[2] - x[1]) / 1e6
            maxchr = max(per_chr.items(), key=lambda kv: kv[1]) if per_chr else ("", 0)
            note = []
            froh = tot / a.genome_size_mb
            if froh > 0.0625:
                note.append("F_ROH>1/16 (consanguinity ≥ first cousins)")
            if maxchr[1] > 30 and froh < 0.05:
                note.append(f"single large ROH on chr{maxchr[0]} – possible UPD")
            if note:
                flags[s].append("ROH: " + "; ".join(note))
            summ.append({"sample": s, "n_roh_ge1Mb": len(segs), "total_Mb": tot, "longest_Mb": longest,
                         "F_ROH": froh, "max_chrom": maxchr[0], "max_chrom_Mb": maxchr[1],
                         "note": "; ".join(note) or ""})
        out["roh"] = {"source": roh_src, "summary": summ}

    out["flags"] = {s: v for s, v in flags.items()}
    write_json(out, os.path.join(a.out, "sample_qc.json"))
    print(f"[sample_qc] done for {len(samples)} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
