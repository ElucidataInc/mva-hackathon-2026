#!/usr/bin/env python3
"""Figure generation for the VCF QC report (matplotlib, PNG). Every function returns
the written path or None when the input data is absent."""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# validated categorical palette (fixed order) + text tokens
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
ALT_CHROM = ["#2a78d6", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9, "axes.titlesize": 10.5,
    "axes.titleweight": "bold", "axes.titlecolor": INK, "axes.titlelocation": "left", "legend.frameon": False,
    "figure.dpi": 110,
})


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def bar(labels, values, path, title, xlabel="", ylabel="count", horizontal=False, color=SERIES[0], log=False,
        annotate=True, figsize=(7, 3.4)):
    if not labels:
        return None
    fig, ax = plt.subplots(figsize=figsize)
    if horizontal:
        y = np.arange(len(labels))
        ax.barh(y, values, color=color, height=0.7)
        ax.set_yticks(y, labels)
        ax.invert_yaxis()
        ax.set_xlabel(ylabel)
        ax.grid(axis="y", visible=False)
        if log:
            ax.set_xscale("log")
        if annotate:
            for yi, v in zip(y, values):
                ax.text(v, yi, f" {v:,}", va="center", fontsize=7.5, color=INK2)
    else:
        x = np.arange(len(labels))
        ax.bar(x, values, color=color, width=0.75)
        ax.set_xticks(x, labels, rotation=45 if len(labels) > 8 else 0, ha="right" if len(labels) > 8 else "center")
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", visible=False)
        if log:
            ax.set_yscale("log")
        if annotate and len(labels) <= 12:
            for xi, v in zip(x, values):
                ax.text(xi, v, f"{v:,}", ha="center", va="bottom", fontsize=7.5, color=INK2)
    ax.set_title(title)
    return _save(fig, path)


def hist_from_bins(h, path, title, xlabel, ylabel="count", color=SERIES[0], vline=None, logy=False):
    if not h or not h.get("counts") or sum(h["counts"]) == 0:
        return None
    e = np.array(h["edges"])
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.bar(e[:-1], h["counts"], width=np.diff(e), align="edge", color=color, edgecolor=SURF, linewidth=0.4)
    if vline is not None:
        ax.axvline(vline, color=INK2, ls="--", lw=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")
    ax.set_title(title)
    ax.grid(axis="x", visible=False)
    return _save(fig, path)


def multi_hist(hists, path, title, xlabel, max_series=6):
    hists = {k: v for k, v in hists.items() if v and v.get("counts") and sum(v["counts"])}
    if not hists:
        return None
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    for i, (name, h) in enumerate(list(hists.items())[:max_series]):
        e = np.array(h["edges"])
        c = np.array(h["counts"], float)
        c = c / c.sum()
        ax.step(e[:-1], c, where="post", lw=2, color=SERIES[i % 8], label=name)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("fraction")
    ax.set_title(title)
    ax.legend(fontsize=7.5, ncol=2)
    return _save(fig, path)


def stacked_chrom(per_chrom, path, cols, title, normalize=False):
    if not per_chrom:
        return None
    chroms = [r["chrom"] for r in per_chrom]
    data = np.array([[float(r.get(c) or 0) for c in cols] for r in per_chrom])
    if normalize:
        s = data.sum(1, keepdims=True)
        s[s == 0] = 1
        data = data / s
    fig, ax = plt.subplots(figsize=(max(7, len(chroms) * 0.3), 3.6))
    bottom = np.zeros(len(chroms))
    x = np.arange(len(chroms))
    for i, c in enumerate(cols):
        if data[:, i].sum() == 0:
            continue
        ax.bar(x, data[:, i], bottom=bottom, color=SERIES[i % 8], label=c, width=0.8, edgecolor=SURF, linewidth=1)
        bottom += data[:, i]
    ax.set_xticks(x, chroms, rotation=90, fontsize=7.5)
    ax.set_ylabel("fraction" if normalize else "records")
    ax.set_title(title)
    ax.grid(axis="x", visible=False)
    ax.legend(fontsize=7.5, ncol=min(7, len(cols)), loc="upper right", bbox_to_anchor=(1, 1.18))
    return _save(fig, path)


def manhattan(bins, path, value, title, ylabel, chrom_order, log=False):
    """bins: list of dict(chrom, bin_start, value). Alternating-colour Manhattan-style genome map."""
    rows = [b for b in bins if b.get(value) not in (None, "")]
    if not rows:
        return None
    offsets, ticks, labels = {}, [], []
    off = 0
    by_chrom = {}
    for b in rows:
        by_chrom.setdefault(b["chrom"], []).append(b)
    order = [c for c in chrom_order if c in by_chrom]
    fig, ax = plt.subplots(figsize=(11, 3.4))
    for i, c in enumerate(order):
        bs = by_chrom[c]
        xs = np.array([float(b["bin_start"]) for b in bs]) + off
        ys = np.array([float(b[value]) for b in bs])
        ax.scatter(xs, ys, s=6, color=ALT_CHROM[i % 2], linewidths=0)
        end = max(float(b["bin_end"]) for b in bs)
        ticks.append(off + end / 2)
        labels.append(c.replace("chr", ""))
        offsets[c] = off
        off += end
    ax.set_xticks(ticks, labels, fontsize=7, rotation=90 if len(labels) > 30 else 0)
    ax.set_xlim(0, off)
    ax.grid(axis="x", visible=False)
    if log:
        ax.set_yscale("symlog", linthresh=1)
    vals = np.array([float(b[value]) for b in rows])
    ax.axhline(np.median(vals), color=SERIES[1], lw=1, ls="--")
    ax.text(off, np.median(vals), " median", color=SERIES[1], fontsize=7.5, va="center")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return _save(fig, path)


def ideogram_heat(bins, path, value, title, chrom_order, contig_len, bin_size):
    """Chromosome map: one row per chromosome, bins coloured by value (sequential blue)."""
    rows = [b for b in bins if b.get(value) not in (None, "")]
    if not rows:
        return None
    by_chrom = {}
    for b in rows:
        by_chrom.setdefault(b["chrom"], {})[int(b["bin_start"]) // bin_size] = float(b[value])
    order = [c for c in chrom_order if c in by_chrom][:30]
    nb = max(int((contig_len.get(c) or (max(by_chrom[c]) + 1) * bin_size) // bin_size) + 1 for c in order)
    M = np.full((len(order), nb), np.nan)
    for i, c in enumerate(order):
        L = contig_len.get(c)
        if L:
            M[i, : L // bin_size + 1] = 0
        for k, v in by_chrom[c].items():
            if k < nb:
                M[i, k] = v
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    cmap.set_bad(SURF)
    fig, ax = plt.subplots(figsize=(11, 0.28 * len(order) + 1.2))
    vmax = np.nanpercentile(M, 99) if np.isfinite(M).any() else 1
    im = ax.imshow(M, aspect="auto", cmap=cmap, interpolation="nearest", vmin=0, vmax=max(vmax, 1e-9))
    ax.set_yticks(range(len(order)), [c.replace("chr", "") for c in order], fontsize=7)
    xt = np.linspace(0, nb - 1, 6)
    ax.set_xticks(xt, [f"{x * bin_size / 1e6:.0f}" for x in xt])
    ax.set_xlabel("position (Mb)")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cb.ax.tick_params(labelsize=7)
    ax.set_title(title)
    return _save(fig, path)


def scatter(x, y, labels, path, title, xlabel, ylabel, groups=None, annotate=True):
    if not x:
        return None
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    if groups:
        ug = sorted(set(groups))
        for i, g in enumerate(ug):
            idx = [k for k, gg in enumerate(groups) if gg == g]
            ax.scatter([x[k] for k in idx], [y[k] for k in idx], s=40, color=SERIES[i % 8], label=str(g),
                       edgecolors=SURF, linewidths=1.5)
        ax.legend(fontsize=7.5)
    else:
        ax.scatter(x, y, s=40, color=SERIES[0], edgecolors=SURF, linewidths=1.5)
    if annotate and len(x) <= 40:
        for xi, yi, l in zip(x, y, labels):
            ax.annotate(l, (xi, yi), fontsize=6.5, color=INK2, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return _save(fig, path)


def heatmap(M, labels, path, title, vmin=0, vmax=0.5):
    if M is None or len(labels) < 2:
        return None
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    n = len(labels)
    fig, ax = plt.subplots(figsize=(min(10, 2 + n * 0.3), min(9, 1.6 + n * 0.3)))
    im = ax.imshow(np.array(M), cmap=cmap, vmin=vmin, vmax=vmax)
    fs = 7 if n <= 30 else 5
    ax.set_xticks(range(n), labels, rotation=90, fontsize=fs)
    ax.set_yticks(range(n), labels, fontsize=fs)
    ax.grid(False)
    fig.colorbar(im, ax=ax, fraction=0.04, label="KING kinship")
    ax.set_title(title)
    return _save(fig, path)


def per_sample_bars(samples, series, path, title, ylabel, stacked=True):
    """series: dict name -> list values (aligned with samples)."""
    if not samples:
        return None
    fig, ax = plt.subplots(figsize=(max(6, len(samples) * 0.35), 3.6))
    x = np.arange(len(samples))
    bottom = np.zeros(len(samples))
    w = 0.8 if stacked else 0.8 / len(series)
    for i, (name, vals) in enumerate(series.items()):
        v = np.array([float(a or 0) for a in vals])
        if stacked:
            ax.bar(x, v, bottom=bottom, color=SERIES[i % 8], label=name, width=w, edgecolor=SURF, linewidth=1)
            bottom += v
        else:
            ax.bar(x - 0.4 + w * (i + 0.5), v, width=w, color=SERIES[i % 8], label=name)
    ax.set_xticks(x, samples, rotation=90, fontsize=7 if len(samples) < 60 else 4)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="x", visible=False)
    if len(series) > 1:
        ax.legend(fontsize=7.5, ncol=len(series), loc="upper right", bbox_to_anchor=(1, 1.18))
    return _save(fig, path)


def indel_hist(indel_len, path, max_len=30):
    if not indel_len:
        return None
    d = {int(k): v for k, v in indel_len.items()}
    xs = list(range(-max_len, max_len + 1))
    ys = [d.get(x, 0) for x in xs]
    ys[0] += sum(v for k, v in d.items() if k < -max_len)
    ys[-1] += sum(v for k, v in d.items() if k > max_len)
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    cols = [SERIES[1] if x < 0 else SERIES[0] for x in xs]
    ax.bar(xs, ys, color=cols, width=0.8)
    ax.set_yscale("log")
    ax.set_xlabel(f"indel length (bp; deletions < 0 < insertions; ends pooled at ±{max_len})")
    ax.set_ylabel("alleles (log)")
    ax.set_title("Indel length distribution")
    ax.text(0.02, 0.92, "deletions", transform=ax.transAxes, color=SERIES[1], fontsize=8, fontweight="bold")
    ax.text(0.88, 0.92, "insertions", transform=ax.transAxes, color=SERIES[0], fontsize=8, fontweight="bold")
    ax.grid(axis="x", visible=False)
    return _save(fig, path)


def missing_hist(counts, path, title, xlabel):
    if not counts or sum(counts) == 0:
        return None
    xs = np.arange(len(counts)) / 20
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.bar(xs, counts, width=0.045, color=SERIES[0])
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("sites (log)")
    ax.set_title(title)
    ax.grid(axis="x", visible=False)
    return _save(fig, path)


def load_json(p):
    return json.load(open(p)) if os.path.exists(p) else None
