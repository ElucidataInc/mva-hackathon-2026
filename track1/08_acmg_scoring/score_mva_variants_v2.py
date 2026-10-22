#!/usr/bin/env python3

import argparse
import math
import os
import re
import pandas as pd

# ------------------------------------------------------------------ gene panel
G_WEIGHT = {
    "BUB1B": 3, "CEP57": 3, "TRIP13": 3,
    "CENATAC": 2, "MAD2L1BP": 2, "MAD1L1": 2, "BUB1": 2, "SMC5": 2, "SLF2": 2,
    "CEP192": 1, "TRIM37": 1, "BUB3": 1, "CDC20": 1,
}
PANEL = set(G_WEIGHT)
# PVS1 only where loss of function is an established disease mechanism (G >= 2)
PVS1_GENES = {g for g, w in G_WEIGHT.items() if w >= 2}

NULL_TERMS = {"stop_gained", "frameshift_variant", "splice_donor_variant",
              "splice_acceptor_variant"}
PROTEIN_TERMS = {"missense_variant", "start_lost", "stop_lost",
                 "inframe_insertion", "inframe_deletion", "protein_altering_variant"}
SPLICE_REGION_TERMS = {"splice_region_variant", "splice_polypyrimidine_tract_variant",
                       "splice_donor_region_variant", "splice_donor_5th_base_variant"}


def num(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def acmg_class(pts):
    if pts >= 10:
        return "Pathogenic"
    if pts >= 6:
        return "Likely pathogenic"
    if pts >= 0:
        return "VUS"
    if pts >= -6:
        return "Likely benign"
    return "Benign"


# ------------------------------------------------------------------ SNV scoring
def norm_gt(gt):
    g = str(gt).replace("|", "/")
    if g in ("1/1",):
        return "hom"
    if g in ("0/1", "1/0"):
        return "het"
    if g in ("1/2", "2/1"):
        return "het-alt"
    return "unknown"


def score_snv(r, keep_lowqc):
    gene = str(r["marker_gene"])
    sym = str(r["symbol"])
    cons = set(str(r["consequence"]).split("&"))
    af = num(r["MAX_AF"])
    am = num(r["am_path"])
    gt = norm_gt(r["gt"])
    flags = []
    basis = []

    out = dict(gt_called=gt, tier=4, acmg_pts=None, acmg_class="Unscored",
               G=0, total=0, basis="", flags="")

    # G applies only when the variant is in the panel gene itself
    in_gene = (sym == gene)
    G = G_WEIGHT.get(gene, 0) if in_gene else 0
    out["G"] = G
    out["total"] = G
    if not in_gene:
        out["flags"] = f"in neighbouring gene {sym}: gene relevance not applied"
        return out

    qc_ok = (str(r["filter"]) == "PASS") and (str(r.get("in_high_conf", "yes")) != "no")
    if not qc_ok and not keep_lowqc:
        out["flags"] = f"fails QC ({r['filter']}): not scored"
        return out

    is_null = bool(cons & NULL_TERMS)
    is_prot = bool(cons & PROTEIN_TERMS)
    is_splreg = bool(cons & SPLICE_REGION_TERMS)
    if not (is_null or is_prot or is_splreg):
        out["flags"] = "non-coding or synonymous: no scorable evidence"
        return out

    pts = 0
    # PVS1 with NMD logic from the exon field, e.g. '17/23'
    if is_null and gene in PVS1_GENES:
        m = re.match(r"^(\d+)/(\d+)$", str(r.get("exon", "")))
        if m:
            n, tot = int(m.group(1)), int(m.group(2))
            if n == tot:
                pts += 4
                basis.append("PVS1_Moderate +4 (last exon, escapes NMD)")
            elif n == tot - 1:
                pts += 8
                basis.append("PVS1 +8")
                flags.append("penultimate exon: check stop is >50 nt from the exon end")
            else:
                pts += 8
                basis.append(f"PVS1 +8 (exon {n}/{tot}, NMD expected)")
        else:
            pts += 8
            basis.append("PVS1 +8")
            flags.append("exon position unavailable: NMD not checked")
    elif is_null:
        flags.append("null variant in a gene without established LoF mechanism: PVS1 not applied")

    # Population frequency
    if af is None or af < 0.001:
        pts += 1
        basis.append("PM2_Supporting +1")
    elif af < 0.01:
        pts -= 1
        basis.append("BS1_Supporting -1")
    elif af < 0.05:
        pts -= 4
        basis.append("BS1 -4")
    else:
        pts -= 8
        basis.append("BA1 -8")

    # AlphaMissense raw score
    if is_prot and am is not None:
        if am >= 0.990:
            pts += 4; basis.append("PP3 +4 (AM)")
        elif am >= 0.906:
            pts += 2; basis.append("PP3 +2 (AM)")
        elif am >= 0.792:
            pts += 1; basis.append("PP3 +1 (AM)")
        elif am <= 0.070:
            pts -= 2; basis.append("BP4 -2 (AM)")
        elif am <= 0.099:
            pts -= 1; basis.append("BP4 -1 (AM)")

    if is_splreg and not is_null:
        flags.append("splice region: SpliceAI needed, scored 0")

    cv = str(r.get("clin_sig", ".")).lower()
    has_p = "pathogenic" in cv and "likely_benign" not in cv and "benign" not in cv
    has_b = "benign" in cv
    if "pathogenic" in cv and has_b:
        flags.append("ClinVar conflicting/co-located: use allele-matched ClinVar")
    elif has_p:
        flags.append("ClinVar P/LP: PP5 retired, check underlying evidence")
    elif has_b:
        flags.append("ClinVar B/LB: BP6 retired, check evidence")
    if gt in ("hom", "het-alt"):
        flags.append("GT " + gt + ": check for biallelic effect")

    cls = acmg_class(pts)
    if cls in ("Pathogenic", "Likely pathogenic"):
        tier = 1
    elif cls in ("Likely benign", "Benign"):
        tier = 5
    elif is_splreg and not is_prot and not is_null:
        tier = 3
    else:
        tier = 2
    out.update(tier=tier, acmg_pts=pts, acmg_class=cls, total=pts + G,
               basis="; ".join(basis), flags="; ".join(flags))
    return out


def mark_compound_het(df):
    sub = df[(df["tier"].isin([1, 2])) & (df["gt_called"] == "het")]
    for gene, grp in sub.groupby("marker_gene"):
        if len(grp) >= 2:
            for i in grp.index:
                f = df.at[i, "flags"]
                df.at[i, "flags"] = (f + "; " if f else "") + \
                    "compound-het candidate, phase unresolved"
    return df


# ------------------------------------------------------------------- SV scoring
def load_panel_bed(path):
    if not path:
        return {}
    genes = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            gene = p[3]
            if gene in PANEL:
                chrom = p[0] if p[0].startswith("chr") else "chr" + p[0]
                genes[gene] = (chrom, int(p[1]), int(p[2]))
    return genes


def chrom_norm(c):
    c = str(c)
    return c if c.startswith("chr") else "chr" + c


def parse_info(info):
    d = {}
    for kv in str(info).split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def recip_overlap(a, b):
    s = max(a[0], b[0]); e = min(a[1], b[1])
    if e <= s:
        return 0.0
    ov = e - s
    return min(ov / (a[1] - a[0]), ov / (b[1] - b[0]))


def score_sv(r, panel_xy, count_source):
    typ = str(r["sv_type"]).upper()
    chrom = r["chrom"]
    start, end = int(r["start"]), int(r["end"])
    size = end - start
    genes = [g for g in str(r["mane_genes"]).split(",") if g and g != "nan"]
    panel_hit = [g for g in genes if g in PANEL]

    if count_source == "annotsv" and num(r["gene_count"]) is not None:
        n_genes = int(num(r["gene_count"]))
    else:
        n_genes = len(genes)

    s1 = s2 = s3 = s4 = 0.0
    notes = []

    if typ in ("DEL", "DUP") and not genes:
        s1 = -0.60

    # S2: deletion fully containing a panel gene, or inversion breakpoint inside one
    contained = []
    for g in panel_hit:
        if g in panel_xy:
            gc, gs, ge = panel_xy[g]
            if typ == "DEL" and gc == chrom and start <= gs and end >= ge:
                contained.append(g)
            if typ == "INV" and gc == chrom and (gs <= start <= ge or gs <= end <= ge):
                contained.append(g)
        elif typ in ("DEL", "INV"):
            contained.append(g)
            notes.append(f"{g} containment not verified (no panel bed)")
    if contained:
        s2 = 1.00

    if typ == "DEL":
        s3 = 0.90 if n_genes >= 35 else 0.45 if n_genes >= 25 else 0.0
    elif typ == "DUP":
        s3 = 0.90 if n_genes >= 50 else 0.45 if n_genes >= 35 else 0.0

    # S4: distance to nearest panel gene for SVs that do not overlap one
    dist = None
    if not panel_hit and panel_xy:
        best = None
        for g, (gc, gs, ge) in panel_xy.items():
            if gc != chrom:
                continue
            d = max(gs - end, start - ge, 0)
            best = d if best is None else min(best, d)
        dist = best
        if typ in ("DEL", "DUP") and dist is not None:
            s4 = 0.45 if dist <= 100_000 else 0.15 if dist <= 1_000_000 else 0.0
    elif panel_hit:
        dist = 0

    # population frequency from the best gnomAD-SV match
    af = num(r["best_af"]); ofp = num(r["best_ofp"])
    if af is not None and af >= 0.01 and (ofp is None or ofp >= 0.5):
        s1 -= 1.00
        notes.append("common in gnomAD-SV (best match)")

    if typ == "INV":
        s1 = s3 = 0.0
        notes.append("copy-neutral: Riggs CNV sections do not apply")

    ev = round(s1 + s2 + s3 + s4, 2)
    return dict(n_genes=n_genes, panel_genes=",".join(panel_hit), dist_to_panel_gene=dist,
                S1=s1, S2=s2, S3=s3, S4=s4, evidence=ev, notes="; ".join(notes))


def priority(ev, conf):
    if ev >= 0.90 and conf == "Low":
        return "VALIDATE FIRST"
    if ev >= 0.90:
        return "Review as candidate"
    if ev >= 0.15:
        return "Context only"
    return "Low / deprioritise"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snv", required=True)
    ap.add_argument("--sv", required=True)
    ap.add_argument("--panel-bed")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--keep-lowqc", action="store_true",
                    help="score rows that fail the filter / in_high_conf check")
    ap.add_argument("--count-source", choices=["mane", "annotsv"], default="mane",
                    help="gene count for S3: MANE gene list (default) or AnnotSV Gene_count")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    # ---- SNVs
    snv = pd.read_csv(a.snv, sep="\t", dtype=str).fillna(".")
    res = snv.apply(lambda r: score_snv(r, a.keep_lowqc), axis=1, result_type="expand")
    snv = pd.concat([snv, res], axis=1)
    snv = mark_compound_het(snv)
    snv["total"] = pd.to_numeric(snv["total"])
    snv = snv.sort_values(["tier", "total"], ascending=[True, False])
    keep = ["tier", "marker_gene", "chrom", "pos", "ref", "alt", "consequence", "hgvsc",
            "hgvsp", "gt", "gt_called", "filter", "MAX_AF", "am_path", "clin_sig",
            "acmg_pts", "acmg_class", "G", "total", "basis", "flags"]
    snv[keep].to_csv(os.path.join(a.outdir, "MVA_scored_SNVs_v2.tsv"),
                     sep="\t", index=False)

    # ---- SVs
    raw = pd.read_csv(a.sv, sep="\t", dtype=str, low_memory=False)
    if "Annotation_mode" in raw.columns:
        raw = raw[raw["Annotation_mode"] == "full"]
    info = raw["INFO"].apply(parse_info)
    sv = pd.DataFrame({
        "id": raw["ID"],
        "chrom": raw["SV_chrom"].apply(chrom_norm),
        "start": pd.to_numeric(raw["SV_start"]),
        "end": pd.to_numeric(raw["SV_end"]),
        "sv_type": raw["SV_type"],
        "caller": info.apply(lambda d: d.get("SRC", "")).str.lower(),
        "mane_genes": raw["MANE_GENE"] if "MANE_GENE" in raw.columns
        else info.apply(lambda d: d.get("MANE_GENE", "")),
        "gene_count": raw["Gene_count"],
        "best_af": raw["Best_gnomAD_AF"],
        "best_ofp": raw["Best_gnomAD_OFP"],
        "annotsv_score": raw["AnnotSV_ranking_score"],
        "annotsv_class": raw["ACMG_class"],
        "annotsv_criteria": raw["AnnotSV_ranking_criteria"],
    }).reset_index(drop=True)

    panel_xy = load_panel_bed(a.panel_bed)
    if a.panel_bed and not panel_xy:
        print("WARNING: no panel genes found in the bed; check column order")

    res = sv.apply(lambda r: score_sv(r, panel_xy, a.count_source),
                   axis=1, result_type="expand")
    sv = pd.concat([sv, res], axis=1)

    # caller concordance: other caller, same type, >= 50% reciprocal overlap
    conc = []
    for i, r in sv.iterrows():
        other = sv[(sv["chrom"] == r["chrom"]) & (sv["sv_type"] == r["sv_type"]) &
                   (sv["caller"] != r["caller"])]
        hit = any(recip_overlap((r["start"], r["end"]), (o.start, o.end)) >= 0.5
                  for o in other.itertuples())
        conc.append(hit)
    sv["both_callers"] = conc
    size = sv["end"] - sv["start"]
    sv["confidence"] = ["High" if b else ("Low" if (t == "INV" or s >= 1_000_000) else "Medium")
                        for b, t, s in zip(sv["both_callers"], sv["sv_type"], size)]
    sv["size_bp"] = size
    sv["priority"] = [priority(e, c) for e, c in zip(sv["evidence"], sv["confidence"])]
    order = {"VALIDATE FIRST": 0, "Review as candidate": 1,
             "Context only": 2, "Low / deprioritise": 3}
    sv["_o"] = sv["priority"].map(order)
    sv = sv.sort_values(["_o", "evidence"], ascending=[True, False]).drop(columns="_o")

    cols = ["priority", "panel_genes", "caller", "both_callers", "sv_type", "chrom", "start",
            "end", "size_bp", "n_genes", "dist_to_panel_gene", "S1", "S2", "S3", "S4",
            "evidence", "confidence", "best_af", "best_ofp", "annotsv_score",
            "annotsv_class", "annotsv_criteria", "notes", "id"]
    sv[cols].to_csv(os.path.join(a.outdir, "all_SVs_scored_v2.tsv"), sep="\t", index=False)
    focus = sv[(sv["panel_genes"] != "") | (sv["dist_to_panel_gene"].fillna(1e12) <= 1_000_000)]
    focus[cols].to_csv(os.path.join(a.outdir, "MVA_scored_SVs_v2.tsv"), sep="\t", index=False)

    print(f"SNV rows: {len(snv)}  tiers: {snv['tier'].value_counts().sort_index().to_dict()}")
    print(f"SV rows: {len(sv)}  panel-linked or within 1 Mb: {len(focus)}")
    print(sv["priority"].value_counts().to_string())


if __name__ == "__main__":
    main()
