#!/usr/bin/env python3
"""09: FINAL AlphaGenome variant table = score_variant summary + predict_variant local summary
+ call quality + known-hypomorph check + interpretation. Reads saved files only (no API call)."""
import os, glob, numpy as np, pandas as pd
BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/alphagenome/v2_team11"))
AG = f"{BASE}/alphagenome_model"; P = f"{AG}/predict"; HALF = 500
FINAL = f"{AG}/FINAL_alphagenome_variant_table.tsv"; DICT = f"{AG}/FINAL_columns_README.tsv"
S = pd.read_csv(f"{AG}/rare_alphagenome_summary.tsv", sep="\t", low_memory=False)
V = pd.read_csv(f"{BASE}/marker_region_vep.tsv", sep="\t", low_memory=False, dtype=str)

# 1. known BUB1B upstream hypomorph (rs576524605, GRCh38 15:40117088 G>A) anywhere in the patient's calls
h = V[(V.chrom == "chr15") & (V.pos == "40117088")]
hyp = "absent from patient calls" if len(h) == 0 else "PRESENT: " + "; ".join(h.ref + ">" + h.alt + " " + h.zygosity + " " + h["filter"])

# 2. predict_variant: local (+/-500 bp) summary per variant folder (fid)
pv = []
for d in sorted(glob.glob(f"{P}/*/")):
    fid = os.path.basename(d.rstrip("/"))
    if not os.path.exists(d + "diff_summary.tsv"): continue
    rows = []
    for f in glob.glob(d + "tracks_*.npz"):
        mod = os.path.basename(f)[7:-4]; z = np.load(f, allow_pickle=True)
        if "resolution" not in z.files: continue
        ref, alt, res = z["ref"].astype(float), z["alt"].astype(float), int(z["resolution"])
        c = ref.shape[0] // 2; w = max(1, HALF // res); lo, hi = max(0, c - w), c + w + 1
        md = pd.read_csv(d + f"metadata_{mod}.tsv", sep="\t")
        for t in range(ref.shape[1]):
            rows.append((mod, md["biosample_name"].iloc[t] if "biosample_name" in md else "", ref[lo:hi, t].sum(), alt[lo:hi, t].sum()))
    df = pd.DataFrame(rows, columns=["mod", "bio", "ref", "alt"]); spl = df["mod"].str.startswith("splice")
    ns = df[~spl].copy(); ns["l2"] = np.log2((ns.alt + 1) / (ns.ref + 1))
    sig = ns[(ns.ref >= 5) & (ns.l2.abs() >= 0.5)]
    base = ns[ns.ref >= 5] if (ns.ref >= 5).any() else ns
    top = base.loc[base.l2.abs().idxmax()]
    sp = df[spl]; spad = (sp.alt - sp.ref).abs().max() if len(sp) else np.nan
    pv.append({"fid": fid, "pv_done": True, "pv_top_modality": top["mod"], "pv_top_biosample": top.bio,
               "pv_top_local_log2": round(top.l2, 3), "pv_top_local_ref": round(top.ref, 2),
               "pv_n_tracks_sig": len(sig), "pv_n_biosamples_sig": sig.bio.nunique(),
               "pv_sig_modalities": ";".join(sorted(sig["mod"].unique())),
               "pv_sig_direction": ("" if len(sig) == 0 else "up" if (sig.l2 > 0).all() else "down" if (sig.l2 < 0).all() else "mixed"),
               "pv_splice_max_absdiff": round(spad, 4), "pv_plot_dir": f"predict/{fid}/"})
F = S.merge(pd.DataFrame(pv), on="fid", how="left")
F["pv_done"] = F["pv_done"].fillna(False).astype(bool)

# 3. score_variant: strongest regulatory signal (DIFF_LOG2_SUM scorers, already in summary)
reg = ["dnase", "atac", "cage", "procap", "chip_tf", "chip_histone"]
def topreg(r):
    best = max(reg, key=lambda c: abs(r[f"{c}_q"]) if pd.notna(r[f"{c}_q"]) else -1)
    return pd.Series({"sv_top_reg_modality": best, "sv_top_reg_q": r[f"{best}_q"], "sv_top_reg_log2": r[f"{best}_raw"], "sv_top_reg_where": r[f"{best}_where"]})
F = pd.concat([F, F.apply(topreg, axis=1)], axis=1)

# 4. quality, evidence, interpretation, Track 1 candidate
F["quality_ok"] = (F["filter"] == "PASS") & (F["in_high_conf"] == "yes")
cv_plp = F["clinvar_exact_clnsig"].astype(str).str.match(r"^(Pathogenic|Likely_pathogenic)")
def evidence(r):
    nb = r.pv_n_biosamples_sig if pd.notna(r.pv_n_biosamples_sig) else 0
    if r.ag_splicing or r.ag_expression or nb >= 2: s = "strong" if (r.ag_splicing or r.ag_expression) else "moderate"
    elif nb == 1 or r.ag_regulatory: s = "weak"
    else: s = "none"
    return s
F["ag_evidence"] = F.apply(evidence, axis=1)
def interp(r):
    p = []
    if r.impact in ("HIGH", "MODERATE"): p.append("protein-level variant: judge by VEP/ClinVar/AlphaMissense (AlphaGenome models RNA/chromatin, not protein)")
    if not r.quality_ok: p.append("LOW-QUALITY CALL: verify in IGV before use")
    nb = r.pv_n_biosamples_sig if pd.notna(r.pv_n_biosamples_sig) else 0
    if nb >= 1: p.append(f"local {r.pv_sig_modalities} change ({r.pv_sig_direction}) in {int(nb)} cell type(s)")
    if r.ag_splicing: p.append("predicted splicing effect")
    if r.ag_expression: p.append("predicted own-gene expression change")
    if r.ag_regulatory and nb == 0: p.append("track-level regulatory flag only")
    if r.pv_sig_direction == "up": p.append("direction = increase (MVA is loss-of-function, so weak relevance)")
    if r.ag_evidence == "none" and r.impact not in ("HIGH", "MODERATE"): p.append("no predicted effect")
    return "; ".join(p)
F["interpretation"] = F.apply(interp, axis=1)
F["track1_candidate"] = F["quality_ok"] & ((F.impact == "HIGH") | cv_plp | ((F.impact == "MODERATE") & (F.am_class == "likely_pathogenic")))
order = {"strong": 0, "moderate": 1, "weak": 2, "none": 3}
F["_e"] = F.ag_evidence.map(order)
F = F.sort_values(["track1_candidate", "quality_ok", "_e", "rank_key", "ag_max_q"], ascending=[False, False, True, True, False]).drop(columns="_e")
F.insert(0, "final_rank", range(1, len(F) + 1))
F = F.rename(columns={"expr_own_q": "sv_expr_own_q", "expr_own_raw": "sv_expr_own_log2fc", "expr_own_where": "sv_expr_own_where",
                      "splice_site_q": "sv_splice_site_q", "splice_site_raw": "sv_splice_site_delta",
                      "splice_usage_q": "sv_splice_usage_q", "splice_usage_raw": "sv_splice_usage_delta",
                      "splice_junction_q": "sv_splice_junction_q", "splice_junction_raw": "sv_splice_junction_delta"})
cols = ["final_rank", "track1_candidate", "marker_gene", "region_part", "consequence", "impact", "hgvsc", "hgvsp",
        "chrom", "pos", "ref", "alt", "zygosity", "filter", "in_high_conf", "quality_ok",
        "MAX_AF", "clinvar_exact_clnsig", "clinvar_exact_clndn", "am_class", "am_path", "rank_key",
        "sv_expr_own_q", "sv_expr_own_log2fc", "sv_expr_own_where", "sv_splice_site_q", "sv_splice_site_delta",
        "sv_splice_usage_q", "sv_splice_usage_delta", "sv_splice_junction_q", "sv_splice_junction_delta",
        "sv_top_reg_modality", "sv_top_reg_q", "sv_top_reg_log2", "sv_top_reg_where",
        "ag_expression", "ag_splicing", "ag_regulatory", "ag_class",
        "pv_done", "pv_top_modality", "pv_top_biosample", "pv_top_local_log2", "pv_top_local_ref",
        "pv_n_tracks_sig", "pv_n_biosamples_sig", "pv_sig_modalities", "pv_sig_direction", "pv_splice_max_absdiff", "pv_plot_dir",
        "ag_evidence", "interpretation", "final_tier", "fid"]
F[cols].to_csv(FINAL, sep="\t", index=False)

D = [
 ("final_rank","09","Order: Track1 candidates, then good-quality calls, then AlphaGenome evidence strength, then rank_key, then ag_max_q"),
 ("track1_candidate","09","True if good-quality call AND (HIGH impact OR allele-exact ClinVar P/LP OR MODERATE with AlphaMissense likely_pathogenic)"),
 ("marker_gene, region_part","01/02","Team marker gene region (gene +/-50 kb, Ensembl 116); upstream / gene_body / downstream by strand"),
 ("consequence, impact, hgvsc, hgvsp","VEP (step 05)","VEP 116 consequence, impact class and HGVS on the marker gene's MANE transcript"),
 ("chrom, pos, ref, alt, zygosity","VCF","Patient variant (GRCh38, chr-prefixed); zygosity from GT"),
 ("filter, in_high_conf, quality_ok","VCF / high-confidence VCF / 09","Caller filter; presence in high-confidence VCF; quality_ok = PASS and in high-confidence"),
 ("MAX_AF","VEP cache","Highest population frequency (gnomAD exomes/genomes, 1000 Genomes); '.' = novel; all rows < 1%"),
 ("clinvar_exact_clnsig, clinvar_exact_clndn","clinvar.norm.vcf.gz","ClinVar significance and condition for the exact allele (chrom,pos,ref,alt)"),
 ("am_class, am_path","VEP AlphaMissense plugin","Missense pathogenicity class and score 0-1"),
 ("rank_key","04/07","1 HIGH, 2 MODERATE, 3 ClinVar P/LP or AlphaMissense LP, 4 novel, 5 other rare"),
 ("sv_expr_own_q / _log2fc / _where","AlphaGenome score_variant (07)","Strongest predicted change in the marker gene's own expression (GeneMaskLFCScorer): quantile, log2 fold change, tissue"),
 ("sv_splice_site/usage/junction _q / _delta","AlphaGenome score_variant (07)","Strongest predicted splice-site, splice-site-usage and junction change (marker gene if available): quantile and raw change"),
 ("sv_top_reg_modality / _q / _log2 / _where","AlphaGenome score_variant (07, 09)","Strongest local regulatory change (DNase/ATAC/CAGE/PRO-cap/TF/histone, DIFF_LOG2_SUM scorers)"),
 ("ag_expression, ag_splicing, ag_regulatory, ag_class","07 (our thresholds)","Flags: expression |q|>=0.95 & |log2FC|>=0.1; splicing site q>=0.95 & delta>=0.2 OR usage |q|>=0.95 & |delta|>=0.1 OR junction |q|>=0.99 & |delta|>=0.5; regulatory |q|>=0.99 & |log2|>=0.5"),
 ("pv_done","08b","True if predict_variant raw tracks were generated (tier A + B variants)"),
 ("pv_top_modality / _biosample / _local_log2 / _local_ref","AlphaGenome predict_variant (08b, 09)","Strongest LOCAL change within +/-500 bp: log2((ALT+1)/(REF+1)) of window sums, tracks with REF signal >=5 preferred"),
 ("pv_n_tracks_sig, pv_n_biosamples_sig","09","Number of tracks / cell types with |local log2|>=0.5 and REF signal >=5"),
 ("pv_sig_modalities, pv_sig_direction","09","Which modalities changed; up / down / mixed"),
 ("pv_splice_max_absdiff","09","Largest absolute local change in splice tracks (probabilities); >=0.2 would be meaningful"),
 ("pv_plot_dir","08b","Folder with REF vs ALT plots (plot_*.png), cropped arrays (.npz) and track metadata"),
 ("ag_evidence","09","strong = splicing/expression flag; moderate = local change in >=2 cell types; weak = 1 cell type or track-level flag; none"),
 ("interpretation","09","Plain-language reading combining variant type, call quality, AlphaGenome signal and direction"),
 ("final_tier","07","A coding/ClinVar; B non-coding with AlphaGenome signal; C non-coding no signal"),
 ("fid","07","Link to raw score file scores/<fid>.tsv.gz and predict/<fid>/"),
]
pd.DataFrame(D, columns=["column(s)", "made_by", "description"]).to_csv(DICT, sep="\t", index=False)

print("Known BUB1B upstream hypomorph (rs576524605):", hyp)
print(f"rows: {len(F)} -> {FINAL} (server only)\ncolumn guide: {DICT}")
print("\n== Track 1 candidates (by gene) ==")
print(F[F.track1_candidate].groupby("marker_gene").size().to_string() if F.track1_candidate.any() else "none")
print("\n== AlphaGenome evidence x call quality ==")
print(pd.crosstab(F.ag_evidence, F.quality_ok.map({True: "good_call", False: "low_quality"})).to_string())
print("\n== top 10 (no coordinates) ==")
print(F.head(10)[["final_rank", "track1_candidate", "marker_gene", "consequence", "impact", "quality_ok", "ag_evidence", "interpretation"]].to_string(index=False))
