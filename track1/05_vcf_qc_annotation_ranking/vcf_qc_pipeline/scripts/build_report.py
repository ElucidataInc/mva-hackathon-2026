#!/usr/bin/env python3
"""Assemble the self-contained HTML QC report from all module outputs."""

import argparse
import base64
import csv
import datetime as dt
import html
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_plots as P  # noqa: E402
from vcf_common import chrom_sort_key  # noqa: E402

E = html.escape
REPORT_WARNINGS = []


def _safe(fn):
    """Wrap a plotting function so a failure on unusual data skips that figure instead of the report."""
    import functools
    import traceback

    @functools.wraps(fn)
    def w(*args, **kw):
        try:
            return fn(*args, **kw)
        except Exception as e:  # noqa
            msg = f"figure {fn.__name__} skipped: {type(e).__name__}: {e}"
            REPORT_WARNINGS.append(msg)
            print("[report] WARNING " + msg + "\n" + traceback.format_exc(limit=3), file=sys.stderr)
            P.plt.close("all")
            return None
    return w


for _n in ("bar", "hist_from_bins", "multi_hist", "stacked_chrom", "manhattan", "ideogram_heat", "scatter",
           "heatmap", "per_sample_bars", "indel_hist", "missing_hist"):
    setattr(P, _n, _safe(getattr(P, _n)))


# --------------------------------------------------------------------------- helpers
def fnum(v, nd=3):
    if v is None or v == "":
        return "–"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, int):
        return f"{v:,}"
    if isinstance(v, float):
        if v != v:
            return "–"
        if abs(v) >= 1000:
            return f"{v:,.0f}"
        if abs(v) < 0.001 and v != 0:
            return f"{v:.2e}"
        return f"{v:.{nd}g}" if abs(v) < 1 else f"{v:,.{nd}f}".rstrip("0").rstrip(".")
    try:
        x = float(v)
        return fnum(int(x) if x.is_integer() and "." not in str(v) else x, nd)
    except (TypeError, ValueError):
        return E(str(v))


def pct(a, b):
    return f"{100 * a / b:.2f}%" if b else "–"


def table(rows, cols=None, headers=None, max_rows=None, cls=""):
    if not rows:
        return '<p class="na">No data.</p>'
    cols = cols or list(rows[0].keys())
    headers = headers or cols
    body = []
    shown = rows if not max_rows else rows[:max_rows]
    for r in shown:
        tds = []
        for c in cols:
            v = r.get(c) if isinstance(r, dict) else r[cols.index(c)]
            raw = "" if v is None else str(v)
            tds.append(f'<td data-v="{E(raw)}">{v if isinstance(v, str) and v.startswith("<span") else fnum(v)}</td>')
        body.append("<tr>" + "".join(tds) + "</tr>")
    more = f'<p class="note">Showing {len(shown):,} of {len(rows):,} rows.</p>' if max_rows and len(rows) > max_rows else ""
    return (f'<div class="tw"><table class="sortable {cls}"><thead><tr>' +
            "".join(f"<th>{E(str(h))}</th>" for h in headers) + "</tr></thead><tbody>" + "".join(body) +
            "</tbody></table></div>" + more)


def kv_table(pairs):
    return '<div class="tw"><table class="kv">' + "".join(
        f"<tr><th>{E(str(k))}</th><td>{v if isinstance(v, str) and v.startswith('<') else fnum(v)}</td></tr>"
        for k, v in pairs) + "</table></div>"


def badge(status, text=None):
    icon = {"ok": "✓", "warn": "!", "fail": "✕", "na": "–"}[status]
    label = text or {"ok": "OK", "warn": "Warning", "fail": "Fail", "na": "Not run"}[status]
    return f'<span class="badge {status}"><b>{icon}</b> {E(label)}</span>'


class Figs:
    def __init__(self, d):
        self.d = d
        os.makedirs(d, exist_ok=True)

    def p(self, name):
        return os.path.join(self.d, name + ".png")


def img(path, caption="", wide=False):
    if not path or not os.path.exists(path):
        return ""
    b64 = base64.b64encode(open(path, "rb").read()).decode()
    return (f'<figure class="{"wide" if wide else ""}"><img src="data:image/png;base64,{b64}" alt="{E(caption)}">'
            f'<figcaption>{E(caption)}</figcaption></figure>')


TEXT_COLS = {"chrom", "sample", "sample1", "sample2", "module", "function", "status", "note"}


def read_tsv(p):
    if not os.path.exists(p):
        return []
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    for r in rows:
        for k, v in r.items():
            if k in TEXT_COLS:          # chromosome "1" or sample "1001" must stay text
                continue
            try:
                r[k] = int(v) if v.lstrip("-").isdigit() else (float(v) if v not in ("", None) else None)
            except (ValueError, AttributeError):
                pass
    return rows


def details(summary, inner, open_=False):
    return f'<details{" open" if open_ else ""}><summary>{E(summary)}</summary>{inner}</details>'


# --------------------------------------------------------------------------- report
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--title", default="VCF Quality Control Report")
    a = ap.parse_args()
    o = a.outdir
    prof = P.load_json(os.path.join(o, "04_profile", "profile.json"))
    if prof is None:
        sys.exit("profile.json missing – run the profile module first")
    tools = P.load_json(os.path.join(o, "09_report", "tool_outputs.json")) or {}
    sq = P.load_json(os.path.join(o, "05_sample_qc", "sample_qc.json")) or {}
    per_chrom = read_tsv(os.path.join(o, "04_profile", "per_chrom.tsv"))
    per_sample = read_tsv(os.path.join(o, "04_profile", "per_sample.tsv"))
    bins = read_tsv(os.path.join(o, "04_profile", "density_bins.tsv"))
    F = Figs(os.path.join(o, "09_report", "figures"))
    C = prof["counts"]
    nrec = C.get("records", 0)
    samples = prof["samples"]
    chrom_order = sorted({r["chrom"] for r in per_chrom}, key=chrom_sort_key)
    sections = []   # (id, title, html)

    # ---------------------------------------------------------------- KPIs
    ac = prof["allele_class"]
    sc = prof["site_class"]
    n_snv = ac.get("SNV", 0)
    n_indel = sum(ac.get(k, 0) for k in ("insertion", "deletion", "complex_indel"))
    tstv = C.get("ts", 0) / C["tv"] if C.get("tv") else None
    known = prof["known"]
    ext_known = (tools.get("external_annotation") or {}).get("known")
    known_rate = (ext_known["known_any"] / ext_known["total"]) if ext_known and ext_known.get("total") else (
        known.get("known_any", 0) / nrec if nrec else None)
    kpis = [
        ("Records", fnum(nrec)), ("Samples", fnum(len(samples))), ("SNV alleles", fnum(n_snv)),
        ("Indel alleles", fnum(n_indel)), ("Ts/Tv", fnum(tstv)), ("PASS", pct(C.get("pass_records", 0), nrec)),
        ("Multiallelic sites", pct(C.get("multiallelic_sites", 0), nrec)),
        ("Known (rsID/dbSNP/gnomAD)", f"{100 * known_rate:.1f}%" if known_rate is not None else "–"),
        ("File size", f"{prof['input']['size_bytes'] / 1e6:,.1f} MB"), ("Genome build", E(prof["build"].split(" ")[0])),
    ]
    kpi_html = '<div class="kpis">' + "".join(
        f'<div class="kpi"><div class="kv">{v}</div><div class="kl">{E(k)}</div></div>' for k, v in kpis) + "</div>"

    # ---------------------------------------------------------------- 1. provenance
    h = prof["header"]
    inp = prof["input"]
    prov = kv_table([
        ("Path", E(inp["path"])), ("Size", f"{inp['size_bytes']:,} bytes ({inp['size_bytes'] / 1e6:,.2f} MB)"),
        ("Compression", "BGZF/gzip" if inp["compressed"] else "uncompressed"),
        ("Index", ", ".join(inp["index"]) or "none found"), ("VCF version", h["fileformat"] or "–"),
        ("Genome build (inferred)", E(prof["build"])), ("Reference in header", E(h["reference"] or "–")),
        ("Call set type", E(prof["mode"])), ("Detected callers / tools", E(", ".join(prof["callers"]) or "unknown")),
        ("##source lines", E("; ".join(h["sources"]) or "–")),
        ("Header lines", h["n_header_lines"]), ("INFO / FORMAT / FILTER definitions",
                                                f"{h['n_info']} / {h['n_format']} / {h['n_filter_defs']}"),
        ("Contigs in header", h["n_contigs"]), ("Samples", len(samples)),
        ("Tumor / normal (header)", E(f"{h.get('tumor_sample') or '–'} / {h.get('normal_sample') or '–'}")),
        ("Profiler runtime", f"{prof['runtime_sec']} s"),
    ])
    cmds = table(h["commands"], ["key", "value"], ["Header key", "Value (truncated)"]) if h["commands"] else \
        '<p class="na">No command-line / version lines in header.</p>'
    sample_list = ", ".join(E(s) for s in samples[:200]) + (" …" if len(samples) > 200 else "")
    sections.append(("provenance", "File, size & provenance", prov +
                     "<h3>How it was produced: header command lines</h3>" + cmds +
                     details(f"INFO fields ({len(h['info_ids'])})", "<p class='mono'>" + E(", ".join(h["info_ids"])) + "</p>") +
                     details(f"FORMAT fields ({len(h['format_ids'])})", "<p class='mono'>" + E(", ".join(h["format_ids"])) + "</p>") +
                     details(f"Sample names ({len(samples)})", f"<p class='mono'>{sample_list or '(sites-only VCF)'}</p>")))

    # ---------------------------------------------------------------- 2. integrity
    unsorted = C.get("unsorted_positions", 0) + C.get("noncontiguous_chrom_blocks", 0)
    dups = C.get("duplicate_records_exact", 0)
    rows = [
        ("Sort order", badge("ok", "sorted") if not unsorted else badge("fail", f"{unsorted:,} unsorted events"),
         f"{C.get('unsorted_positions', 0):,} positions go backwards; {C.get('noncontiguous_chrom_blocks', 0):,} "
         f"non-contiguous chromosome blocks; {C.get('chrom_out_of_header_order', 0):,} chromosome transitions out of "
         f"header order"),
        ("Exact duplicate records", badge("ok", "none") if not dups else badge("warn", f"{dups:,}"),
         "same CHROM/POS/REF/ALT repeated"),
        ("Records sharing a position", badge("ok", "0") if not C.get("records_sharing_position") else
         badge("warn", fnum(C["records_sharing_position"])), "split multi-allelics or overlapping calls"),
        ("Chromosomes absent from header", badge("ok", "0") if not C.get("chrom_not_in_header") else
         badge("warn", fnum(C["chrom_not_in_header"])), ""),
        ("REF with non-ACGTN characters", badge("ok", "0") if not C.get("ref_non_acgtn") else
         badge("warn", fnum(C["ref_non_acgtn"])), ""),
        ("Malformed lines", badge("ok", "0") if not C.get("malformed_lines") else badge("fail", fnum(C["malformed_lines"])), ""),
    ]
    if prof.get("ref_fasta"):
        chk, mm = C.get("ref_checked", 0), C.get("ref_mismatch", 0)
        st = "ok" if mm == 0 else ("warn" if mm / max(chk, 1) < 0.001 else "fail")
        rows.append(("REF allele concordance (FASTA)", badge(st, f"{100 * (chk - mm) / chk:.4f}% concordant" if chk else "n/a"),
                     f"{chk:,} checked, {mm:,} mismatching, {C.get('ref_contig_missing_in_fasta', 0):,} on contigs "
                     f"missing from FASTA, {C.get('ref_out_of_bounds', 0):,} beyond contig end"))
    else:
        rows.append(("REF allele concordance (FASTA)", badge("na"), "supply --ref to check REF alleles"))
    nc = tools.get("norm_check")
    if nc:
        rows.append(("bcftools norm --check-ref", badge("ok", "0 mismatches") if not nc["ref_mismatch_lines"]
                     else badge("warn", f"{nc['ref_mismatch_lines']:,} REF_MISMATCH"), E(nc.get("norm_summary", ""))))
    vv = tools.get("validate_variants")
    if vv:
        rows.append(("GATK ValidateVariants", badge("ok", "passed") if vv["passed"] else
                     badge("fail", f"{vv['n_error_lines']} error lines"), f"{vv['n_warn_lines']} warnings"))
    else:
        rows.append(("GATK ValidateVariants", badge("na"), "requires gatk + --ref"))
    integ = table([{"Check": r[0], "Result": r[1], "Detail": r[2]} for r in rows])
    ex = prof["examples"]
    for k, lab in (("sort_issues", "Sort-order examples"), ("duplicates", "Duplicate examples"),
                   ("ref_mismatch", "REF mismatch examples")):
        if ex.get(k):
            integ += details(lab, "<pre>" + E("\n".join(ex[k])) + "</pre>")
    if nc and nc.get("examples"):
        integ += details("bcftools norm REF_MISMATCH lines", "<pre>" + E("\n".join(nc["examples"])) + "</pre>")
    if vv and (vv["errors"] or vv["warnings"]):
        integ += details("ValidateVariants messages", "<pre>" + E("\n".join(vv["errors"] + vv["warnings"])) + "</pre>")
    if tools.get("bcftools_index_check"):
        integ += details("Index / contig summary (bcftools index -s)", "<pre>" + E(tools["bcftools_index_check"][:4000]) + "</pre>")
    sections.append(("integrity", "Sort order, duplicates & REF concordance", integ))

    # ---------------------------------------------------------------- 3. composition
    comp = kv_table([
        ("Records (lines)", nrec), ("ALT alleles", sum(ac.values())),
        ("Biallelic sites", f"{C.get('biallelic_sites', 0):,} ({pct(C.get('biallelic_sites', 0), nrec)})"),
        ("Multiallelic sites", f"{C.get('multiallelic_sites', 0):,} ({pct(C.get('multiallelic_sites', 0), nrec)})"),
        ("Sites with no real ALT (ref blocks / '.')", C.get("no_alt_sites", 0)),
        ("Sites with spanning deletion '*'", C.get("sites_with_spanning_deletion", 0)),
        ("Transitions / transversions", f"{C.get('ts', 0):,} / {C.get('tv', 0):,}  (Ts/Tv = {fnum(tstv)})"),
        ("Insertions / deletions", f"{ac.get('insertion', 0):,} / {ac.get('deletion', 0):,}  (ratio "
                                   f"{fnum(ac.get('insertion', 0) / ac['deletion']) if ac.get('deletion') else '–'})"),
        ("Singleton alleles (AC=1)", C.get("singleton_alleles", 0)),
    ])
    sc_sorted = sorted(sc.items(), key=lambda kv: -kv[1])
    ac_sorted = sorted(ac.items(), key=lambda kv: -kv[1])
    f1 = P.bar([k for k, _ in sc_sorted], [v for _, v in sc_sorted], F.p("site_class"),
               "Variant classes (per record)", horizontal=True, ylabel="records", figsize=(6, 3))
    f2 = P.bar([k for k, _ in ac_sorted], [v for _, v in ac_sorted], F.p("allele_class"),
               "Variant kinds (per ALT allele)", horizontal=True, ylabel="alleles", figsize=(6, 3))
    nah = prof["n_alt_hist"]
    f3 = P.bar([f"{k} ALT" for k in nah], list(nah.values()), F.p("n_alt"), "Biallelic vs multiallelic sites",
               ylabel="records", log=True, figsize=(5, 3))
    subst = prof["substitutions"]
    order6 = [("C>A", "G>T"), ("C>G", "G>C"), ("C>T", "G>A"), ("T>A", "A>T"), ("T>C", "A>G"), ("T>G", "A>C")]
    f4 = P.bar([a_ for a_, _ in order6], [subst.get(a_, 0) + subst.get(b_, 0) for a_, b_ in order6],
               F.p("spectrum"), "SNV substitution spectrum (strand-collapsed)", ylabel="SNVs", figsize=(5, 3))
    f5 = P.indel_hist(prof["indel_length"], F.p("indel_len"))
    il = {int(k): v for k, v in prof["indel_length"].items()}
    ind_tab = ""
    if il:
        tot = sum(il.values())
        frame = sum(v for k, v in il.items() if k % 3 == 0)
        cum, med_len = 0, None
        for k_, v_ in sorted(((abs(k), v) for k, v in il.items())):
            cum += v_
            if med_len is None and cum >= tot / 2:
                med_len = k_
        ind_tab = kv_table([("Indel alleles", tot), ("Median |length| (bp)", med_len),
            ("1-bp indels", f"{il.get(1, 0) + il.get(-1, 0):,} ({pct(il.get(1, 0) + il.get(-1, 0), tot)})"),
            ("Length multiple of 3 (in-frame if coding)", f"{frame:,} ({pct(frame, tot)})"),
            ("Longest insertion / deletion", f"{max(il)} / {min(il)} bp")])
    sections.append(("composition", "Variant kinds, classes & indel lengths", comp +
                     '<div class="grid">' + img(f1, "Per-record class") + img(f2, "Per-allele kind") + img(f3, "Sites by number of ALT alleles") +
                     img(f4, "Six substitution classes") + "</div>" +
                     table([{"class": k, "records": v, "fraction": pct(v, nrec)} for k, v in sc_sorted]) +
                     "<h3>Indel length distribution</h3>" + img(f5, "Indel lengths", wide=True) + ind_tab))

    # ---------------------------------------------------------------- 4. chromosomes
    cls_cols = ["SNV", "MNV", "insertion", "deletion", "complex_indel", "indel_mixed", "mixed", "other"]
    cls_cols = [c for c in cls_cols if any((r.get(c) or 0) for r in per_chrom)]
    f6 = P.stacked_chrom(per_chrom, F.p("chrom_types"), cls_cols, "Variant types per chromosome")
    f7 = P.stacked_chrom(per_chrom, F.p("chrom_types_frac"), cls_cols, "Variant-type composition per chromosome",
                         normalize=True)
    dens = [r for r in per_chrom if r.get("density_per_Mb") not in (None, "")]
    f8 = P.bar([r["chrom"] for r in dens], [round(r["density_per_Mb"], 1) for r in dens], F.p("chrom_density"),
               "Variants per Mb by chromosome", ylabel="records / Mb", annotate=False,
               figsize=(max(7, len(dens) * .3), 3.2)) if dens else None
    bs = prof["bin_size"]
    f9 = P.manhattan(bins, F.p("manhattan_density"), "n_variants",
                     f"Genome-wide variant density ({bs / 1e6:g} Mb bins)", "variants / bin", chrom_order)
    f10 = P.ideogram_heat(bins, F.p("ideogram_density"), "n_variants", "Chromosome map: variant density",
                          chrom_order, prof["contig_lengths"], bs)
    f11 = P.manhattan(bins, F.p("manhattan_depth"), "mean_site_dp",
                      f"Genome-wide mean site depth ({bs / 1e6:g} Mb bins)", "mean DP", chrom_order)
    f12 = P.ideogram_heat(bins, F.p("ideogram_depth"), "mean_site_dp", "Chromosome map: mean site depth",
                          chrom_order, prof["contig_lengths"], bs)
    chrom_cols = ["chrom", "length", "records", "pass", "density_per_Mb"] + [c for c in cls_cols] + \
                 ["het_gt", "homalt_gt", "het_homalt_ratio", "missing_rate"]
    sections.append(("chromosomes", "Chromosomes: counts, types, density & depth maps",
                     "<h3>Number of variants and types per chromosome</h3>" +
                     table(per_chrom, chrom_cols) + img(f6, "Records per chromosome by type", wide=True) +
                     img(f7, "Type fractions per chromosome", wide=True) + img(f8, "Density normalised by contig length", wide=True) +
                     "<h3>Variant density maps</h3>" + img(f9, "Manhattan-style density; dashed = median bin", wide=True) +
                     img(f10, "Chromosome heat map of variant counts", wide=True) +
                     "<h3>Depth maps</h3>" + (img(f11, "Manhattan-style mean INFO/DP (or summed FORMAT/DP) per bin", wide=True) +
                                             img(f12, "Chromosome heat map of mean depth", wide=True) if f11 else
                                             '<p class="na">No DP information in this VCF.</p>')))

    # ---------------------------------------------------------------- 5. filters & quality
    flt = prof["filters"]
    fdefs = h.get("filter_defs", {})
    f13 = P.bar([k for k, _ in flt[:20]], [v for _, v in flt[:20]], F.p("filters"), "FILTER values",
                horizontal=True, ylabel="records", log=len(flt) > 1, figsize=(6.5, max(2.5, .3 * min(20, len(flt)))))
    ftab = table([{"FILTER": k, "records": v, "fraction": pct(v, nrec), "SNV": dict(prof["filter_by_type"].get("SNV", [])).get(k, 0),
                   "INDEL": dict(prof["filter_by_type"].get("INDEL", [])).get(k, 0), "description": fdefs.get(k, "")}
                  for k, v in flt])
    q = prof.get("qual") or {}
    f14 = P.hist_from_bins(prof.get("qual_hist"), F.p("qual"), "QUAL distribution", "QUAL (upper tail pooled)")
    f15 = P.hist_from_bins(prof.get("site_dp_hist"), F.p("site_dp"), "Site depth distribution", "DP (upper tail pooled)")
    gi_figs = ""
    gi_rows = []
    for k, d in (prof.get("gatk_info") or {}).items():
        gi_figs += img(P.hist_from_bins(d["hist"], F.p("info_" + k), k, k, color=P.SERIES[2]), k)
        gi_rows.append({"field": k, **{x: d["summary"].get(x) for x in ("n", "mean", "p05", "median", "p95")}})
    sections.append(("filters", "Filters & call quality", img(f13, "FILTER column") + ftab +
                     '<div class="grid">' + img(f14, "QUAL") + img(f15, "Site DP") + "</div>" +
                     kv_table([("QUAL median (p05–p95)", f"{fnum(q.get('median'))} ({fnum(q.get('p05'))}–{fnum(q.get('p95'))})" if q else "–"),
                               ("Site DP median (p05–p95)", f"{fnum((prof.get('site_dp') or {}).get('median'))} "
                                                           f"({fnum((prof.get('site_dp') or {}).get('p05'))}–{fnum((prof.get('site_dp') or {}).get('p95'))})")]) +
                     (("<h3>GATK site annotations (hard-filter metrics)</h3>" + table(gi_rows) + '<div class="grid">' + gi_figs + "</div>")
                      if gi_rows else "")))

    # ---------------------------------------------------------------- 6. genotypes, missingness, phasing
    if samples:
        S = [r["sample"] for r in per_sample]
        f16 = P.per_sample_bars(S, {"het": [r["het"] for r in per_sample], "hom-alt": [r["hom_alt"] for r in per_sample]},
                                F.p("gt_counts"), "Non-reference genotypes per sample", "genotypes")
        f16b = P.per_sample_bars(S, {"hom-ref": [r["hom_ref"] for r in per_sample], "het": [r["het"] for r in per_sample],
                                     "hom-alt": [r["hom_alt"] for r in per_sample], "missing": [r["missing"] for r in per_sample]},
                                 F.p("gt_all"), "All genotype calls per sample", "genotypes")
        tot = {k: sum(r[k] or 0 for r in per_sample) for k in ("hom_ref", "het", "hom_alt", "het_alt_alt", "missing")}
        gtab = kv_table([("Hom-REF (0/0)", tot["hom_ref"]), ("Heterozygous REF/ALT", tot["het"] - tot["het_alt_alt"]),
                         ("Heterozygous ALT/ALT (1/2)", tot["het_alt_alt"]), ("Hom-ALT", tot["hom_alt"]),
                         ("Missing (./.)", tot["missing"]), ("Overall Het/HomAlt", fnum(tot["het"] / tot["hom_alt"]) if tot["hom_alt"] else "–")])
        f17 = P.per_sample_bars(S, {"missing rate": [r["missing_rate"] for r in per_sample]}, F.p("miss_sample"),
                                "Genotype missingness per sample", "fraction missing")
        f18 = P.missing_hist(prof["site_missing_hist"], F.p("miss_site"), "Per-site genotype missingness",
                             "fraction of samples missing")
        vt = tools.get("vcftools") or {}
        f18b = P.missing_hist(vt.get("lmiss_hist"), F.p("miss_site_vcftools"), "Per-site missingness (vcftools --missing-site)",
                              "F_MISS") if vt.get("lmiss_hist") else None
        mc = [r for r in per_chrom if r.get("missing_rate") not in (None, "")]
        f18c = P.bar([r["chrom"] for r in mc], [r["missing_rate"] for r in mc], F.p("miss_chrom"),
                     "Genotype missingness per chromosome", ylabel="fraction missing", annotate=False,
                     figsize=(max(7, len(mc) * .3), 3)) if mc else None
        ph = [{"sample": r["sample"], "het": r["het"], "phased_het": r["phased_het"], "phased_het_frac": r["phased_het_frac"],
               "phased_gt": r["phased_gt"], "phase_sets(PS)": r["phase_sets"]} for r in per_sample]
        any_phased = any((r["phased_gt"] or 0) for r in per_sample)
        f19 = P.per_sample_bars(S, {"phased het fraction": [r["phased_het_frac"] for r in per_sample]}, F.p("phasing"),
                                "Fraction of heterozygous calls phased", "fraction") if any_phased else None
        sections.append(("genotypes", "Genotypes, missingness & phasing",
                         "<h3>Homozygous / heterozygous REF/ALT counts</h3>" + gtab + img(f16, "", wide=True) + img(f16b, "", wide=True) +
                         "<h3>Genotype missingness</h3>" + '<div class="grid">' + img(f17, "Per sample") + img(f18, "Per site (profiler)") +
                         img(f18b, "Per site (vcftools)") + "</div>" + img(f18c, "Per chromosome", wide=True) +
                         "<h3>Phasing</h3>" + (img(f19, "Phased hets ('|' genotypes)", wide=True) + table(ph) if any_phased else
                                               '<p class="na">No phased genotypes (no "|" separators) found.</p>')))

    # ---------------------------------------------------------------- 7. known/novel & annotation
    kn = prof["known"]
    krows = [{"evidence": "ID column not '.'", "records": kn.get("with_ID", 0), "fraction": pct(kn.get("with_ID", 0), nrec)},
             {"evidence": "rsID in ID column", "records": kn.get("with_rsID", 0), "fraction": pct(kn.get("with_rsID", 0), nrec)},
             {"evidence": "INFO/DB flag (dbSNP)", "records": kn.get("DB_flag", 0), "fraction": pct(kn.get("DB_flag", 0), nrec)},
             {"evidence": "gnomAD AF present", "records": kn.get("in_gnomAD", 0), "fraction": pct(kn.get("in_gnomAD", 0), nrec)},
             {"evidence": "Known by any evidence", "records": kn.get("known_any", 0), "fraction": pct(kn.get("known_any", 0), nrec)}]
    tvk = lambda p: kn.get(p + "_ts", 0) / kn[p + "_tv"] if kn.get(p + "_tv") else None
    ktype = [{"type": t, "total": kn.get(f"{t}_total", 0), "known": kn.get(f"{t}_known", 0),
              "known_rate": pct(kn.get(f"{t}_known", 0), kn.get(f"{t}_total", 0))} for t in ("SNV", "INDEL", "OTHER") if kn.get(f"{t}_total")]
    ann_html = "<h3>Known vs novel (from the input VCF)</h3>" + table(krows) + table(ktype) + kv_table(
        [("Ts/Tv known SNVs", fnum(tvk("known"))), ("Ts/Tv novel SNVs", fnum(tvk("novel")))])
    if ext_known:
        ann_html += "<h3>Known vs novel after annotation against dbSNP / gnomAD</h3>" + table(
            [{"evidence": k, "records": v, "fraction": pct(v, ext_known["total"])} for k, v in ext_known.items()
             if not k.endswith(("_total", "_known")) and k != "total"] +
            [{"evidence": f"{t} known", "records": ext_known.get(f"{t}_known", 0),
              "fraction": pct(ext_known.get(f"{t}_known", 0), ext_known.get(f"{t}_total", 0))}
             for t in ("SNV", "INDEL") if ext_known.get(f"{t}_total")])
    for src, ann in (("input VCF", prof["annotation"]), ("external annotation (VEP/SnpEff)",
                                                          (tools.get("external_annotation") or {}).get("annotation"))):
        if not ann or not ann.get("n_annotated_records"):
            continue
        tag = "in" if src == "input VCF" else "ext"
        ann_html += f"<h3>Functional annotation – {E(src)}</h3>"
        ann_html += kv_table([("Annotation sources", E(", ".join(ann["sources"].keys()))),
                              ("Annotated records", ann["n_annotated_records"])])
        g = ""
        if ann["consequence"]:
            g += img(P.bar([k for k, _ in ann["consequence"][:20]], [v for _, v in ann["consequence"][:20]],
                           F.p(f"{tag}_csq"), "Consequences (records)", horizontal=True, log=True,
                           figsize=(6.5, max(2.5, .28 * min(20, len(ann['consequence'])))), color=P.SERIES[0]), "Consequence terms")
        if ann["impact"]:
            ks = [k for k in ("HIGH", "MODERATE", "LOW", "MODIFIER", "NA") if k in ann["impact"]]
            g += img(P.bar(ks, [ann["impact"][k] for k in ks], F.p(f"{tag}_impact"), "Most severe impact per record",
                           figsize=(5, 3), color=P.SERIES[1]), "Impact")
        if ann["gnomad_af_summary"]:
            b_ = ann["gnomad_af_bins"]
            g += img(P.bar([x[0] for x in b_], [x[1] for x in b_], F.p(f"{tag}_gnomad"), "gnomAD allele-frequency bins",
                           figsize=(6, 3), color=P.SERIES[2]), "gnomAD AF")
        if ann["clinsig"]:
            g += img(P.bar([k for k, _ in ann["clinsig"][:12]], [v for _, v in ann["clinsig"][:12]], F.p(f"{tag}_clin"),
                           "ClinVar significance", horizontal=True, figsize=(6, 3), color=P.SERIES[6]), "ClinVar")
        ann_html += '<div class="grid">' + g + "</div>"
        ann_html += details("Top genes / biotypes", table([{"gene": k, "records": v} for k, v in ann["top_genes"]]) +
                            table([{"biotype": k, "records": v} for k, v in ann["biotype"]]))
    if not prof["annotation"].get("n_annotated_records") and not tools.get("external_annotation"):
        ann_html += '<p class="na">No functional annotation in the VCF and the annotation module did not run ' \
                    '(provide --vep-cache, --snpeff-db, --dbsnp or --gnomad).</p>'
    sections.append(("annotation", "Known vs novel & functional annotation", ann_html))

    # ---------------------------------------------------------------- 8. sample QC
    if samples:
        flags = sq.get("flags", {})
        cols = ["sample", "n_nonref", "n_snv", "n_indel", "het", "hom_alt", "het_homalt_ratio", "ts_tv", "missing_rate",
                "mean_dp", "mean_gq", "singletons", "private"]
        rows = []
        for r in per_sample:
            rr = {c: r.get(c) for c in cols}
            fl = flags.get(r["sample"], [])
            rr["flags"] = (f'<span class="badge warn"><b>!</b> {E("; ".join(fl))}</span>' if fl else
                           '<span class="badge ok"><b>✓</b> OK</span>')
            rows.append(rr)
        S = [r["sample"] for r in per_sample]
        s_html = kv_table([("Number of samples", len(samples)), ("Samples with ≥1 flag", sum(1 for v in flags.values() if v))])
        s_html += table(rows, cols + ["flags"], ["sample", "variants", "SNVs", "indels", "het", "hom-alt", "Het/HomAlt",
                                                  "Ts/Tv", "missing", "mean DP", "mean GQ", "singletons", "private", "flags"],
                        max_rows=500)
        s_html += '<div class="grid">' + img(P.per_sample_bars(S, {"Ts/Tv": [r["ts_tv"] for r in per_sample]}, F.p("s_tstv"),
                                                              "Ts/Tv per sample", "Ts/Tv"), "Ts/Tv") + \
            img(P.per_sample_bars(S, {"Het/HomAlt": [r["het_homalt_ratio"] for r in per_sample]}, F.p("s_hethom"),
                                  "Het/HomAlt per sample", "ratio"), "Het/HomAlt") + \
            img(P.per_sample_bars(S, {"singletons": [r["singletons"] for r in per_sample],
                                      "private (non-singleton)": [max(0, (r["private"] or 0) - (r["singletons"] or 0)) for r in per_sample]},
                                  F.p("s_single"), "Singletons and private variants", "alleles"), "Singleton / private") + \
            img(P.per_sample_bars(S, {"mean DP": [r["mean_dp"] for r in per_sample]}, F.p("s_dp"), "Mean DP per sample", "DP"), "Depth") + \
            img(P.per_sample_bars(S, {"mean GQ": [r["mean_gq"] for r in per_sample]}, F.p("s_gq"), "Mean GQ per sample", "GQ"), "GQ") + \
            img(P.multi_hist(prof.get("per_sample_dp_hist", {}), F.p("s_dp_hist"), "Per-sample FORMAT/DP distribution", "DP"),
                "DP distributions (first 6 samples)") + "</div>"

        # sex check
        sx = sq.get("sex_check") or []
        s_html += "<h3>Sex check</h3>"
        if sx:
            s_html += table(sx, ["sample", "chrX_nonPAR_called", "chrX_het_rate", "x_auto_het_ratio", "chrY_call_rate",
                                 "chrY_nonref", "declared", "inferred", "status"])
            s_html += img(P.scatter([r["x_auto_het_ratio"] or 0 for r in sx], [r["chrY_call_rate"] or 0 for r in sx],
                                    [r["sample"] for r in sx], F.p("sexcheck"), "chrX heterozygosity vs chrY call rate",
                                    "chrX/autosome het-rate ratio", "chrY call rate",
                                    groups=[r["declared"] for r in sx]), "Colour = declared sex; males expected bottom-right of X axis ≈ 0")
        for k, lab in (("guess_ploidy", "bcftools +guess-ploidy"), ("peddy_sex_check", "peddy sex check")):
            if tools.get(k):
                s_html += f"<h4>{lab}</h4>" + table(tools[k], max_rows=200)
        if tools.get("check_ploidy"):
            s_html += details("bcftools +check-ploidy", "<pre>" + E(tools["check_ploidy"]) + "</pre>")

        # relatedness
        s_html += "<h3>Relatedness & duplicates</h3>"
        kp = sq.get("kinship_pairs") or []
        if kp:
            s_html += f'<p class="note">KING-robust kinship on {sq.get("gt_matrix_sites", 0):,} down-sampled autosomal PASS biallelic SNVs. ' \
                      'Thresholds: &gt;0.354 duplicate/MZ, 0.177–0.354 first degree, 0.088–0.177 second, 0.044–0.088 third.</p>'
            rel_counts = {}
            for p_ in kp:
                rel_counts[p_["relationship"]] = rel_counts.get(p_["relationship"], 0) + 1
            s_html += kv_table(list(rel_counts.items()))
            kcols = ["sample1", "sample2", "kinship", "ibs0", "n_sites", "relationship"]
            related = [p_ for p_ in kp if p_["relationship"] != "unrelated"]
            s_html += table(related, kcols, max_rows=200) if related else '<p class="note">No related pairs (kinship &gt; 0.044).</p>'
            if len(related) < len(kp):
                s_html += details(f"All {len(kp):,} pairs", table(kp, kcols, max_rows=1000))
            s_html += '<div class="grid">' + img(P.scatter([p_["ibs0"] for p_ in kp], [p_["kinship"] for p_ in kp], [""] * len(kp),
                                                           F.p("kin_scatter"), "Kinship vs IBS0 (all pairs)", "IBS0", "KING kinship",
                                                           groups=[p_["relationship"] for p_ in kp], annotate=False), "Pairs") + \
                img(P.heatmap(sq.get("kinship_matrix"), samples, F.p("kin_heat"), "Kinship matrix"), "Kinship heat map") + "</div>"
        else:
            s_html += '<p class="na">Needs ≥2 samples and ≥100 eligible sites.</p>'
        if tools.get("somalier_pairs"):
            s_html += "<h4>somalier relate</h4>" + table(tools["somalier_pairs"], max_rows=100)
        if tools.get("peddy_ped_check"):
            s_html += "<h4>peddy ped check</h4>" + table(tools["peddy_ped_check"], max_rows=100)

        # ancestry
        s_html += "<h3>Ancestry</h3>"
        pc = sq.get("pca")
        anc = tools.get("somalier_ancestry")
        if anc:
            pc_cols = [c for c in anc[0] if c.startswith("PC")][:2]
            grp = [r.get("predicted_ancestry", "?") for r in anc]
            if len(pc_cols) == 2:
                s_html += img(P.scatter([r[pc_cols[0]] for r in anc], [r[pc_cols[1]] for r in anc],
                                        [r.get("sample_id", "") for r in anc], F.p("somalier_pca"),
                                        "somalier: projection onto 1000G reference", pc_cols[0], pc_cols[1], groups=grp,
                                        annotate=False), "1000G reference + query samples")
            s_html += table(anc, max_rows=200)
        if pc:
            pcs = pc["pcs"]
            ve = pc["var_explained"]
            s_html += img(P.scatter([p_[0] for p_ in pcs], [p_[1] if len(p_) > 1 else 0 for p_ in pcs], samples, F.p("pca"),
                                    "Within-cohort PCA", f"PC1 ({100 * ve[0]:.1f}%)", f"PC2 ({100 * ve[1]:.1f}%)" if len(ve) > 1 else "PC2"),
                          f"{pc['n_sites']:,} sites; reference-population projection needs somalier (--somalier-sites/--somalier-labels)")
        if not pc and not anc:
            s_html += '<p class="na">PCA needs ≥3 samples; reference projection needs somalier with 1000G labelled sites.</p>'

        # contamination
        s_html += "<h3>Contamination screen</h3>"
        cont = sq.get("contamination") or []
        s_html += '<p class="note">VCF-level proxies only (elevated heterozygosity, het allele balance away from 0.5); ' \
                  'confirm with VerifyBamID2 on the BAM/CRAM.</p>'
        s_html += table(cont, max_rows=300)
        abh = prof.get("per_sample_ab_hist", {})
        abh_bins = {k: {"edges": [i / 20 for i in range(21)], "counts": v} for k, v in abh.items()}
        s_html += img(P.multi_hist(abh_bins, F.p("ab_hist"), "Allele balance of heterozygous calls (AD)", "ALT read fraction"),
                      "First 6 samples; a clean diploid sample peaks at 0.5")
        if tools.get("peddy_het_check"):
            s_html += "<h4>peddy het check</h4>" + table(tools["peddy_het_check"], max_rows=100)

        # ROH
        roh = sq.get("roh")
        s_html += "<h3>Runs of homozygosity</h3>"
        if roh:
            s_html += f'<p class="note">Source: {E(roh["source"])}. F_ROH = total ROH ≥1 Mb / autosomal length.</p>'
            s_html += table(roh["summary"], max_rows=300)
            s_html += img(P.per_sample_bars([r["sample"] for r in roh["summary"]], {"F_ROH": [r["F_ROH"] for r in roh["summary"]]},
                                            F.p("froh"), "Genomic inbreeding (F_ROH)", "F_ROH"), "", wide=True)
        else:
            s_html += '<p class="na">Not available.</p>'

        # Mendelian
        s_html += "<h3>Mendelian errors</h3>"
        md = prof.get("mendelian") or []
        if md:
            s_html += table(md)
        if tools.get("mendelian2"):
            s_html += details("bcftools +mendelian2 output", "<pre>" + E(tools["mendelian2"]) + "</pre>", open_=not md)
        if not md and not tools.get("mendelian2"):
            s_html += '<p class="na">Provide --ped with trios whose members are in the VCF.</p>'
        sections.append(("samples", "Sample QC", s_html))

    # ---------------------------------------------------------------- 9. somatic
    som = prof.get("somatic_info") or {}
    if prof["mode"] == "somatic" or som:
        t, n = h.get("tumor_sample"), h.get("normal_sample")
        s_html = kv_table([("Tumor sample", E(t or "not declared")), ("Normal sample", E(n or "not declared"))])
        s_html += '<div class="grid">' + img(P.multi_hist(prof.get("per_sample_vaf_hist", {}), F.p("vaf_samples"),
                                                          "VAF distribution per sample", "VAF (FORMAT/AF or AD)"), "VAF") + \
            img(P.hist_from_bins(prof.get("vaf_hist"), F.p("vaf_all"), "VAF (all non-ref calls)", "VAF"), "VAF pooled")
        for k in ("TLOD", "NLOD", "NALOD", "GERMQ", "POPAF", "ROQ", "STRANDQ", "SEQQ", "CONTQ"):
            if k in som:
                s_html += img(P.hist_from_bins(som[k]["hist"], F.p("som_" + k), k, k, color=P.SERIES[1]), k)
        s_html += img(P.multi_hist(prof.get("per_sample_dp_hist", {}), F.p("som_dp"), "Tumor and normal depth", "FORMAT/DP"),
                      "DP per sample")
        s_html += img(P.multi_hist(prof.get("orientation_bias_hist", {}), F.p("som_f1r2"),
                                   "Orientation bias: ALT F1R2 / (F1R2+F2R1)", "fraction"), "Expected ≈0.5; skew → OxoG/FFPE artefacts")
        s_html += img(P.multi_hist(prof.get("strand_bias_hist", {}), F.p("som_sb"), "Strand bias: ALT forward fraction (SB)",
                                   "fraction"), "Expected ≈0.5") + "</div>"
        s_html += table([{"metric": k, **{x: v["summary"].get(x) for x in ("n", "mean", "p05", "median", "p95")}} for k, v in som.items()])
        pon = [{"FILTER": k, "records": v, "fraction": pct(v, nrec), "description": fdefs.get(k, "")} for k, v in flt
               if any(x in k.lower() for x in ("panel", "pon", "germline", "normal_artifact", "strand", "orientation",
                                               "contamination", "weak_evidence", "clustered", "multiallelic", "slippage",
                                               "haplotype", "fragment", "base_qual", "map_qual", "position"))]
        s_html += "<h4>Panel-of-normals, germline-risk and artefact filters</h4>" + table(pon)
        sections.append(("somatic", "Somatic metrics (VAF, TLOD/NLOD, depth, bias, PON/germline)", s_html))

    # ---------------------------------------------------------------- 10. bcftools stats
    bst = tools.get("bcftools_stats")
    if bst:
        b_html = '<p class="note">From <code>bcftools stats -s - [-F ref]</code>; raw file and plot-vcfstats output are in ' \
                 '<code>03_bcftools_stats/</code>.</p>'
        b_html += "<h3>Summary numbers (SN)</h3>" + table([{"metric": r.get("key"), "value": r.get("value")} for r in bst.get("SN", [])])
        if bst.get("TSTV"):
            b_html += "<h3>Ts/Tv (TSTV)</h3>" + table(bst["TSTV"], [c for c in bst["TSTV"][0] if c not in ("TSTV", "id")])
        if bst.get("PSC"):
            b_html += "<h3>Per-sample counts (PSC)</h3>" + table(bst["PSC"], [c for c in bst["PSC"][0] if c not in ("PSC", "id")], max_rows=500)
        if bst.get("PSI"):
            b_html += details("Per-sample indels (PSI)", table(bst["PSI"], [c for c in bst["PSI"][0] if c not in ("PSI", "id")], max_rows=500))
        g = ""
        if bst.get("AF"):
            g += img(P.bar([fnum(r["allele frequency"]) for r in bst["AF"]], [r.get("number of SNPs", 0) for r in bst["AF"]],
                           F.p("bst_af"), "Non-reference allele frequency spectrum (SNPs)", annotate=False, log=True), "AF spectrum")
        if bst.get("IDD"):
            idd = {str(r["length (deletions negative)"]): r.get("number of sites", r.get("count", 0)) for r in bst["IDD"]}
            g += img(P.indel_hist(idd, F.p("bst_idd")), "IDD")
        if bst.get("QUAL"):
            g += img(P.bar([str(r["Quality"]) for r in bst["QUAL"]][:80], [r.get("number of SNPs", 0) for r in bst["QUAL"]][:80],
                           F.p("bst_qual"), "SNPs by QUAL", annotate=False), "QUAL")
        if bst.get("DP"):
            g += img(P.bar([str(r["bin"]) for r in bst["DP"]][:100], [r.get("number of sites", r.get("number of genotypes", 0)) for r in bst["DP"]][:100],
                           F.p("bst_dp"), "Sites by depth", annotate=False), "DP")
        b_html += '<div class="grid">' + g + "</div>"
        if tools.get("bcftools_stats_pass"):
            b_html += details("PASS-only summary (SN)", table([{"metric": r.get("key"), "value": r.get("value")}
                                                               for r in tools["bcftools_stats_pass"].get("SN", [])]))
        pv = os.path.join(o, "03_bcftools_stats", "plots")
        if os.path.isdir(pv):
            pngs = sorted(f for f in os.listdir(pv) if f.endswith(".png"))[:12]
            if pngs:
                b_html += details(f"plot-vcfstats figures ({len(pngs)})", '<div class="grid">' +
                                  "".join(img(os.path.join(pv, f), f) for f in pngs) + "</div>")
        sections.append(("bcftools", "bcftools stats", b_html))

    # ---------------------------------------------------------------- 11. vcftools / GATK / truth
    x_html = ""
    vt = tools.get("vcftools")
    if vt:
        x_html += "<h3>vcftools</h3>" + kv_table([("HWE sites tested", (vt.get("hwe") or {}).get("sites_tested")),
                                                  ("HWE p<1e-4", (vt.get("hwe") or {}).get("p_lt_1e-4")),
                                                  ("HWE p<1e-6", (vt.get("hwe") or {}).get("p_lt_1e-6")),
                                                  ("Median site mean depth", (vt.get("site_mean_depth") or {}).get("median"))])
        if vt.get("het"):
            x_html += "<h4>--het (inbreeding F)</h4>" + table(vt["het"], max_rows=300)
        if vt.get("imiss"):
            x_html += details("--missing-indv", table(vt["imiss"], max_rows=500))
        if vt.get("tstv_summary"):
            x_html += details("--TsTv-summary", table(vt["tstv_summary"]))
    for k, lab in (("gatk_cvcm_summary", "GATK CollectVariantCallingMetrics – summary"),
                   ("gatk_cvcm_detail", "GATK CollectVariantCallingMetrics – per sample")):
        if tools.get(k):
            rows = tools[k]
            x_html += f"<h3>{lab}</h3>" + (kv_table(list(rows[0].items())) if len(rows) == 1 else table(rows, max_rows=500))
    if tools.get("happy"):
        x_html += "<h3>Truth-set concordance (hap.py)</h3>" + table(tools["happy"])
    if tools.get("vcfeval"):
        x_html += "<h3>Truth-set concordance (RTG vcfeval)</h3>" + table(tools["vcfeval"])
    if x_html:
        sections.append(("external", "vcftools, GATK metrics & truth concordance", x_html))

    # ---------------------------------------------------------------- 12. methods
    st = read_tsv(os.path.join(o, "logs", "module_status.tsv"))
    m_html = "<h3>Module log</h3>" + (table(st) if st else '<p class="na">No module log.</p>')
    if REPORT_WARNINGS:
        m_html += "<h3>Report warnings</h3><pre>" + E("\n".join(REPORT_WARNINGS)) + "</pre>"
    tv = os.path.join(o, "logs", "tool_versions.txt")
    if os.path.exists(tv):
        m_html += "<h3>Tool versions</h3><pre>" + E(open(tv).read()) + "</pre>"
    pr = os.path.join(o, "logs", "params.txt")
    if os.path.exists(pr):
        m_html += "<h3>Parameters</h3><pre>" + E(open(pr).read()) + "</pre>"
    m_html += """<h3>Interpretation guide</h3><ul>
<li><b>Ts/Tv</b>: WGS germline ≈ 2.0–2.1, WES ≈ 2.8–3.3; values near 0.5 indicate random (artefactual) SNVs.</li>
<li><b>Het/HomAlt</b>: ≈1.5–2.0 for outbred WGS samples; high values suggest contamination, low values consanguinity or haploid regions.</li>
<li><b>Known rate</b>: &gt;99% of SNVs in a germline WGS call set are typically in dbSNP; lower known rates suggest artefacts or rare/somatic calls.</li>
<li><b>Indels</b>: in coding regions in-frame (multiple of 3) indels should dominate; insertion/deletion ratio ≈ 1 in modern callers.</li>
<li><b>Missingness</b>: per-sample &gt;5% or a long tail of high-missingness sites deserves inspection.</li></ul>"""
    sections.append(("methods", "Methods, tool versions & module log", m_html))

    # ---------------------------------------------------------------- write HTML
    nav = "".join(f'<a href="#{sid}">{i + 1}. {E(t)}</a>' for i, (sid, t, _) in enumerate(sections))
    body = "".join(f'<section id="{sid}"><h2>{i + 1}. {E(t)}</h2>{c}</section>' for i, (sid, t, c) in enumerate(sections))
    html_doc = TEMPLATE.format(title=E(a.title), vcf=E(os.path.basename(inp["path"])),
                               date=dt.datetime.now().strftime("%Y-%m-%d %H:%M"), kpis=kpi_html, nav=nav, body=body,
                               mode=E(prof["mode"]), callers=E(", ".join(prof["callers"][:4]) or "unknown caller"))
    with open(a.report, "w") as fh:
        fh.write(html_doc)
    print(f"[report] wrote {a.report} ({os.path.getsize(a.report) / 1e6:.1f} MB)", file=sys.stderr)


TEMPLATE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title>
<style>
:root{{--bg:#f5f5f3;--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#8a8984;--line:#e4e3df;--accent:#2a78d6;
--ok:#0b7a3e;--okbg:#e3f4ea;--warn:#8a5a00;--warnbg:#fdf1d6;--fail:#b3261e;--failbg:#fbe4e2;--na:#52514e;--nabg:#efeeea}}
@media (prefers-color-scheme:dark){{:root:not([data-theme="light"]){{--bg:#121211;--surface:#1a1a19;--ink:#fff;--ink2:#c3c2b7;
--muted:#8f8e86;--line:#33332f;--accent:#3987e5;--ok:#6fd49b;--okbg:#173524;--warn:#f2c46b;--warnbg:#3a2e12;--fail:#ff8a80;
--failbg:#3d1a17;--na:#c3c2b7;--nabg:#2a2a27}}}}
:root[data-theme="dark"]{{--bg:#121211;--surface:#1a1a19;--ink:#fff;--ink2:#c3c2b7;--muted:#8f8e86;--line:#33332f;--accent:#3987e5;
--ok:#6fd49b;--okbg:#173524;--warn:#f2c46b;--warnbg:#3a2e12;--fail:#ff8a80;--failbg:#3d1a17;--na:#c3c2b7;--nabg:#2a2a27}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}}
header{{padding:28px 32px 18px;border-bottom:1px solid var(--line);background:var(--surface)}}
header h1{{margin:0 0 4px;font-size:24px}}header .sub{{color:var(--ink2)}}
.layout{{display:flex;align-items:flex-start}}
nav{{position:sticky;top:0;width:250px;flex:none;max-height:100vh;overflow:auto;padding:18px 12px;border-right:1px solid var(--line)}}
nav a{{display:block;padding:5px 10px;color:var(--ink2);text-decoration:none;border-radius:6px;font-size:13px}}
nav a:hover{{background:var(--surface);color:var(--accent)}}
main{{flex:1;min-width:0;padding:8px 32px 60px;max-width:1300px}}
section{{background:var(--surface);border:1px solid var(--line);border-radius:12px;padding:20px 24px;margin:20px 0}}
h2{{margin:0 0 12px;font-size:19px}}h3{{margin:22px 0 8px;font-size:15px}}h4{{margin:16px 0 6px;font-size:13.5px;color:var(--ink2)}}
.kpis{{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:10px;margin-top:16px}}
.kpi{{background:var(--bg);border:1px solid var(--line);border-radius:10px;padding:10px 14px}}
.kpi .kv{{font-size:21px;font-weight:650;font-variant-numeric:tabular-nums}}.kpi .kl{{color:var(--ink2);font-size:12px}}
.tw{{overflow-x:auto;margin:8px 0}}table{{border-collapse:collapse;font-size:12.5px;font-variant-numeric:tabular-nums}}
th,td{{padding:5px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}
thead th{{position:sticky;top:0;background:var(--bg);cursor:pointer;white-space:nowrap;color:var(--ink2);font-weight:600}}
thead th:hover{{color:var(--accent)}}table.kv th{{color:var(--ink2);font-weight:500;width:300px}}
tbody tr:hover{{background:var(--bg)}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(420px,1fr));gap:14px}}
figure{{margin:10px 0;background:#fcfcfb;border:1px solid var(--line);border-radius:8px;padding:6px;max-width:760px}}
figure.wide{{max-width:1180px}}.grid figure{{max-width:none}}
figure img{{width:100%;display:block}}figcaption{{font-size:11.5px;color:#52514e;padding:4px 4px 0}}
.badge{{display:inline-block;padding:1px 8px;border-radius:999px;font-size:12px;white-space:nowrap}}
.badge.ok{{color:var(--ok);background:var(--okbg)}}.badge.warn{{color:var(--warn);background:var(--warnbg);white-space:normal}}
.badge.fail{{color:var(--fail);background:var(--failbg)}}.badge.na{{color:var(--na);background:var(--nabg)}}
.na{{color:var(--muted);font-style:italic}}.note{{color:var(--ink2);font-size:12.5px}}
pre{{background:var(--bg);padding:10px;border-radius:6px;overflow:auto;font-size:11.5px;max-height:420px}}
.mono{{font-family:ui-monospace,Menlo,monospace;font-size:12px;word-break:break-word}}
details{{margin:8px 0}}summary{{cursor:pointer;color:var(--accent)}}
code{{font-family:ui-monospace,Menlo,monospace;font-size:12px}}
@media (max-width:900px){{nav{{display:none}}main{{padding:8px 16px}}.grid{{grid-template-columns:1fr}}header{{padding:20px 16px}}
section{{padding:16px}}table.kv th{{width:auto}}}}
</style></head><body>
<header><h1>{title}</h1><div class="sub"><b>{vcf}</b> · {mode} · {callers} · generated {date}</div>{kpis}</header>
<div class="layout"><nav>{nav}</nav><main>{body}</main></div>
<script>
document.querySelectorAll("table.sortable").forEach(t=>{{t.querySelectorAll("thead th").forEach((th,i)=>{{let asc=true;
th.addEventListener("click",()=>{{const tb=t.tBodies[0];const rows=[...tb.rows];rows.sort((a,b)=>{{
const x=a.cells[i].dataset.v,y=b.cells[i].dataset.v;const nx=parseFloat(x),ny=parseFloat(y);
const r=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y);return asc?r:-r}});asc=!asc;rows.forEach(r=>tb.appendChild(r))}})}})}});
</script></body></html>"""

if __name__ == "__main__":
    main()
