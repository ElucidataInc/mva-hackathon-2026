#!/usr/bin/env python3
"""
Turn structural-variant calls into per-gene SV candidates for rank_variants.py.

Input (auto-detected, repeatable):
  * an AnnotSV TSV (e.g. merged Manta + cue calls annotated by AnnotSV), per-gene "split" lines:
      effect     Location/Location2, Overlapped_CDS_length/_percent, Frameshift (in-frame exon deletion
                 -> SV_LOF_POSSIBLE), Dist_nearest_SS/Nearest_SS_type (intronic SV near a splice site -> SV_SPLICE),
                 RE_gene (regulatory element of a gene -> SV_REGULATORY)
      frequency  B_*_AFmax for the SV type, gnomAD AF columns and gnomAD AF keys inside INFO
      genotype   sample GT, unless it is a placeholder (same GT-only value for every SV): then "not measured",
                 counted as one copy; exon depth (--exon-depth) can show a homozygous deletion
      confidence Gap/ENCODE blacklist/SegDup/Repeat at the breakpoints -> high/medium/low (low: score -1)
      evidence   ACMG_class + AnnotSV_ranking_score/criteria (whole SV: raise only coding effects),
                 P_*/po_P_* known pathogenic CNVs, po_B_*_allG benign overlap (score -1), Gene_count,
                 gene-level HI/TS/pLI/LOEUF/OMIM/GenCC inheritance (gene_info, gene_flag)
      transcript the MANE transcript line when MANE IDs are in INFO, else the most severe line
      duplicates records of one event (same gene, type, >=50% reciprocal overlap) are merged
  * a plain SV VCF: gene effects from --gtf (MANE) or the bundled MVA genes, frequency from --pop
Optional: --exon-depth (exon depth screen) adds depth support; flagged exons without an SV call
become depth-only CNV candidates.
Output: one row per SV per gene (sv_candidates.tsv).
"""
import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict

CLASS = {"SV_LOF": 6,            # deletion removing coding sequence (frameshift / large) or the whole gene
         "SV_LOF_POSSIBLE": 4,   # in-frame exon deletion, intragenic duplication, inversion/BND/insertion breakpoint in a gene
         "CNV_DEPTH": 3,         # exon depth suggests a CNV, no SV call
         "SV_SPLICE": 3,         # intronic SV next to a splice site (Dist_nearest_SS)
         "SV_REGULATORY": 2,     # upstream (promoter), UTR only, or a regulatory element of the gene (RE_gene)
         "SV_INTRONIC": 1,
         "SV_DOSAGE": 1}         # whole-gene duplication
csv.field_size_limit(sys.maxsize)


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def norm_chrom(c, style_chr):
    b = c[3:] if c.startswith("chr") else c
    b = "M" if b in ("MT", "M") else b
    return ("chr" + b) if style_chr else ("MT" if b == "M" else b)


def type_class(t, alt=""):
    t, a = (t or "").upper(), (alt or "").upper()
    if "BND" in t or "TRA" in t or "[" in a or "]" in a:
        return "BND"
    for k in ("DEL", "DUP", "INV", "INS"):          # INV-DEL / INVDEL -> DEL-like
        if k in t:
            return k
    for k in ("DEL", "DUP", "INV", "INS"):
        if f"<{k}" in a:
            return k
    return "CNV" if "CNV" in t else "OTHER"


def fnum(x):
    vals = []
    for v in re.split(r"[;,|&]", str(x or "")):
        try:
            vals.append(float(v))
        except ValueError:
            pass
    return max(vals) if vals else None


def gt_from(fmt, sample):
    if not fmt or not sample:
        return "./."
    d = dict(zip(fmt.split(":"), sample.split(":")))
    return d.get("GT", "./.")


def empty(x):
    return (x or "").strip() in ("", "NA", ".", "-")


def split_top(s, sep=";"):
    """split on sep outside parentheses: 'A (HI=3;morbid);B' -> ['A (HI=3;morbid)', 'B']"""
    out, depth, cur = [], 0, ""
    for ch in s or "":
        depth += (ch == "(") - (ch == ")")
        if ch == sep and depth == 0:
            out.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


def info_dict(s):
    return parse_info(s) if s and s not in (".", "NA") else {}


# AnnotSV columns that describe one gene / transcript: never copied from the full line to a split line
GENE_LEVEL = re.compile(r"^(Gene_name|Gene_count|Tx.*|Overlapped_.*|Frameshift|Exon_count|Location2?|Dist_nearest_SS|"
                        r"Nearest_SS_type|Intersect_.*|HI|TS|DDD_.*|ExAC_.*|GenCC_.*|NCBI_gene_ID|HGNC_gene_ID|"
                        r"OMIM_.*|LOEUF_bin|GnomAD_pLI|Closest_.*)$")
TRICKY_REPEAT = re.compile(r"simple|low_complex|satellite|sat|\(", re.I)


def breakpoint_confidence(r):
    """high / medium / low from AnnotSV breakpoint context (gaps, ENCODE blacklist, segdups, repeats)"""
    why = []
    for side in ("left", "right"):
        if not empty(r.get(f"Gap_{side}")):
            why.append(f"gap at {side} breakpoint")
        if not empty(r.get(f"ENCODE_blacklist_{side}")):
            why.append(f"ENCODE blacklist at {side} breakpoint" +
                       (f" ({r.get(f'ENCODE_blacklist_characteristics_{side}')})"
                        if not empty(r.get(f"ENCODE_blacklist_characteristics_{side}")) else ""))
    seg = [side for side in ("left", "right") if not empty(r.get(f"SegDup_{side}"))]
    rep = [f"{r.get(f'Repeat_type_{side}')} at {side}" for side in ("left", "right")
           if not empty(r.get(f"Repeat_type_{side}")) and TRICKY_REPEAT.search(r.get(f"Repeat_type_{side}"))]
    if why or len(seg) == 2:
        return "low", why + (["segmental duplication at both breakpoints"] if len(seg) == 2 else [])
    if seg or rep:
        return "medium", ([f"segmental duplication at {seg[0]} breakpoint"] if seg else []) + rep
    return "high", []


def gene_info(r):
    """compact gene-level evidence from AnnotSV split line + dosage flag (HI for losses, TS for gains)"""
    moi = " ".join(x for x in (r.get("OMIM_inheritance"), r.get("GenCC_moi")) if not empty(x))
    parts = []
    if moi:
        ms = sorted({m for m in ("AR", "AD", "XLR", "XLD", "XL") if re.search(rf"\b{m}\b", moi)} |
                    ({"AR"} if re.search("recessive", moi, re.I) else set()) |
                    ({"AD"} if re.search("dominant", moi, re.I) else set()))
        parts.append("MOI=" + ("/".join(ms) if ms else moi[:30]))
    for k, lab in (("HI", "HI"), ("TS", "TS"), ("GnomAD_pLI", "pLI"), ("LOEUF_bin", "LOEUF_bin")):
        if not empty(r.get(k)):
            parts.append(f"{lab}={r.get(k)}")
    dis = next((r.get(k) for k in ("OMIM_phenotype", "GenCC_disease") if not empty(r.get(k))), "")
    if dis:
        parts.append("disease=" + dis.replace("\t", " ")[:60])
    if not empty(r.get("GenCC_classification")):
        parts.append("GenCC=" + r.get("GenCC_classification")[:30])
    hi = fnum(r.get("HI"))
    ts = fnum(r.get("TS"))
    pli = fnum(r.get("GnomAD_pLI"))
    flag = []
    if hi == 3 or (pli is not None and pli >= 0.9) or "AD" in moi.replace("dominant", "AD"):
        flag.append("loss-sensitive")
    if ts == 3:
        flag.append("gain-sensitive")
    return "; ".join(parts) or "-", ",".join(flag)


# older / alternative symbols of the MVA genes that AnnotSV gene lists may use
MVA_ALIASES = {"CCDC84": "CENATAC", "FAM178A": "SLF2", "C10orf6": "SLF2", "BUBR1": "BUB1B", "BUB1L": "BUB1",
               "CMT2": "MAD2L1BP", "MUL": "TRIM37", "SMC5L1": "SMC5", "PCH2": "TRIP13", "TSP57": "CEP57", "MAD1": "MAD1L1"}


def mva_near(mva, chrom, start, end, flank=50_000):
    """MVA genes within flank bp of an SV: [(gene, relation)]"""
    out = []
    for g, d in (mva or {}).items():
        if d["chrom"] != chrom:
            continue
        if start <= d["end"] and end >= d["start"]:
            out.append((g, "overlaps gene"))
        else:
            dist = d["start"] - end if end < d["start"] else start - d["end"]
            if dist <= flank:
                out.append((g, f"{dist / 1000:.1f} kb from gene"))
    return out


def is_annotsv(path):
    with opener(path) as fh:
        first = fh.readline()
    return "AnnotSV_ID" in first and "Annotation_mode" in first


# ---------------------------------------------------------------- AnnotSV TSV
def classify_annotsv(r, t, ss_dist=50):
    loc = (r.get("Location") or "").strip()
    loc2 = (r.get("Location2") or "").strip()
    cds = fnum(r.get("Overlapped_CDS_length")) or 0
    cds_pct = fnum(r.get("Overlapped_CDS_percent"))
    fs = (r.get("Frameshift") or "").strip().lower()
    frameshift = fs == "yes"
    ss = fnum(r.get("Dist_nearest_SS"))
    ss_type = (r.get("Nearest_SS_type") or "").strip()
    near_ss = ss is not None and ss <= ss_dist
    whole = loc == "txStart-txEnd"
    m = re.fullmatch(r"intron(\d+)-intron(\d+)", loc)
    intronic = bool(m and m.group(1) == m.group(2))
    up = "upstream" in loc2.lower() or (not loc and "upstream" in (r.get("Location2") or "").lower())
    down = "downstream" in loc2.lower()
    exonic = cds > 0 or "exon" in loc or "txStart" in loc or "txEnd" in loc
    utr_only = exonic and cds == 0 and "UTR" in loc2.upper()
    if t in ("DEL", "CNV"):
        if whole:
            return "SV_LOF", "whole-gene deletion"
        if cds > 0:
            pct = f", {cds_pct:.0f}% of CDS" if cds_pct is not None else ""
            # in-frame loss of internal exon(s) of a small part of the protein: weaker than a frameshift (PVS1 strength reduced)
            if fs == "no" and "txStart" not in loc and "txEnd" not in loc and (cds_pct is None or cds_pct < 50):
                return "SV_LOF_POSSIBLE", f"in-frame deletion of coding sequence ({loc}; {int(cds)} bp CDS{pct})"
            return "SV_LOF", f"deletes coding sequence ({loc}; {int(cds)} bp CDS{pct}{', frameshift' if frameshift else ''})"
        if utr_only:
            return "SV_REGULATORY", f"deletes UTR only ({loc})"
        if intronic and near_ss:
            return "SV_SPLICE", f"intronic deletion {int(ss)} bp from a splice {ss_type or 'site'} ({loc})"
        if intronic:
            return "SV_INTRONIC", f"intronic deletion ({loc})"
        if up:
            return "SV_REGULATORY", "deletion upstream (promoter region)"
        if down:
            return "SV_INTRONIC", "deletion downstream of gene"
        return ("SV_LOF_POSSIBLE", f"deletion in gene ({loc})") if exonic else ("SV_INTRONIC", f"deletion ({loc or loc2})")
    if t == "DUP":
        if whole:
            return "SV_DOSAGE", "whole-gene duplication"
        if exonic and not up and not down:
            return "SV_LOF_POSSIBLE", f"partial/intragenic duplication ({loc}{', frameshift' if frameshift else ''})"
        if intronic and near_ss:
            return "SV_SPLICE", f"intronic duplication {int(ss)} bp from a splice {ss_type or 'site'} ({loc})"
        if intronic:
            return "SV_INTRONIC", f"intronic duplication ({loc})"
        return "SV_REGULATORY", f"duplication {('upstream' if up else 'near gene')} ({loc or loc2})"
    # INV, BND, INS
    if whole and t == "INV":
        return "SV_REGULATORY", "inversion spanning whole gene"
    if exonic and not up and not down:
        return "SV_LOF_POSSIBLE", f"{t.lower()} breakpoint/insertion in coding gene region ({loc})"
    if intronic and near_ss:
        return "SV_SPLICE", f"intronic {t.lower()} {int(ss)} bp from a splice {ss_type or 'site'} ({loc})"
    if intronic:
        return "SV_INTRONIC", f"intronic {t.lower()} ({loc})"
    return "SV_REGULATORY", f"{t.lower()} {('upstream' if up else 'near gene')} ({loc or loc2})"


def read_annotsv(path, style_chr, max_af, sample, gt_mode="auto", ss_dist=50, models=None, models_name="",
                 re_max_len=1_000_000, mva=None, fate=None):
    rows = []
    dropped = defaultdict(int)
    with opener(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        cols = rdr.fieldnames or []
        fmt_i = cols.index("FORMAT") if "FORMAT" in cols else None
        samp_col = sample if sample in cols else (cols[fmt_i + 1] if fmt_i is not None and fmt_i + 1 < len(cols)
                                                   and cols[fmt_i + 1] != "Annotation_mode" else None)
        gnomad_cols = [c for c in cols if re.search("gnomad", c, re.I) and re.search("AF", c)]
        supp_cols = [c for c in cols if re.fullmatch(r"(SUPP|SUPP_VEC|callers?|CALLERS?|SOURCE|Caller)", c)]
        all_rows = list(rdr)
    split = [r for r in all_rows if r.get("Annotation_mode") == "split"]
    # split lines, plus the full line of any SV without split lines (intergenic SVs: RE_gene can still name a gene)
    with_split = {r.get("AnnotSV_ID") for r in split}
    use = split + [r for r in all_rows if r.get("Annotation_mode") != "split" and r.get("AnnotSV_ID") not in with_split]
    # SV-wide values (ACMG class, ranking score, INFO, breakpoint context, benign/pathogenic overlaps) are on the
    # "full" line of each SV: copy them to its split lines, but never gene/transcript-level fields
    full = {r.get("AnnotSV_ID"): r for r in all_rows if r.get("Annotation_mode") == "full"}
    n_acmg_from_full = 0
    for r in split:
        f = full.get(r.get("AnnotSV_ID"))
        if not f:
            continue
        for c, v in f.items():
            if c and not GENE_LEVEL.match(c) and empty(r.get(c)) and not empty(v):
                r[c] = v
                n_acmg_from_full += c == "ACMG_class"

    # genotype: AnnotSV needs a sample column, so pipelines often insert GT=0/1 for every SV. That is a
    # placeholder, not a measurement (a homozygous deletion would look heterozygous).
    gts = {gt_from(r.get("FORMAT"), r.get(samp_col)) for r in all_rows} if samp_col else set()
    fmts = {(r.get("FORMAT") or "") for r in all_rows}
    placeholder = gt_mode == "placeholder" or (gt_mode == "auto" and (not samp_col or (len(gts) <= 1 and fmts <= {"GT", ""})))
    print(f"[sv] AnnotSV: {len(all_rows):,} lines ({len(split):,} split, {len(full):,} full); sample column: {samp_col or 'none'}; "
          f"genotype: {'PLACEHOLDER (not measured; ' + (','.join(sorted(gts)) or 'none') + ') - SVs counted as one copy unless depth shows two' if placeholder else 'measured'}; "
          f"frequency: B_*_AFmax{(' + ' + ', '.join(gnomad_cols)) if gnomad_cols else ''} + gnomAD AF keys in INFO; "
          f"support column: {', '.join(supp_cols) or 'none (SUPP/SUPP_VEC/IDLIST in INFO used if present)'}", file=sys.stderr)
    if split:
        print(f"[sv] SV-wide fields copied from full lines (ACMG class on {n_acmg_from_full:,} split line(s))", file=sys.stderr)
    n_fullonly = sum(1 for r in use if r.get("Annotation_mode") != "split")
    if n_fullonly:
        # full lines carry no per-gene Location / CDS overlap: compute the effect on each gene from gene models
        print(f"[sv] {n_fullonly:,} SV(s) have only a full line (no per-gene effect columns): effect computed from "
              f"{models_name}; genes without a model are skipped. Best fix: run AnnotSV with split lines "
              f"(-annotationMode both, the default){'' if 'GTF' in models_name else ' or give --gene-gtf (MANE GTF)'}",
              file=sys.stderr)

    n_conf = defaultdict(int)
    fate = {} if fate is None else fate        # AnnotSV_ID -> record of SVs near MVA genes (what happened to them)
    for r in use:
        sid = r.get("AnnotSV_ID", "")
        try:
            c0 = norm_chrom(r.get("SV_chrom", ""), style_chr)
            s0 = int(float(r.get("SV_start") or 0))
            e0 = int(float(r.get("SV_end") or s0 + 1))
        except ValueError:
            c0, s0, e0 = "", 0, 0
        near = mva_near(mva, c0, s0, e0)
        if near and sid not in fate:
            fate[sid] = {"chrom": c0, "start": s0, "end": e0, "svtype": r.get("SV_type", ""), "near": near,
                         "reason": "", "gene_name": r.get("Gene_name", ""), "id": r.get("ID", "")}

        def drop(reason):
            dropped[reason] += 1
            if sid in fate and not fate[sid]["reason"]:
                fate[sid]["reason"] = reason
        t = type_class(r.get("SV_type"), r.get("ALT"))
        if t == "OTHER":
            drop("unknown SV type")
            continue
        info = info_dict(r.get("INFO"))
        genes = [MVA_ALIASES.get(g, g) for g in re.split(r"[;,]", r.get("Gene_name") or "") if g]
        # MVA genes hit by coordinates (gene body or 10 kb upstream) even if AnnotSV does not name them
        coord_mva = [g for g, rel in near if g not in genes and annotate_gene({"start": s0, "end": e0, "type": t}, mva[g])]
        # regulatory elements (RE_gene, full line): genes whose enhancer/promoter is hit but not their body
        f = full.get(r.get("AnnotSV_ID")) or r
        re_genes = []
        for item in split_top(f.get("RE_gene") or ""):
            m = re.match(r"([A-Za-z0-9][\w.\-]*)\s*(\(([^)]*)\))?", item)
            if m and not empty(m.group(1)):
                re_genes.append((m.group(1), m.group(3) or ""))
        if not genes and not re_genes and not coord_mva:
            drop("no gene")
            continue

        # population frequency: AnnotSV benign-SV AFmax for the SV type + gnomAD AF columns + gnomAD keys in INFO
        key = {"DEL": "B_loss_AFmax", "CNV": "B_loss_AFmax", "DUP": "B_gain_AFmax",
               "INS": "B_ins_AFmax", "INV": "B_inv_AFmax"}.get(t)
        afs = [fnum(r.get(key))] if key else []
        afs += [fnum(r.get(c)) for c in gnomad_cols]
        afs += [fnum(v) for k, v in info.items() if v is not True and re.search("gnomad", k, re.I) and re.search("AF", k)]
        afs = [x for x in afs if x is not None]
        af = max(afs) if afs else None
        if af is not None and af >= max_af:
            drop(f"common SV (population AF {af:.3g} >= {max_af})" if sid in fate else f"common SV (population AF >= {max_af})")
            continue

        start = int(float(r.get("SV_start") or 0))
        end = int(float(r.get("SV_end") or start + 1))
        gt = "./." if placeholder else (gt_from(r.get("FORMAT"), r.get(samp_col)) if samp_col else "./.")
        acmg = (r.get("ACMG_class") or "").strip()
        acmg = "" if empty(acmg) else acmg.replace("full=", "")
        conf, conf_why = breakpoint_confidence(r)
        n_conf[conf] += 1

        # SV-wide notes (same for every gene of this SV)
        sv_notes, adj = [], 0
        if acmg:
            crit = r.get("AnnotSV_ranking_criteria") or r.get("ACMG") or ""
            sv_notes.append(f"AnnotSV ACMG class {acmg} for the whole SV" +
                            (f" (score {r.get('AnnotSV_ranking_score')}" + (f"; {crit[:80]}" if not empty(crit) else "") + ")"
                             if not empty(r.get("AnnotSV_ranking_score")) else ""))
        if conf != "high":
            sv_notes.append(f"breakpoint confidence {conf}: " + ", ".join(conf_why))
            adj -= conf == "low"
        if "IMPRECISE" in info:
            sv_notes.append("imprecise breakpoints")
        p_key = {"DEL": "loss", "CNV": "loss", "DUP": "gain", "INS": "ins"}.get(t)
        if p_key:
            for pre, what in (("P", "overlaps a known pathogenic"), ("po_P", "partly overlaps a known pathogenic")):
                phen = r.get(f"{pre}_{p_key}_phen")
                if not empty(phen):
                    pct = r.get(f"{pre}_{p_key}_percent")
                    sv_notes.append(f"{what} {p_key}: {phen[:60]} [{r.get(f'{pre}_{p_key}_source') or '?'}"
                                    + (f"; {pct}%" if not empty(pct) else "") + "]")
                    break
            allg = r.get(f"po_B_{p_key}_allG_source") if p_key in ("loss", "gain") else ""
            if not empty(allg):
                sv_notes.append(f"benign {p_key}es ({allg}) cover all genes of this SV")
                adj -= 1
        n_genes = fnum(r.get("Gene_count"))
        if n_genes and n_genes > 5:
            sv_notes.append(f"multi-gene SV ({int(n_genes)} genes, {r.get('CytoBand') or ''}): also apply ACMG CNV scoring")
        if placeholder:
            sv_notes.append("genotype not measured (AnnotSV placeholder)")

        # caller support: columns, else SUPP / SUPP_VEC / IDLIST in INFO (SURVIVOR / Jasmine merges)
        supp = ";".join(f"{c}={r.get(c)}" for c in supp_cols if r.get(c))
        if not supp:
            supp = ";".join(f"{k}={info[k]}" for k in ("SUPP", "SUPP_VEC", "IDLIST") if k in info and info[k] is not True)
        supp_raw = r.get("SUPP") or (info.get("SUPP") if info.get("SUPP") is not True else "")
        n_supp = fnum(supp_raw) if not empty(supp_raw) else None
        if n_supp is None:
            vec = (r.get("SUPP_VEC") or (info.get("SUPP_VEC") if info.get("SUPP_VEC") is not True else "") or "").strip()
            n_supp = vec.count("1") if re.fullmatch(r"[01]+", vec) else None
        if n_supp and n_supp >= 2:
            sv_notes.append(f"called by {int(n_supp)} callers")
        elif n_supp == 1:
            sv_notes.append("called by one caller only")
        # MANE transcripts named in INFO (prefer those annotation lines)
        mane = {re.sub(r"\.\d+$", "", x) for k, v in info.items() if v is not True and re.search("mane", k, re.I)
                for x in re.findall(r"(?:NM_|ENST)\d+(?:\.\d+)?", v)}
        flt = r.get("FILTER", "PASS") or "PASS"
        base = {"chrom": norm_chrom(r.get("SV_chrom", ""), style_chr), "start": start, "end": end, "svtype": t,
                "svlen": abs(int(float(r.get("SV_length") or (end - start)))), "gt": gt,
                "gt_source": "placeholder" if placeholder else "measured", "callers": supp or "AnnotSV input",
                "filter": flt, "qual": r.get("QUAL", "."), "pr": "", "sr": "", "pop_af": "" if af is None else f"{af:.3g}",
                "confidence": conf, "acmg": acmg or "-", "cytoband": r.get("CytoBand") or "-",
                "ids": r.get("AnnotSV_ID", "") + (f"|{r.get('ID')}" if not empty(r.get("ID")) else "")}

        fullonly = r.get("Annotation_mode") != "split"
        gene_effects = {}
        if genes and fullonly:
            for g in genes:
                d = (models or {}).get(g)
                if not d or norm_chrom(d["chrom"], style_chr) != base["chrom"]:
                    dropped["gene without a model (full-line-only SV)"] += 1
                    continue
                res = annotate_gene({"start": start, "end": end, "type": t}, d)
                if res:
                    gene_effects[g] = (res[0], res[1] + " [from gene model]")
        elif genes:
            cls, what = classify_annotsv(r, t, ss_dist)
            gene_effects = {g: (cls, what) for g in genes}
        for g in coord_mva:                            # MVA gene not named by AnnotSV: effect from the bundled model
            if g not in gene_effects:
                res = annotate_gene({"start": start, "end": end, "type": t}, mva[g])
                if res:
                    gene_effects[g] = (res[0], res[1] + " [MVA gene by coordinates; not named by AnnotSV]")
        for g, (cls, what) in gene_effects.items():
            gnotes = list(sv_notes)
            if acmg:
                # the ACMG class describes the whole SV (it may be driven by another gene), so it only
                # strengthens the call where THIS gene's coding sequence is affected
                if acmg.startswith("5") and t in ("DEL", "CNV") and cls == "SV_LOF_POSSIBLE":
                    cls = "SV_LOF"
                    gnotes.append("raised to SV_LOF: ACMG class 5 deletion affecting this gene")
                elif acmg.startswith(("4", "5")) and cls in ("SV_REGULATORY", "SV_INTRONIC", "SV_DOSAGE", "SV_SPLICE"):
                    gnotes.append(f"class not raised: effect on this gene is {what}")
            # gene-level columns on a full line are for all genes together: only use them for a single-gene SV
            ginfo, gflag = gene_info(r) if (not fullonly or len(genes) == 1) else ("-", "")
            tx = re.sub(r"\.\d+$", "", r.get("Tx") or "")
            rows.append(dict(base, gene=g, **{"class": cls}, score=max(1, CLASS[cls] + adj), effect=what,
                             exons=r.get("Location", "") if not fullonly else "", notes=gnotes, gene_info=ginfo,
                             gene_flag=gflag or "-", _mane=bool(tx and tx in mane)))
        for g, det in re_genes:
            # regulatory hits of a large SV add little next to its gene content and flood the gene list
            if g in genes or t not in ("DEL", "DUP", "CNV", "INV") or base["svlen"] > re_max_len:
                if base["svlen"] > re_max_len and g not in genes:
                    dropped[f"RE_gene of SV > {re_max_len // 1000} kb (not used)"] += 1
                continue
            rows.append(dict(base, gene=g, **{"class": "SV_REGULATORY"}, score=max(1, CLASS["SV_REGULATORY"] + adj),
                             effect=f"{t.lower()} hits a regulatory element of the gene" + (f" ({det[:60]})" if det else ""),
                             exons="", notes=list(sv_notes), gene_info="-", gene_flag="-", _mane=False))

    # one row per SV and gene: the MANE transcript line if named in INFO, else the most severe line
    best = {}
    for x in rows:
        k = (x["gene"], x["chrom"], x["start"], x["end"], x["svtype"])
        if k not in best or (x["_mane"], x["score"]) > (best[k]["_mane"], best[k]["score"]):
            best[k] = x
    rows = sorted(best.values(), key=lambda x: (x["gene"], x["chrom"], x["start"]))
    # the same event reported twice (e.g. Manta and cue records not merged): same gene, compatible type,
    # >=50% reciprocal overlap -> one row, otherwise it would look like two alleles (false compound het)
    out = []
    for x in rows:
        hit = next((y for y in out[-20:] if y["gene"] == x["gene"] and y["chrom"] == x["chrom"] and
                    compatible({"type": x["svtype"]}, {"type": y["svtype"]}) and
                    (recip(x, y) >= 0.5 if x["svtype"] not in ("INS", "BND") else abs(x["start"] - y["start"]) <= 500)), None)
        if hit:
            keep, other = (hit, x) if hit["score"] >= x["score"] else (x, hit)
            keep["ids"] = hit["ids"] + "," + x["ids"]
            keep["notes"].append(f"same event as {other['svtype']} {other['start']}-{other['end']} (overlapping records merged)")
            if other["filter"] == "PASS":
                keep["filter"] = "PASS"
            if keep is x:
                out[out.index(hit)] = x
            dropped["duplicate SV records merged"] += 1
        else:
            out.append(x)
    for x in out:
        x.pop("_mane", None)
    print(f"[sv] breakpoint confidence: " + ", ".join(f"{k} {v:,}" for k, v in sorted(n_conf.items())), file=sys.stderr)
    for sid, fx in fate.items():
        if not fx["reason"]:
            fx["reason"] = "no effect on the MVA gene (outside gene body and 10 kb upstream)"
    return out, dropped, placeholder


# ---------------------------------------------------------------- plain SV VCF
def parse_info(s):
    d = {}
    for kv in s.split(";"):
        k, _, v = kv.partition("=")
        d[k] = v if _ else True
    return d


def read_sv_vcf(label, path, style_chr):
    out = []
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = parse_info(f[7])
            t = type_class(info.get("SVTYPE"), f[4])
            if t == "OTHER":
                continue
            if t == "BND" and ("INV3" in info or "INV5" in info):
                t = "INV"
            if t == "BND" and "MATEID" in info and f[2] > str(info["MATEID"]):
                continue
            pos = int(f[1])
            try:
                end = int(info.get("END", pos))
            except ValueError:
                end = pos
            try:
                svlen = abs(int(str(info.get("SVLEN")).split(",")[0]))
            except (TypeError, ValueError):
                svlen = abs(end - pos)
            if t in ("DEL", "DUP", "INV", "CNV") and end <= pos and svlen:
                end = pos + svlen
            if t in ("INS", "BND"):
                end = pos + 1
            gt, pr, sr = "./.", "", ""
            if len(f) > 9:
                d = dict(zip(f[8].split(":"), f[9].split(":")))
                gt, pr, sr = d.get("GT", "./."), d.get("PR", ""), d.get("SR", "")
            out.append({"caller": label, "chrom": norm_chrom(f[0], style_chr), "start": pos, "end": max(end, pos + 1),
                        "type": t, "svlen": svlen or (end - pos), "filter": f[6], "qual": f[5], "gt": gt,
                        "pr": pr, "sr": sr, "id": f[2], "imprecise": "IMPRECISE" in info})
    return out


def recip(a, b):
    ov = min(a["end"], b["end"]) - max(a["start"], b["start"])
    return 0.0 if ov <= 0 else min(ov / max(1, a["end"] - a["start"]), ov / max(1, b["end"] - b["start"]))


def compatible(a, b):
    return a["type"] == b["type"] or {a["type"], b["type"]} <= {"DEL", "CNV"} or {a["type"], b["type"]} <= {"DUP", "CNV"}


def merge(calls):
    calls.sort(key=lambda c: (c["chrom"], c["start"]))
    merged = []
    for c in calls:
        hit = None
        for m in merged[-50:]:
            if m["chrom"] == c["chrom"] and compatible(m, c) and \
               ((c["type"] in ("INS", "BND") and abs(m["start"] - c["start"]) <= 500) or recip(m, c) >= 0.5):
                hit = m
                break
        if hit:
            hit["callers"].add(c["caller"])
            hit["members"].append(c)
            if c["filter"] == "PASS":
                hit["filter"] = "PASS"
            hit["imprecise"] = hit["imprecise"] and c["imprecise"]
            if hit["gt"] in ("./.", ".") and c["gt"] not in ("./.", "."):
                hit["gt"] = c["gt"]
        else:
            m = dict(c)
            m["callers"], m["members"] = {c["caller"]}, [c]
            merged.append(m)
    return merged


def load_genes(gtf, style_chr):
    genes, rows, mane = {}, [], False
    with opener(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            g = re.search(r'gene_name "([^"]+)"', f[8])
            if not g:
                continue
            is_mane = 'tag "MANE_Select"' in f[8]
            mane = mane or is_mane
            rows.append((g.group(1), f[2], norm_chrom(f[0], style_chr), int(f[3]), int(f[4]), f[6], is_mane))
    for g, kind, c, s, e, strand, is_mane in rows:
        if mane and not is_mane:
            continue
        d = genes.setdefault(g, {"chrom": c, "start": s, "end": e, "strand": strand, "exons": [], "cds": []})
        if d["chrom"] != c:
            continue
        d["start"], d["end"] = min(d["start"], s), max(d["end"], e)
        (d["exons"] if kind == "exon" else d["cds"]).append((s, e))
    for d in genes.values():
        d["exons"], d["cds"] = sorted(set(d["exons"])), sorted(set(d["cds"]))
    return genes


def load_bundled(here, style_chr):
    genes = {}
    for line in open(os.path.join(here, "mva_genes_GRCh38.bed")):
        if not line.startswith("#"):
            c, s, e, g, _, strand = line.rstrip("\n").split("\t")[:6]
            genes[g] = {"chrom": norm_chrom(c, style_chr), "start": int(s) + 1, "end": int(e), "strand": strand,
                        "exons": [], "cds": []}
    for line in open(os.path.join(here, "mva_exons_GRCh38.bed")):
        if not line.startswith("#"):
            c, s, e, name = line.rstrip("\n").split("\t")[:4]
            g = name.split("_exon")[0]
            if g in genes:
                genes[g]["exons"].append((int(s) + 1, int(e)))
    for d in genes.values():
        d["exons"].sort()
        d["cds"] = list(d["exons"])
    return genes


def annotate_gene(sv, d, upstream=10000):
    s, e, t = sv["start"], sv["end"], sv["type"]
    if t in ("INS", "BND"):
        e = s
    in_gene = s <= d["end"] and e >= d["start"]
    up_s, up_e = (d["start"] - upstream, d["start"] - 1) if d["strand"] == "+" else (d["end"] + 1, d["end"] + upstream)
    in_up = s <= up_e and e >= up_s
    whole = s <= d["start"] and e >= d["end"]
    if not in_gene and not in_up and not whole:
        return None
    ex = d["exons"] if d["strand"] == "+" else list(reversed(d["exons"]))
    exons = [i + 1 for i, (a, b) in enumerate(ex) if a <= e and b >= s]
    cds_hit = any(a <= e and b >= s for a, b in d["cds"])
    span = (f"exon {exons[0]}" if len(exons) == 1 else f"exons {exons[0]}-{exons[-1]}") + f" of {len(ex)}" if exons else ""
    if t in ("DEL", "CNV"):
        if whole:
            return "SV_LOF", "whole-gene deletion", exons
        if cds_hit:
            return "SV_LOF", f"deletes {span}", exons
        if exons:
            return "SV_REGULATORY", "deletes UTR exon(s) only", exons
        return ("SV_INTRONIC", "intronic deletion", exons) if in_gene else ("SV_REGULATORY", "deletion upstream (promoter)", exons)
    if t == "DUP":
        if whole:
            return "SV_DOSAGE", "whole-gene duplication", exons
        if exons:
            return "SV_LOF_POSSIBLE", f"intragenic/partial duplication of {span}", exons
        return ("SV_INTRONIC", "intronic duplication", exons) if in_gene else ("SV_REGULATORY", "duplication upstream", exons)
    if t == "INV" and whole:
        return "SV_REGULATORY", "inversion spanning whole gene", exons
    if in_gene and (exons or t in ("INV", "BND")):
        return "SV_LOF_POSSIBLE", f"{t.lower()} breakpoint/insertion in gene" + (f" ({span})" if span else ""), exons
    return ("SV_INTRONIC", f"intronic {t.lower()}", exons) if in_gene else ("SV_REGULATORY", f"{t.lower()} upstream", exons)


def pop_af(sv, pop, max_dist=500):
    if not pop:
        return None
    c, s, e, best = sv["chrom"], sv["start"], sv["end"], None
    for chrom in (c, c[3:] if c.startswith("chr") else "chr" + c):
        try:
            r = subprocess.run(["bcftools", "query", "-r", f"{chrom}:{max(1, s - max_dist)}-{e + max_dist}", "-f",
                                "%POS\t%INFO/END\t%INFO/SVTYPE\t%ALT\t%INFO/AF\t%INFO/GRPMAX_AF\n", pop],
                               capture_output=True, text=True, timeout=120)
        except Exception:
            continue
        if r.returncode != 0:
            continue
        for line in r.stdout.splitlines():
            p, en, ty, alt, af, gaf = (line.split("\t") + [""] * 6)[:6]
            try:
                rec = {"start": int(p), "end": int(en) if en not in (".", "") else int(p) + 1, "type": type_class(ty, alt)}
            except ValueError:
                continue
            if not (compatible(rec, sv) or (sv["type"] in ("DEL", "DUP") and rec["type"] == "CNV")):
                continue
            if (abs(rec["start"] - s) <= max_dist) if sv["type"] in ("INS", "BND") else recip(rec, sv) >= 0.5:
                v = [x for x in (fnum(gaf), fnum(af)) if x is not None]
                if v:
                    best = max(best or 0.0, max(v))
        if r.stdout:
            break
    return best


def read_vcf_svs(specs, gtf, pop, style_chr, max_af, here):
    calls = []
    for label, path in specs:
        calls += read_sv_vcf(label, path, style_chr)
    merged = merge(calls)
    genes = load_genes(gtf, style_chr) if gtf else load_bundled(here, style_chr)
    print(f"[sv] VCF input: {len(merged):,} merged SVs; gene models: {len(genes):,} ({'GTF' if gtf else 'bundled MVA genes'})",
          file=sys.stderr)
    by_chrom = defaultdict(list)
    for g, d in genes.items():
        by_chrom[d["chrom"]].append((d["start"] - 10000, d["end"] + 10000, g))
    rows, dropped = [], defaultdict(int)
    for sv in merged:
        cand = [g for (s, e, g) in by_chrom.get(sv["chrom"], []) if s <= sv["end"] and e >= sv["start"]]
        if not cand:
            dropped["no gene within 10 kb"] += 1
            continue
        af = pop_af(sv, pop)
        if af is not None and af >= max_af:
            dropped[f"common in gnomAD-SV (AF >= {max_af})"] += 1
            continue
        for g in cand:
            res = annotate_gene(sv, genes[g])
            if not res:
                continue
            cls, what, exons = res
            notes = []
            if len(sv["callers"]) > 1:
                notes.append("called by " + "+".join(sorted(sv["callers"])))
            if sv["filter"] != "PASS":
                notes.append(f"FILTER={sv['filter']}")
            if sv["imprecise"]:
                notes.append("imprecise breakpoints")
            m0 = sv["members"][0]
            rows.append({"gene": g, "chrom": sv["chrom"], "start": sv["start"], "end": sv["end"], "svtype": sv["type"],
                         "svlen": sv["svlen"], "gt": sv["gt"], "callers": ",".join(sorted(sv["callers"])),
                         "filter": sv["filter"], "qual": m0["qual"], "pr": m0["pr"], "sr": m0["sr"],
                         "pop_af": "" if af is None else f"{af:.3g}", "class": cls, "score": CLASS[cls], "effect": what,
                         "exons": ",".join(map(str, exons)), "notes": notes,
                         "ids": ",".join(f"{x['caller']}:{x['id']}" for x in sv["members"])})
    return rows, dropped


# ---------------------------------------------------------------- read depth inside SVs (BAM)
def sv_depth(rows, bam, ref=None, min_len=300, threads=1):
    """Depth inside each DEL/DUP relative to its own flanks (samtools bedcov, MAPQ>=1) -> genotype.
    Local flanks cancel most GC / mappability bias and need no genome-wide mean."""
    def bam_style():
        r = subprocess.run(["samtools", "view", "-H", bam] + (["--reference", ref] if ref else []),
                           capture_output=True, text=True)
        names = {f[3:] for line in r.stdout.splitlines() if line.startswith("@SQ")
                 for f in line.split("\t") if f.startswith("SN:")}
        return names
    names = bam_style()
    if not names:
        print("[sv] depth: could not read the BAM header; SV depth genotypes skipped", file=sys.stderr)
        return 0
    def to_bam(c):
        for cand in (c, c[3:] if c.startswith("chr") else "chr" + c, "MT" if c in ("chrM", "M") else "chrM" if c == "MT" else c):
            if cand in names:
                return cand
        return None
    svs = {}
    for x in rows:
        if x["svtype"] in ("DEL", "DUP", "CNV") and int(x["svlen"] or 0) >= min_len and x.get("callers") != "depth":
            svs.setdefault((x["chrom"], int(x["start"]), int(x["end"]), x["svtype"]), []).append(x)
    if not svs:
        return 0
    bed = []
    for i, (c, s_, e, t) in enumerate(svs):
        bc = to_bam(c)
        if not bc:
            continue
        ln = e - s_
        inner = [(s_, e)] if ln <= 100_000 else [(s_ + k * ln // 10, s_ + k * ln // 10 + 10_000) for k in range(10)]
        f, gap = min(max(ln, 1000), 10_000), 200
        parts = [("in", a, b) for a, b in inner] + [("L", max(0, s_ - gap - f), max(1, s_ - gap)), ("R", e + gap, e + gap + f)]
        bed += [f"{bc}\t{a}\t{b}\t{i}|{p}" for p, a, b in parts if b > a]
    cmd = ["samtools", "bedcov", "-Q", "1"] + (["--reference", ref] if ref else []) + ["/dev/stdin", bam]
    r = subprocess.run(cmd, input="\n".join(bed) + "\n", capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[sv] depth: samtools bedcov failed ({r.stderr.strip()[:200]}); SV depth genotypes skipped", file=sys.stderr)
        return 0
    acc = defaultdict(lambda: [0.0, 0])
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 5:
            acc[f[3]][0] += float(f[4])
            acc[f[3]][1] += int(f[2]) - int(f[1])
    mean = {k: (v[0] / v[1] if v[1] else None) for k, v in acc.items()}
    n = 0
    for i, (key, xs) in enumerate(svs.items()):
        inside, left, right = mean.get(f"{i}|in"), mean.get(f"{i}|L"), mean.get(f"{i}|R")
        fl = [v for v in (left, right) if v is not None and v >= 5]          # a flank in a gap/telomere is ignored
        if inside is None or not fl:
            for x in xs:
                x["notes"].append("depth: flanks not covered, no depth genotype")
            continue
        ratio = inside / (sum(fl) / len(fl))
        t = key[3]
        gt, call, adj = None, "", 0
        if t in ("DEL", "CNV") and (t == "DEL" or ratio < 1):
            if ratio < 0.15:
                gt, call = "1/1", "HOMOZYGOUS deletion"
            elif ratio < 0.70:
                gt, call = "0/1", "heterozygous deletion"
            elif ratio < 0.85:
                call = "partial depth loss (mosaic, or a false/smaller call)"
            else:
                call, adj = "deletion NOT supported by depth", -1 if key[2] - key[1] >= 1000 else 0
        else:
            if ratio >= 1.8:
                gt, call = "1/1", "duplication, >=2 extra copies (homozygous or amplified)"
            elif ratio >= 1.3:
                gt, call = "0/1", "heterozygous duplication"
            elif ratio >= 1.15:
                call = "partial depth gain (mosaic, or a false/smaller call)"
            else:
                call, adj = "duplication NOT supported by depth", -1 if key[2] - key[1] >= 1000 else 0
        for x in xs:
            x["depth_ratio"] = f"{ratio:.2f}"
            x["notes"].append(f"depth {ratio:.2f}x of flanks: {call}")
            if gt and x.get("gt_source") != "measured":
                x["gt"], x["gt_source"] = gt, "depth"
                x["notes"] = [nt for nt in x["notes"] if not nt.startswith("genotype not measured")]
            elif gt and x.get("gt_source") == "measured" and gt != x["gt"].replace("|", "/"):
                x["notes"].append(f"caller GT {x['gt']} differs from depth ({gt})")
            if adj:
                x["score"] = max(1, int(x["score"]) + adj)
        n += 1
    return n


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sv", action="append", default=[], help="[label=]path: AnnotSV TSV or SV VCF (repeatable)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--gtf", help="gene models for plain SV VCFs (MANE GTF)")
    ap.add_argument("--pop", help="gnomAD-SV sites VCF for plain SV VCFs")
    ap.add_argument("--exon-depth", help="exon depth table (02_coverage/exon_depth.tsv)")
    ap.add_argument("--sample", help="sample column name in the AnnotSV TSV [first sample column]")
    ap.add_argument("--bam", help="sample BAM/CRAM: depth inside each DEL/DUP vs its flanks -> genotype (samtools)")
    ap.add_argument("--bam-ref", help="FASTA for a CRAM")
    ap.add_argument("--depth-min-len", type=int, default=300, help="shortest DEL/DUP to genotype by depth [300]")
    ap.add_argument("--re-max-len", type=int, default=1_000_000, help="use RE_gene (regulatory element) hits only for SVs up to this length [1000000]")
    ap.add_argument("--sv-gt", choices=["auto", "placeholder", "measured"], default="auto",
                    help="AnnotSV sample genotype: auto = placeholder if every SV has the same GT-only genotype [auto]")
    ap.add_argument("--ss-dist", type=int, default=50, help="intronic SV this close (bp) to a splice site -> SV_SPLICE [50]")
    ap.add_argument("--max-af", type=float, default=0.01)
    ap.add_argument("--style-chr", type=int, default=1)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))

    rows, dropped = [], defaultdict(int)
    vcf_specs = []
    models, models_name = None, ""
    mva_models = load_bundled(here, a.style_chr)
    mva_fate = {}
    depth_dropped = {}
    for spec in a.sv:
        label, _, path = spec.partition("=")
        if not path:
            label, path = os.path.basename(spec).split(".")[0], spec
        if not os.path.exists(path):
            sys.exit(f"[sv] ERROR: SV file not found: {path}")
        if is_annotsv(path):
            if models is None:
                models = load_genes(a.gtf, a.style_chr) if a.gtf else load_bundled(here, a.style_chr)
                models_name = f"GTF {os.path.basename(a.gtf)} ({len(models):,} genes)" if a.gtf else \
                    f"bundled MVA gene models only ({len(models)} genes)"
            r, d, _ = read_annotsv(path, a.style_chr, a.max_af, a.sample, a.sv_gt, a.ss_dist, models, models_name,
                                   a.re_max_len, mva_models, mva_fate)
            rows += r
            for k, v in d.items():
                dropped[k] += v
        else:
            vcf_specs.append((label, path))
    if vcf_specs:
        r, d = read_vcf_svs(vcf_specs, a.gtf, a.pop, a.style_chr, a.max_af, here)
        rows += r
        for k, v in d.items():
            dropped[k] += v

    # depth support + depth-only CNVs
    depth = defaultdict(list)
    if a.exon_depth and os.path.exists(a.exon_depth):
        with open(a.exon_depth) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 7 and re.search("deletion|duplication", f[6]):
                    depth[f[3].split("_exon")[0]].append((norm_chrom(f[0], a.style_chr), int(f[1]) + 1, int(f[2]), f[3], f[6], f[5]))
    used = set()
    for x in rows:
        sup, ratios = [], []
        for (c, s, e, ex, flag, ratio) in depth.get(x["gene"], []):
            if c == x["chrom"] and s <= x["end"] and e >= x["start"] and \
               ((x["svtype"] in ("DEL", "CNV") and "deletion" in flag) or (x["svtype"] == "DUP" and "duplication" in flag)):
                sup.append(f"{ex}({ratio}x)")
                used.add((x["gene"], ex))
                if fnum(ratio) is not None:
                    ratios.append(fnum(ratio))
        if sup:
            x["notes"].append("depth-supported: " + ",".join(sup))
        # deletion genotype from depth when the SV file has none: ~0x both copies lost, ~0.5x one copy
        if x["svtype"] in ("DEL", "CNV") and ratios and x.get("gt_source") != "measured":
            x["notes"] = [n for n in x["notes"] if not n.startswith("genotype not measured")]
            if max(ratios) < 0.2:
                x["gt"], x["gt_source"] = "1/1", "depth"
                x["notes"].append(f"HOMOZYGOUS deletion by depth (max ratio {max(ratios):.2f}x)")
            elif max(ratios) < 0.65:
                x["gt"], x["gt_source"] = "0/1", "depth"
        if x["svlen"] and int(x["svlen"]) > 5_000_000:
            x["notes"].append("large CNV: also apply ACMG CNV scoring")
    for g, lst in depth.items():
        for kind in ("deletion", "duplication"):
            grp = sorted([y for y in lst if kind in y[4] and (g, y[3]) not in used], key=lambda y: y[1])
            if grp:
                exons = [y[3].split("_exon")[1] for y in grp]
                rows.append({"gene": g, "chrom": grp[0][0], "start": grp[0][1], "end": grp[-1][2],
                             "svtype": "DEL" if kind == "deletion" else "DUP", "svlen": grp[-1][2] - grp[0][1],
                             "gt": "1/1" if kind == "deletion" and all((fnum(y[5]) or 1) < 0.2 for y in grp) else "0/1",
                             "callers": "depth", "filter": "PASS", "qual": ".", "pr": "", "sr": "", "pop_af": "",
                             "class": "CNV_DEPTH", "score": CLASS["CNV_DEPTH"], "gt_source": "depth",
                             "effect": f"depth suggests {kind} of exon(s) {','.join(exons)} (no SV call)", "exons": ",".join(exons),
                             "notes": ["ratios " + ",".join(y[5] for y in grp), "confirm by MLPA/qPCR"], "ids": "depth"})

    if a.bam:
        n_dep = sv_depth(rows, a.bam, a.bam_ref, a.depth_min_len)
        n_gt = sum(1 for x in rows if x.get("gt_source") == "depth")
        print(f"[sv] depth from BAM: {n_dep:,} DEL/DUP intervals measured; {n_gt:,} SV-gene rows genotyped by depth "
              f"({sum(1 for x in rows if x.get('gt_source') == 'depth' and x['gt'] == '1/1'):,} homozygous); "
              f"{sum(1 for x in rows if 'NOT supported by depth' in ' '.join(x['notes'])):,} rows not supported by depth",
              file=sys.stderr)
        # weak effects (regulatory / intronic / dosage) of an SV the BAM does not support are dropped;
        # coding effects are kept (score -1, flagged) for review in IGV
        weak = ("SV_REGULATORY", "SV_INTRONIC", "SV_DOSAGE")
        bad = [x for x in rows if x["class"] in weak and "NOT supported by depth" in " ".join(x["notes"])]
        if bad:
            for x in bad:
                for i in str(x.get("ids", "")).split(","):
                    depth_dropped[(i.split("|")[0], x["gene"])] = x.get("depth_ratio", "-")
            ids = set(map(id, bad))
            rows = [x for x in rows if id(x) not in ids]
            dropped["weak-effect rows of SVs not supported by depth"] += len(bad)
    # every input SV within 50 kb of an MVA gene, and what happened to it
    if mva_fate:
        fp = os.path.join(os.path.dirname(os.path.abspath(a.out)), "mva_gene_svs_all.tsv")
        kept = defaultdict(list)
        for x in rows:
            for i in str(x.get("ids", "")).split(","):
                kept[(i.split("|")[0], x["gene"])].append(x)
        n_kept = 0
        with open(fp, "w") as fh:
            fh.write("AnnotSV_ID\tcaller_ID\tchrom\tstart\tend\tsvtype\tlength\tmva_gene\trelation\tin_ranking\t"
                     "class\tscore\tgt\tgt_source\tdepth_ratio\tfate\tAnnotSV_Gene_name\n")
            for sid, fx in sorted(mva_fate.items(), key=lambda kv: (kv[1]["chrom"], kv[1]["start"])):
                for g, rel in fx["near"]:
                    ks = kept.get((sid, g))
                    if ks:
                        x = ks[0]
                        n_kept += 1
                        row = ["yes", x["class"], x["score"], x["gt"], x.get("gt_source", "-"), x.get("depth_ratio", "-"),
                               "kept: " + x["effect"]]
                    else:
                        reason = fx["reason"]
                        if (sid, g) in depth_dropped:
                            reason = f"weak effect (regulatory/intronic) and NOT supported by read depth ({depth_dropped[(sid, g)]}x of flanks)"
                        row = ["no", "-", "-", "-", "-", "-", "removed: " + reason]
                    fh.write("\t".join(str(v) for v in [sid, fx["id"] or "-", fx["chrom"], fx["start"], fx["end"],
                             fx["svtype"], abs(fx["end"] - fx["start"]), g, rel] + row + [fx["gene_name"][:80] or "-"]) + "\n")
        print(f"[sv] MVA genes: {len(mva_fate):,} input SV(s) within 50 kb of an MVA gene; {n_kept:,} SV-gene pair(s) kept "
              f"for ranking; every one with its fate -> {fp}", file=sys.stderr)
    cols = ["gene", "chrom", "start", "end", "svtype", "svlen", "gt", "gt_source", "depth_ratio", "callers", "filter", "qual", "pr", "sr",
            "pop_af", "confidence", "class", "score", "effect", "exons", "acmg", "cytoband", "gene_info", "gene_flag",
            "notes", "ids"]
    defaults = {"depth_ratio": "-", "gt_source": "measured", "confidence": "-", "acmg": "-", "cytoband": "-", "gene_info": "-", "gene_flag": "-"}
    rows.sort(key=lambda r: (-r["score"], r["gene"]))
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            r = dict(defaults, **r)
            r["notes"] = "; ".join(r["notes"]) or "-"
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"[sv] {len(rows):,} SV-gene candidate rows -> {a.out}" +
          ("; excluded: " + "; ".join(f"{k}: {v:,}" for k, v in dropped.items()) if dropped else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
