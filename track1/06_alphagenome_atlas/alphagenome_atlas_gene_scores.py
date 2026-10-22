#!/usr/bin/env python3
"""Fetch pre-computed AlphaGenome Atlas scores around target genes and merge
them with variants of interest.

INPUTS
------
  --genes     TSV of target genes. The gene column is picked from (in order)
              --gene-col, a column named gene / gene_name / gene_symbol /
              symbol / hgnc_symbol / gene_id / ensembl_gene_id / target_gene /
              genes (case-insensitive), else the first column. Single-column
              files without a header are fine. Extra columns are carried into
              the interval output.
  --gtf       Reference GTF (GRCh38; .gtf or .gtf.gz), e.g. GENCODE/Ensembl.
  --variants  Either
                - a VCF (.vcf / .vcf.gz / .vcf.bgz), or
                - an annotated, prioritised variants TSV/CSV with at least
                  CHROM, POS, REF, ALT (e.g. the VEP/ClinVar/AlphaMissense
                  export: CHROM POS ID REF ALT QUAL FILTER ClinVar_* GT GQ DP
                  VAF AD Consequence IMPACT SYMBOL Gene ... priority_reason).
              Format is auto-detected; override with --variant-format.

OUTPUTS (<prefix> = --out-prefix)
---------------------------------
  <prefix>.gene_intervals.tsv       target genes, gene coordinates and the
                                    +/- flank interval of interest (1-based,
                                    inclusive), plus your extra gene columns
  <prefix>.gene_intervals.bed       the same intervals as BED (0-based)
  <prefix>.atlas_query_intervals.bed  the exact intervals requested from Atlas
  <prefix>.atlas_queried_positions.tsv  every queried position, the Atlas
                                    reference base and the alleles returned
  <prefix>.position_alleles.tsv     ALL possible alleles (3 per position) at
                                    every queried position, independent of
                                    your ALTs: AVI (avi_raw, avi_cdf_quantile,
                                    avi_phred, avi_top_percent), AVI feature
                                    attributions (avi_fi_*, avi_top_feature),
                                    per-scorer strongest track score; then
                                    is_input_variant / input_variant_ids /
                                    target_genes mark which alleles are yours
  <prefix>.variants_atlas_summary.tsv  your variants (all original columns),
                                    one row per variant x target gene:
                                      avi_* for your observed ALT, plus
                                        avi_max_phred_any_allele /
                                        avi_top_allele_any_allele
                                      <SCORER>__max_abs_score / __top_track
                                        for your observed ALT allele
                                      <SCORER>__max_abs_score_any_allele /
                                        __top_allele_any_allele across all 3
                                      atlas_ref, atlas_alleles,
                                      ref_matches_atlas, atlas_status
  <prefix>.atlas_scores_all/        (--fetch-mode full --keep-all only)
                                    unmerged scores for every position in
                                    the gene intervals, all 3 alleles:
                                    part-NNNNN.parquet (or .tsv.gz) with
                                    chrom, pos, atlas_ref, atlas_alt, scorer,
                                    atlas_gene_*, track_idx, raw_score,
                                    quantile_score.
  <prefix>.atlas_scores_all.tracks.tsv  maps (scorer, track_idx) -> track
                                    metadata for the dump. Read with:
                                    pd.read_parquet(dir).merge(
                                      pd.read_csv(tracks, sep="\t"),
                                      on=["scorer", "track_idx"])
  <prefix>.variants_atlas_long.parquet (or .tsv.gz; skip with --no-long)
                                    every retrieved score at your variant
                                    positions for ALL THREE alleles: one row
                                    per variant x target gene x atlas_alt x
                                    scorer x (gene) x track, written in
                                    chunks. Carries only identifying columns
                                    (variant_key, CHROM, POS, ID, REF, ALT,
                                    target_gene); join to the summary on
                                    variant_key for your annotations.
                                    is_observed_alt marks your ALT.
  <prefix>.atlas_fetch_cache/       raw fetched scores, saved during the
                                    fetch; --resume reuses them (no API calls)

STEPS
-----
1. Look up each target gene in the GTF ("gene" features; symbol or Ensembl ID).
2. Interval of interest = gene body +/- 1 Mb (--flank). Overlapping intervals
   are merged so each region is fetched once.
3. Query the Atlas BY POSITION ONLY (AtlasClient.query_interval on genomic
   intervals; your REF/ALT are never sent). Each position returns all three
   possible alternate alleles, for the AVI score (AVI_SCORE: one value per
   SNV; AVI_SCORE_FEATURE_IMPORTANCE: SHAP attributions) and for the
   per-track variant scorers. AVI scorers are always added unless --no-avi.
     --fetch-mode targeted (default): a 1-bp interval at each variant
       position (consecutive positions joined), fetched in parallel.
     --fetch-mode full: every position of every gene interval, fetched in
       --chunk-size pieces. Only your variant positions are merged; add
       --keep-all to ALSO write every position's unmerged scores to
       <prefix>.atlas_scores_all/ (very large: ~6 M SNVs per 2 Mb x
       thousands of tracks; narrow with --scorers / --ontology-terms).
4. Match to your variants by position, keep all three alleles, and flag the
   one matching your ALT (is_observed_alt / is_input_variant). REF is only checked
   (ref_matches_atlas), never used for matching.

NOTES
-----
* The Atlas holds single-nucleotide substitutions only. For indels/MNVs the
  three SNV alleles at the anchor position are still reported, but none
  matches ALT (atlas_status = not_snv_position_alleles_only).
* atlas_status: scored | outside_gene_intervals | not_snv_position_alleles_only
  | alt_not_in_atlas_alleles (ALT equals the reference base, e.g. swapped
  REF/ALT) | queried_no_atlas_result.
* Multi-allelic ALT ("A,T") is split into one row per allele.
* Chromosome naming ("1" vs "chr1", "MT" vs "chrM") is normalised.
* VCF input is streamed and only records inside the intervals are kept, so a
  whole-genome VCF is fine. For a TSV, all rows are kept; rows outside every
  interval are reported with atlas_status = outside_gene_intervals.
* VCF sample columns: with one sample, FORMAT fields become columns (GT, GQ,
  DP, AD, ...). With several samples, they become <sample>:<field>. INFO is
  kept as one raw column; add --expand-info to split it into INFO_<key>.

EXAMPLE
-------
  pip install alphagenome pandas pyarrow
  export ALPHAGENOME_API_KEY=...
  python alphagenome_atlas_gene_scores.py \
      --genes target_genes.tsv \
      --gtf gencode.v46.annotation.gtf.gz \
      --variants prioritised_variants.tsv \
      --out-prefix results/run1 \
      --scorers RNA_SEQ CAGE ATAC DNASE SPLICE_SITES SPLICE_SITE_USAGE
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import logging
import os
import re
import sys
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd

log = logging.getLogger("atlas_gene_scores")

TILE = 32  # Atlas sub-request size used by AtlasClient.query_interval
DEFAULT_FLANK = 1_000_000
CACHE_PART_ROWS = 2_000_000  # rows per saved cache part (bounds peak memory)
VARIANT_COLS = ["CHROM", "POS", "REF", "ALT"]
AVI_SCORER = "AVI_SCORE"  # 1 value per SNV: X = raw score, quantiles layer = CDF
AVI_FI_SCORER = "AVI_SCORE_FEATURE_IMPORTANCE"  # SHAP attribution per AVI feature
AVI_SCORERS = (AVI_SCORER, AVI_FI_SCORER)
GENE_COL_CANDIDATES = [
    "gene", "gene_name", "gene_symbol", "symbol", "hgnc_symbol",
    "gene_id", "ensembl_gene_id", "target_gene", "genes",
]


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def norm_chrom(c) -> str:
    """UCSC-style chromosome name ('chr1', 'chrX', 'chrM')."""
    c = str(c).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    if c == "MT":
        c = "M"
    return "chr" + c


def strip_version(gene_id: str) -> str:
    return re.sub(r"\.\d+$", "", str(gene_id))


def variant_key(chrom, pos, ref, alt) -> str:
    return f"{norm_chrom(chrom)}:{int(pos)}:{ref}>{alt}"


def pos_key(chrom, pos) -> str:
    return f"{norm_chrom(chrom)}:{int(pos)}"


def _open(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


# ----------------------------------------------------------------------------
# 1. target genes + GTF
# ----------------------------------------------------------------------------
def read_gene_table(path: str, gene_col: str | None = None) -> pd.DataFrame:
    """Return a DataFrame with a 'query' column (gene) + any extra columns."""
    sep = "," if re.search(r"\.csv(\.gz)?$", path) else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, comment="#")
    lower = {c.lower(): c for c in df.columns}
    if gene_col:
        if gene_col not in df.columns:
            sys.exit(f"--gene-col '{gene_col}' not in gene table columns {list(df.columns)}")
        col = gene_col
    else:
        col = next((lower[c] for c in GENE_COL_CANDIDATES if c in lower), None)
        if col is None and df.shape[1] == 1:
            # single column, no recognisable header -> treat as headerless
            df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, header=None, comment="#")
            df.columns = ["gene"]
            col = "gene"
        elif col is None:
            col = df.columns[0]
            log.warning("No standard gene column found; using first column '%s'", col)
    df = df.rename(columns={col: "query"})
    df["query"] = df["query"].str.strip()
    df = df[df["query"] != ""].drop_duplicates("query").reset_index(drop=True)
    log.info("Read %d target genes from %s (column '%s')", len(df), path, col)
    return df


@dataclass
class Gene:
    query: str
    gene_id: str
    gene_name: str
    gene_type: str
    chrom: str  # 'chrN'
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    strand: str


def _parse_gtf_attrs(attr: str) -> dict:
    return {m.group(1): m.group(2) for m in re.finditer(r'(\w+)\s+"([^"]*)"', attr)}


def _is_primary(chrom: str) -> bool:
    return re.fullmatch(r"chr(\d+|X|Y|M)", chrom) is not None


def genes_from_gtf(gtf_path: str, genes: list[str]) -> list[Gene]:
    """Find genes (symbol or Ensembl ID, case-insensitive, version-agnostic)."""
    lookup = {}
    for g in genes:
        lookup[g.upper()] = g
        lookup[strip_version(g).upper()] = g
    found: dict[str, Gene] = {}
    with _open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = _parse_gtf_attrs(f[8])
            name, gid = attrs.get("gene_name", ""), attrs.get("gene_id", "")
            hit = next(
                (lookup[k] for k in (name.upper(), gid.upper(), strip_version(gid).upper()) if k and k in lookup),
                None,
            )
            if hit is None:
                continue
            g = Gene(hit, gid, name or gid, attrs.get("gene_type", attrs.get("gene_biotype", "")),
                     norm_chrom(f[0]), int(f[3]), int(f[4]), f[6])
            prev = found.get(hit)
            # duplicates (PAR_Y copies, alt contigs): prefer primary chromosome, non-PAR
            if prev is None or (_is_primary(g.chrom) and not _is_primary(prev.chrom)) or (
                "_PAR_Y" in prev.gene_id and "_PAR_Y" not in g.gene_id
            ):
                found[hit] = g
    missing = [g for g in genes if g not in found]
    if missing:
        log.warning("%d gene(s) not found in GTF: %s", len(missing), ", ".join(missing))
    return [found[g] for g in genes if g in found]


def build_windows(genes: list[Gene], flank: int, gene_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for g in genes:
        ws = max(0, g.start - 1 - flank)  # 0-based half-open internally
        we = g.end + flank
        rows.append(dict(
            query=g.query, gene_id=g.gene_id, gene_name=g.gene_name, gene_type=g.gene_type,
            chrom=g.chrom, gene_start=g.start, gene_end=g.end, strand=g.strand,
            interval_start=ws + 1, interval_end=we, interval_length=we - ws,
            upstream_flank=flank, downstream_flank=flank,
            win_start=ws, win_end=we,
        ))
    win = pd.DataFrame(rows)
    extra = gene_table.drop(columns=[c for c in gene_table.columns if c in win.columns and c != "query"])
    return win.merge(extra, on="query", how="left")


def merge_windows(win: pd.DataFrame) -> list[tuple[str, int, int]]:
    """Merge overlapping per-gene windows into non-redundant regions."""
    regions = []
    for chrom, sub in win.groupby("chrom"):
        cur = None
        for s, e in sorted(zip(sub.win_start, sub.win_end)):
            if cur and s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                if cur:
                    regions.append((chrom, *cur))
                cur = [s, e]
        regions.append((chrom, *cur))
    return regions


class RegionIndex:
    """Fast 'is this position inside any region' lookup."""

    def __init__(self, regions):
        self.by_chrom = {}
        for c, s, e in regions:
            self.by_chrom.setdefault(c, []).append((s, e))
        for c in self.by_chrom:
            self.by_chrom[c].sort()
        self.starts = {c: [s for s, _ in v] for c, v in self.by_chrom.items()}

    def contains(self, chrom: str, pos1: int) -> bool:
        iv = self.by_chrom.get(chrom)
        if not iv:
            return False
        i = bisect.bisect_right(self.starts[chrom], pos1 - 1) - 1
        return i >= 0 and iv[i][0] <= pos1 - 1 < iv[i][1]


# ----------------------------------------------------------------------------
# 2. variants (VCF or TSV)
# ----------------------------------------------------------------------------
def detect_variant_format(path: str) -> str:
    if re.search(r"\.vcf(\.b?gz)?$", path):
        return "vcf"
    with _open(path) as fh:
        first = fh.readline()
    return "vcf" if first.startswith("##fileformat=VCF") else "tsv"


def read_vcf(path: str, region_index: RegionIndex | None, expand_info: bool) -> pd.DataFrame:
    samples, rows, n_total = [], [], 0
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line[1:].rstrip("\n").split("\t")[9:]
                continue
            n_total += 1
            f = line.rstrip("\n").split("\t")
            if region_index is not None and not region_index.contains(norm_chrom(f[0]), int(f[1])):
                continue
            rec = dict(CHROM=f[0], POS=f[1], ID=f[2], REF=f[3], ALT=f[4],
                       QUAL=f[5] if len(f) > 5 else "", FILTER=f[6] if len(f) > 6 else "",
                       INFO=f[7] if len(f) > 7 else "")
            if expand_info and rec["INFO"] not in ("", "."):
                for kv in rec["INFO"].split(";"):
                    k, _, v = kv.partition("=")
                    rec[f"INFO_{k}"] = v if _ else "True"
            if len(f) > 9:
                fmt = f[8].split(":")
                for s, val in zip(samples, f[9:]):
                    for k, v in zip(fmt, val.split(":")):
                        rec[k if len(samples) == 1 else f"{s}:{k}"] = v
            rows.append(rec)
    log.info("VCF: kept %d / %d records inside the gene intervals", len(rows), n_total)
    return pd.DataFrame(rows, dtype=str)


def read_variant_table(path: str) -> pd.DataFrame:
    sep = "," if re.search(r"\.csv(\.gz)?$", path) else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def prepare_variants(df: pd.DataFrame) -> pd.DataFrame:
    miss = [c for c in VARIANT_COLS if c not in df.columns]
    if miss:
        sys.exit(f"Variant table is missing required columns: {miss}")
    df = df.copy()
    df["ALT"] = df["ALT"].astype(str).str.split(",")
    df = df.explode("ALT", ignore_index=True)
    df["POS"] = df["POS"].astype(int)
    df["_chrom_norm"] = df["CHROM"].map(norm_chrom)
    df["pos_key"] = [pos_key(c, p) for c, p in zip(df["_chrom_norm"], df["POS"])]
    df["variant_key"] = [variant_key(c, p, r, a) for c, p, r, a in
                         zip(df["_chrom_norm"], df["POS"], df["REF"], df["ALT"])]
    acgt = list("ACGT")
    df["is_snv"] = df["REF"].isin(acgt) & df["ALT"].isin(acgt) & (df["REF"] != df["ALT"])
    return df


def assign_target_genes(variants: pd.DataFrame, win: pd.DataFrame) -> pd.DataFrame:
    """Variant x gene-interval pairs (a variant can fall in several intervals)."""
    parts = []
    for _, g in win.iterrows():
        m = variants[(variants["_chrom_norm"] == g.chrom)
                     & (variants["POS"] - 1 >= g.win_start) & (variants["POS"] - 1 < g.win_end)]
        if m.empty:
            continue
        m = m.copy()
        m["target_gene"] = g.gene_name
        m["target_gene_id"] = g.gene_id
        # signed distance relative to transcription direction: <0 upstream, >0 downstream
        before = m["POS"] < g.gene_start
        after = m["POS"] > g.gene_end
        dist = np.where(before, m["POS"] - g.gene_start, np.where(after, m["POS"] - g.gene_end, 0))
        if g.strand == "-":
            dist = -dist
        m["dist_to_gene_bp"] = dist
        m["position_vs_gene"] = np.where(dist < 0, "upstream", np.where(dist > 0, "downstream", "gene_body"))
        parts.append(m)
    cols = ["target_gene", "target_gene_id", "dist_to_gene_bp", "position_vs_gene"]
    if not parts:
        return variants.iloc[0:0].assign(**{c: None for c in cols})
    return pd.concat(parts, ignore_index=True)


# ----------------------------------------------------------------------------
# 3. Atlas fetching (by position only; all three alternate alleles returned)
# ----------------------------------------------------------------------------
def targeted_intervals(chrom: str, positions: np.ndarray):
    """One interval per variant position (0-based [pos-1, pos)). Runs of
    consecutive positions are joined, up to TILE bp, so each request still
    maps to a single Atlas sub-request."""
    runs = []
    for p in np.unique(positions.astype(int)):
        if runs and p - 1 == runs[-1][1] and (p - runs[-1][0]) <= TILE:
            runs[-1][1] = p
        else:
            runs.append([p - 1, p])
    return [(chrom, int(s), int(e)) for s, e in runs]


def full_intervals(chrom: str, start: int, end: int, step: int = 1_000):
    return [(chrom, s, min(s + step, end)) for s in range(start, end, step)]


TRACK_META_COLS = ("name", "strand", "ontology_curie", "biosample_name", "biosample_type",
                   "gtex_tissue", "transcription_factor", "histone_mark", "data_source", "Assay title")


def adata_to_long(adatas: dict, keep_pos: set[str] | None, gene_filter: set[str] | None):
    """Flatten {scorer: AnnData} into a long DataFrame keyed by position,
    keeping every allele the Atlas returns. Filters early to bound memory.

    Track metadata is NOT repeated on every row: rows carry (scorer,
    track_idx), and the per-scorer track table is returned separately.
    Returns (long_df, {scorer: track_metadata_df})."""
    frames, track_meta = [], {}
    keep_arr = np.array(sorted(keep_pos)) if keep_pos is not None else None
    for scorer, ad in adatas.items():
        obs = ad.obs.reset_index(drop=True)
        if "variant" not in obs.columns:
            continue
        var = ad.var.reset_index(drop=True)
        vmeta = var[[c for c in TRACK_META_COLS if c in var.columns]].add_prefix("track_")
        vmeta.insert(0, "track_idx", np.arange(len(var)))
        track_meta[scorer] = vmeta
        vs = obs["variant"].tolist()
        pkeys = np.array([pos_key(v.chromosome, v.position) for v in vs])
        refs = np.array([v.reference_bases for v in vs])
        alts = np.array([v.alternate_bases for v in vs])
        mask = np.ones(len(obs), dtype=bool)
        if keep_arr is not None:
            mask &= np.isin(pkeys, keep_arr)
        if gene_filter is not None and ("gene_name" in obs or "gene_id" in obs):
            gn = obs["gene_name"].astype(str).str.upper() if "gene_name" in obs else pd.Series("", index=obs.index)
            gi = obs["gene_id"].astype(str).map(strip_version).str.upper() if "gene_id" in obs else pd.Series("", index=obs.index)
            mask &= (gn.isin(gene_filter) | gi.isin(gene_filter)).to_numpy()
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            continue
        X = np.asarray(ad.X)[idx]
        Q = np.asarray(ad.layers["quantiles"])[idx] if "quantiles" in ad.layers else None
        n_tracks = X.shape[1]
        base = pd.DataFrame({
            "pos_key": np.repeat(pkeys[idx], n_tracks),
            "atlas_ref": np.repeat(refs[idx], n_tracks),
            "atlas_alt": np.repeat(alts[idx], n_tracks),
            "scorer": scorer,
            "atlas_gene_id": np.repeat(obs["gene_id"].to_numpy()[idx], n_tracks) if "gene_id" in obs else None,
            "atlas_gene_name": np.repeat(obs["gene_name"].to_numpy()[idx], n_tracks) if "gene_name" in obs else None,
            "track_idx": np.tile(np.arange(n_tracks, dtype=np.int32), len(idx)),
            "raw_score": X.ravel().astype(np.float32),
            "quantile_score": Q.ravel().astype(np.float32) if Q is not None else np.float32(np.nan),
        })
        frames.append(base)
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), track_meta


def track_table(track_meta: dict) -> pd.DataFrame:
    return pd.concat([m.assign(scorer=k) for k, m in track_meta.items()], ignore_index=True) \
        if track_meta else pd.DataFrame(columns=["scorer", "track_idx"])


def _is_transient(err: Exception) -> bool:
    if isinstance(err, (TimeoutError, ConnectionError)):
        return True
    code = getattr(err, "code", None)
    if callable(code):
        try:
            import grpc

            return code() in (grpc.StatusCode.UNAVAILABLE, grpc.StatusCode.DEADLINE_EXCEEDED,
                              grpc.StatusCode.RESOURCE_EXHAUSTED) and "larger than max" not in str(err)
        except Exception:
            return False
    return False


def fetch_interval(client, interval, scorers, ontology_terms, keep_pos, gene_filter, inner_workers):
    from alphagenome.data import genome

    chrom, s, e = interval
    iv = genome.Interval(chrom, s, e)
    for attempt in range(5):
        try:
            res = client.query_interval(iv, requested_scorers=scorers, ontology_terms=ontology_terms,
                                        progress_bar=False, max_workers=inner_workers)
            return adata_to_long(res, keep_pos, gene_filter)
        except Exception as err:
            if not _is_transient(err) or attempt == 4:
                raise RuntimeError(f"Failed fetching {chrom}:{s}-{e}: {err}") from err
            time.sleep(2 ** attempt)


def get_client(api_key: str, timeout: float = 60):
    """Same as alphagenome.atlas.atlas.create(), but with no gRPC message-size
    cap. atlas.create() keeps gRPC's 4 MB default, and Atlas responses are
    often tens of MB ("Received message larger than max"). The DNA-model client
    already sets these limits to -1 (unlimited)."""
    import importlib.resources

    import grpc
    from alphagenome.atlas import atlas
    from alphagenome.protos import atlas_service_pb2_grpc

    service_config = (importlib.resources.files("alphagenome") / "protos/grpc_service_config.json").read_text()
    channel = grpc.secure_channel(
        "dns:///gdmscience.googleapis.com:443",
        grpc.ssl_channel_credentials(),
        options=(
            ("grpc.service_config", service_config),
            ("grpc.max_send_message_length", -1),
            ("grpc.max_receive_message_length", -1),
        ),
    )
    grpc.channel_ready_future(channel).result(timeout)
    stub = atlas_service_pb2_grpc.AtlasServiceStub(channel=channel)
    return atlas.AtlasClient(stub, metadata=[("x-goog-api-key", api_key)])


# ----------------------------------------------------------------------------
# 4. merge (by position, then ALT) + outputs
# ----------------------------------------------------------------------------
def _top_by_abs(df: pd.DataFrame, keys: list[str], cols: list[str]) -> pd.DataFrame:
    idx = df.groupby(keys, dropna=False, sort=False)["_abs"].idxmax().dropna()
    return df.loc[idx, keys + cols]


def compute_top(scores: pd.DataFrame, tracks: pd.DataFrame) -> pd.DataFrame:
    """Strongest track per position x allele x scorer x gene.

    This is the only table the summaries are built from, so memory stays small
    even when `scores` has tens of millions of rows."""
    keys = ["pos_key", "atlas_ref", "atlas_alt", "scorer", "atlas_gene_name", "atlas_gene_id"]
    metric = scores["quantile_score"].where(scores["quantile_score"].notna(), scores["raw_score"])
    idx = metric.abs().groupby([scores[k] for k in keys], dropna=False, sort=False).idxmax()
    idx = idx.dropna().astype(np.int64).to_numpy()
    top = scores.loc[idx, keys + ["track_idx", "raw_score", "quantile_score"]].reset_index(drop=True)
    top["metric"] = metric.to_numpy()[idx]
    return top.merge(tracks, on=["scorer", "track_idx"], how="left")


def summarize_variants(vg: pd.DataFrame, top: pd.DataFrame, all_genes: bool) -> pd.DataFrame:
    """One row per variant x target gene: your original columns + per-scorer
    strongest score for your ALT and for any of the three alleles."""
    base = vg.drop(columns=["_chrom_norm"]).copy()
    if top.empty:
        return base.assign(atlas_available=False)

    pos_info = top.groupby("pos_key").agg(atlas_ref=("atlas_ref", "first"),
                                          atlas_alleles=("atlas_alt", lambda x: ",".join(sorted(set(x)))))
    base = base.merge(pos_info, on="pos_key", how="left")
    base["ref_matches_atlas"] = (base["REF"].astype(str).str[:1] == base["atlas_ref"]).where(base["atlas_ref"].notna())
    have = set(top["pos_key"] + ">" + top["atlas_alt"])
    base["atlas_available"] = (base["pos_key"] + ">" + base["ALT"].astype(str)).isin(have)

    tv = top[~top["scorer"].isin(AVI_SCORERS)]
    m = vg[["variant_key", "pos_key", "ALT", "target_gene", "target_gene_id"]].merge(tv, on="pos_key")
    if not all_genes:
        has_gene = m["atlas_gene_name"].notna() | m["atlas_gene_id"].notna()
        same = (m["atlas_gene_name"].astype(str).str.upper() == m["target_gene"].astype(str).str.upper()) | (
            m["atlas_gene_id"].astype(str).map(strip_version) == m["target_gene_id"].astype(str).map(strip_version))
        m = m[~has_gene | same]
    if m.empty:
        return base
    m = m.assign(_abs=m["metric"].abs())
    keys = ["variant_key", "target_gene", "scorer"]
    obs = _top_by_abs(m[m["atlas_alt"] == m["ALT"]], keys, ["metric", "track_name"])
    anya = _top_by_abs(m, keys, ["metric", "track_name", "atlas_alt"])

    def wide(df, col, suffix):
        w = df.pivot(index=["variant_key", "target_gene"], columns="scorer", values=col)
        w.columns = [f"{c}__{suffix}" for c in w.columns]
        return w.reset_index()

    for df, col, suf in ((obs, "metric", "max_abs_score"), (obs, "track_name", "top_track"),
                         (anya, "metric", "max_abs_score_any_allele"),
                         (anya, "atlas_alt", "top_allele_any_allele")):
        if len(df):
            base = base.merge(wide(df, col, suf), on=["variant_key", "target_gene"], how="left")
    return base


def build_position_alleles(fi_rows: pd.DataFrame, top: pd.DataFrame, tracks: pd.DataFrame) -> pd.DataFrame:
    """One row per queried position x Atlas allele (all three), independent of
    the input variants: AVI raw / CDF / PHRED, AVI feature attributions, and
    per-scorer strongest track score for that allele."""
    keys = ["pos_key", "atlas_ref", "atlas_alt"]
    if top.empty:
        return pd.DataFrame(columns=keys)
    tbl = top[keys].drop_duplicates()

    avi = top[top["scorer"] == AVI_SCORER]  # single track -> top row is the value
    if len(avi):
        a = avi.groupby(keys).agg(avi_raw=("raw_score", "first"),
                                  avi_cdf_quantile=("quantile_score", "first")).reset_index()
        a["avi_tail_quantile"] = 1 - a["avi_cdf_quantile"]
        a["avi_phred"] = -10 * np.log10(a["avi_tail_quantile"].clip(lower=1e-12))  # 20 = top 1 %, 30 = top 0.1 %
        a["avi_top_percent"] = 100 * a["avi_tail_quantile"]
        tbl = tbl.merge(a, on=keys, how="left")

    fi = fi_rows
    if len(fi):
        fi = fi.merge(tracks, on=["scorer", "track_idx"], how="left")
        feat = fi["track_name"].astype(str) if "track_name" in fi else fi["track_idx"].astype(str)
        w = fi.assign(_feat=feat).pivot_table(index=keys, columns="_feat", values="raw_score", aggfunc="first")
        topf = w.abs().idxmax(axis=1)
        w.columns = [f"avi_fi_{c}" for c in w.columns]
        w["avi_top_feature"] = topf
        w["avi_top_feature_value"] = [w.loc[i, f"avi_fi_{t}"] for i, t in zip(w.index, topf)]
        tbl = tbl.merge(w.reset_index(), on=keys, how="left")

    tr = top[~top["scorer"].isin(AVI_SCORERS)]
    if len(tr):
        tr = _top_by_abs(tr.assign(_abs=tr["metric"].abs()), keys + ["scorer"], ["metric", "track_name"])
        for col, suf in (("metric", "max_abs_score"), ("track_name", "top_track")):
            w = tr.pivot(index=keys, columns="scorer", values=col)
            w.columns = [f"{c}__{suf}" for c in w.columns]
            tbl = tbl.merge(w.reset_index(), on=keys, how="left")

    split = tbl["pos_key"].str.rsplit(":", n=1, expand=True)
    tbl.insert(0, "chrom", split[0])
    tbl.insert(1, "pos", split[1].astype(int))
    return tbl.sort_values(["chrom", "pos", "atlas_alt"]).reset_index(drop=True)


def attach_avi(summary: pd.DataFrame, pa: pd.DataFrame) -> pd.DataFrame:
    """Add AVI for the observed ALT and the strongest AVI among all 3 alleles."""
    if "avi_raw" not in pa.columns:
        return summary
    avi_cols = [c for c in pa.columns if c.startswith("avi_")]
    obs = pa[["pos_key", "atlas_alt"] + avi_cols].rename(columns={"atlas_alt": "ALT"})
    summary = summary.merge(obs, on=["pos_key", "ALT"], how="left")
    best = pa.loc[pa.groupby("pos_key")["avi_phred"].idxmax().dropna(), ["pos_key", "atlas_alt", "avi_phred"]]
    best = best.rename(columns={"atlas_alt": "avi_top_allele_any_allele", "avi_phred": "avi_max_phred_any_allele"})
    return summary.merge(best, on="pos_key", how="left")


LONG_ID_COLS = ["variant_key", "CHROM", "POS", "ID", "REF", "ALT", "target_gene", "target_gene_id"]


class LongWriter:
    """Appends every score at your variant positions (all 3 alleles) to the
    long table, one cache part at a time, carrying only identifying variant
    columns (join to the summary on variant_key for annotations)."""

    NUMERIC = {"POS", "raw_score", "quantile_score", "track_idx", "is_observed_alt", "ref_matches_atlas"}

    def __init__(self, prefix: str, fmt: str, vg: pd.DataFrame, tracks: pd.DataFrame, all_genes: bool,
                 max_rows: int = 2_000_000):
        id_cols = [c for c in LONG_ID_COLS if c in vg.columns]
        self.vmin, self.tracks, self.all_genes, self.max_rows = vg[id_cols + ["pos_key"]], tracks, all_genes, max_rows
        self.fmt, self.writer, self.schema, self.n, self.first = fmt, None, None, 0, True
        if fmt == "parquet":
            try:
                import pyarrow.parquet  # noqa: F401
            except ImportError:
                log.warning("pyarrow not installed; writing tsv.gz instead")
                self.fmt = "tsv"
        self.path = f"{prefix}.variants_atlas_long." + ("parquet" if self.fmt == "parquet" else "tsv.gz")
        if os.path.exists(self.path):
            os.remove(self.path)

    def write(self, scores: pd.DataFrame):
        codes, uniq = pd.factorize(scores["pos_key"], sort=True)
        order = np.argsort(codes, kind="stable")
        bounds = np.searchsorted(codes[order], np.arange(len(uniq) + 1))
        counts = np.diff(bounds)
        i = 0
        while i < len(uniq):
            j, rows = i, 0
            while j < len(uniq) and (rows == 0 or rows + counts[j] <= self.max_rows):
                rows += counts[j]
                j += 1
            self._write_chunk(scores.iloc[order[bounds[i]:bounds[j]]])
            i = j

    def _write_chunk(self, ch: pd.DataFrame):
        m = self.vmin.merge(ch, on="pos_key")
        if not self.all_genes:
            has_gene = m["atlas_gene_name"].notna() | m["atlas_gene_id"].notna()
            same = (m["atlas_gene_name"].astype(str).str.upper() == m["target_gene"].astype(str).str.upper()) | (
                m["atlas_gene_id"].astype(str).map(strip_version) == m["target_gene_id"].astype(str).map(strip_version))
            m = m[~has_gene | same]
        if m.empty:
            return
        m = m.assign(is_observed_alt=(m["atlas_alt"] == m["ALT"]),
                     ref_matches_atlas=(m["REF"].astype(str).str[:1] == m["atlas_ref"]))
        m = m.merge(self.tracks, on=["scorer", "track_idx"], how="left").drop(columns=["pos_key"])
        for c in m.columns:
            if c not in self.NUMERIC:
                m[c] = m[c].astype("string")
        if self.fmt == "parquet":
            import pyarrow as pa_
            import pyarrow.parquet as pq
            t = pa_.Table.from_pandas(m, preserve_index=False, schema=self.schema)
            if self.writer is None:
                self.schema = t.schema
                self.writer = pq.ParquetWriter(self.path, self.schema)
            self.writer.write_table(t)
        else:
            m.to_csv(self.path, sep="\t", index=False, mode="a", header=self.first, compression="gzip")
        self.first = False
        self.n += len(m)

    def close(self):
        if self.writer is not None:
            self.writer.close()


# ----------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True, help="TSV with target genes")
    p.add_argument("--gene-col", help="Column holding gene symbols/IDs (default: auto-detect)")
    p.add_argument("--gtf", required=True, help="Reference GTF (GRCh38), .gtf or .gtf.gz")
    p.add_argument("--variants", required=True, help="VCF (.vcf/.vcf.gz) or annotated variants TSV/CSV")
    p.add_argument("--variant-format", choices=["auto", "vcf", "tsv"], default="auto")
    p.add_argument("--expand-info", action="store_true", help="VCF: split INFO into INFO_<key> columns")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--flank", type=int, default=DEFAULT_FLANK, help="bp up/downstream of gene (default 1,000,000)")
    p.add_argument("--scorers", nargs="+", default=None,
                   help="Atlas scorers, e.g. RNA_SEQ CAGE ATAC DNASE CHIP_TF CHIP_HISTONE SPLICE_SITES "
                        "SPLICE_SITE_USAGE SPLICE_JUNCTIONS POLYADENYLATION CONTACT_MAPS PROCAP. "
                        "Default: every scorer the Atlas lists.")
    p.add_argument("--ontology-terms", nargs="+", default=None,
                   help="Restrict tracks to these ontology CURIEs, e.g. UBERON:0000955 CL:0000236")
    p.add_argument("--fetch-mode", choices=["targeted", "full"], default="targeted",
                   help="targeted (default): only the positions of your variants. "
                        "full: every position of every gene interval.")
    p.add_argument("--keep-all", action="store_true",
                   help="(full mode) also write the unmerged scores for EVERY position in the gene "
                        "intervals to <prefix>.atlas_scores_all/ (chunked part files). Very large.")
    p.add_argument("--chunk-size", type=int, default=1_000,
                   help="(full mode) bp per request chunk; bounds memory (default 1000)")
    p.add_argument("--resume", action="store_true",
                   help="Reuse scores saved in <prefix>.atlas_fetch_cache/ by a previous run with the same "
                        "settings, instead of calling the API again")
    p.add_argument("--no-long", action="store_true",
                   help="Skip the long per-track table (fastest; summary/position files still written)")
    p.add_argument("--no-avi", action="store_true",
                   help="Do not add the AVI_SCORE / AVI_SCORE_FEATURE_IMPORTANCE scorers")
    p.add_argument("--all-genes", action="store_true",
                   help="Keep gene-level scores for every gene in the interval, not only the target gene")
    p.add_argument("--max-workers", type=int, default=10, help="parallel Atlas requests")
    p.add_argument("--format", choices=["parquet", "tsv"], default="parquet", help="format of the long output")
    p.add_argument("--api-key", default=os.environ.get("AG_API"))
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv=None, client=None):
    import concurrent.futures as cf
    import functools

    a = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)

    # --- genes -> intervals ------------------------------------------------
    gene_table = read_gene_table(a.genes, a.gene_col)
    genes = genes_from_gtf(a.gtf, gene_table["query"].tolist())
    if not genes:
        sys.exit("None of the target genes were found in the GTF.")
    win = build_windows(genes, a.flank, gene_table)
    win.drop(columns=["win_start", "win_end"]).to_csv(f"{a.out_prefix}.gene_intervals.tsv", sep="\t", index=False)
    win.assign(score=0)[["chrom", "win_start", "win_end", "gene_name", "score", "strand"]].to_csv(
        f"{a.out_prefix}.gene_intervals.bed", sep="\t", index=False, header=False)
    log.info("Resolved %d/%d genes -> %s.gene_intervals.{tsv,bed}", len(win), len(gene_table), a.out_prefix)
    regions = merge_windows(win)

    # --- variants ------------------------------------------------------------
    fmt = detect_variant_format(a.variants) if a.variant_format == "auto" else a.variant_format
    log.info("Variant input format: %s", fmt)
    raw = read_vcf(a.variants, RegionIndex(regions), a.expand_info) if fmt == "vcf" else read_variant_table(a.variants)
    if raw.empty:
        sys.exit("No variants to score (none read, or none inside the gene intervals).")
    variants = prepare_variants(raw)
    vg = assign_target_genes(variants, win)
    log.info("%d variant rows (%d SNVs); %d variant x gene-interval pairs",
             len(variants), int(variants.is_snv.sum()), len(vg))

    # --- Atlas intervals (positions only; REF/ALT are not sent) --------------
    in_win = vg.drop_duplicates("pos_key")
    intervals = []
    for chrom, s, e in regions:
        if a.fetch_mode == "full":
            intervals += full_intervals(chrom, s, e, a.chunk_size)
        else:
            sub = in_win[(in_win["_chrom_norm"] == chrom) & (in_win["POS"] - 1 >= s) & (in_win["POS"] - 1 < e)]
            if len(sub):
                intervals += targeted_intervals(chrom, sub["POS"].to_numpy())
    pd.DataFrame(intervals, columns=["chrom", "start", "end"]).to_csv(
        f"{a.out_prefix}.atlas_query_intervals.bed", sep="\t", index=False, header=False)
    log.info("%d unique variant positions in %d merged region(s) -> %d Atlas request interval(s)",
             len(in_win), len(regions), len(intervals))

    if client is None:
        if not a.api_key:
            sys.exit("Set ALPHAGENOME_API_KEY or pass --api-key")
        client = get_client(a.api_key)
    # query_interval() re-fetches scorer metadata on every call; cache it once.
    client.scorer_metadata = functools.lru_cache(maxsize=None)(client.scorer_metadata)
    available = client.scorer_metadata()
    scorers = list(a.scorers or sorted(available.keys()))
    if not a.no_avi:
        for s_ in AVI_SCORERS:
            if s_ not in scorers:
                if s_ in available:
                    scorers.append(s_)
                else:
                    log.warning("%s is not offered by this Atlas endpoint; AVI columns will be missing", s_)
    unknown = [s_ for s_ in scorers if s_ not in available]
    if unknown:
        log.warning("Scorers not listed by the Atlas (will return nothing): %s", ", ".join(unknown))
    log.info("Scorers: %s", ", ".join(scorers))

    var_pos = set(in_win["pos_key"])
    n_tr = sum(len(available[s_].track_metadata) for s_ in scorers if s_ in available)
    est = len(var_pos) * 3 * n_tr
    log.info("Expected ~%s score rows at your %s variant positions (3 alleles x ~%s tracks; gene-level "
             "scorers add one set per target gene)", f"{est:,}", f"{len(var_pos):,}", f"{n_tr:,}")
    if est > 200_000_000:
        log.warning("That is a lot of data. It is processed in parts and will not exhaust memory, but the "
                    "download and the long table will be large. Consider --scorers / --ontology-terms, "
                    "or --no-long if the summary and position files are enough.")
    dump = a.fetch_mode == "full" and a.keep_all
    if a.keep_all and a.fetch_mode != "full":
        log.warning("--keep-all only applies to --fetch-mode full; ignoring")
    keep_pos = None if dump else var_pos
    dump_dir = None
    if dump:
        dump_dir = f"{a.out_prefix}.atlas_scores_all"
        os.makedirs(dump_dir, exist_ok=True)
        n_pos = sum(e - s for _, s, e in intervals)
        n_tracks = sum(len(available[s_].track_metadata) for s_ in scorers if s_ in available)
        log.warning("--keep-all: dumping all %s positions x 3 alleles x ~%s tracks (~%.1e rows, more for "
                    "gene-level scorers) to %s/. Narrow with --scorers/--ontology-terms if this is too big.",
                    f"{n_pos:,}", f"{n_tracks:,}", n_pos * 3 * n_tracks, dump_dir)
    gene_filter = None if a.all_genes else (
        {x.upper() for x in win["gene_name"]} | {strip_version(x).upper() for x in win["gene_id"]})

    # --- fetch (or resume from cache) ------------------------------------------
    # Fetched scores at your variant positions are saved to a cache folder as
    # they arrive, so a crash later never costs a re-download (--resume).
    import hashlib
    import json
    import pickle
    import shutil

    cache_dir = f"{a.out_prefix}.atlas_fetch_cache"
    fingerprint = hashlib.sha1(json.dumps({
        "scorers": sorted(scorers), "ontology": sorted(a.ontology_terms or []), "mode": a.fetch_mode,
        "all_genes": a.all_genes, "positions": sorted(var_pos), "keep_all": dump,
        "chunk": a.chunk_size if a.fetch_mode == "full" else None,
    }).encode()).hexdigest()
    marker = os.path.join(cache_dir, "_COMPLETE")
    resumed = False
    if a.resume and os.path.exists(marker):
        if open(marker).read().strip() == fingerprint:
            resumed = True
            log.info("--resume: loading fetched scores from %s/ (no API calls)", cache_dir)
        else:
            log.warning("--resume: cache in %s/ was made with different settings; fetching again", cache_dir)

    if not resumed:
        if os.path.isdir(cache_dir):
            shutil.rmtree(cache_dir)
        os.makedirs(cache_dir)
        # targeted: many tiny requests -> parallelise across positions.
        # full: few large intervals -> let query_interval parallelise inside each.
        outer, inner = (a.max_workers, 1) if a.fetch_mode == "targeted" else (1, a.max_workers)
        buf, buf_rows, track_meta, n_dumped, part, cpart = [], 0, {}, 0, 0, 0

        def flush():
            nonlocal buf, buf_rows, cpart
            if buf:
                with open(os.path.join(cache_dir, f"part-{cpart:05d}.pkl"), "wb") as fh:
                    pickle.dump(pd.concat(buf, ignore_index=True), fh, protocol=pickle.HIGHEST_PROTOCOL)
                cpart += 1
            buf, buf_rows = [], 0

        def completed(ex):
            """Yield results as they finish, keeping at most a few requests
            queued and never holding on to finished ones (bounded memory)."""
            it, pending = iter(intervals), set()
            while True:
                while len(pending) < outer * 4:
                    iv = next(it, None)
                    if iv is None:
                        break
                    pending.add(ex.submit(fetch_interval, client, iv, scorers, a.ontology_terms,
                                          keep_pos, gene_filter, inner))
                if not pending:
                    return
                done, pending = cf.wait(pending, return_when=cf.FIRST_COMPLETED)
                for fut in done:
                    yield fut.result()
                del done

        with cf.ThreadPoolExecutor(max_workers=outer) as ex:
            for n, (df, meta) in enumerate(completed(ex), 1):
                track_meta.update(meta)
                if not df.empty and dump_dir:
                    split = df["pos_key"].str.rsplit(":", n=1, expand=True)
                    out = df.assign(chrom=split[0], pos=split[1].astype(int))
                    out = out[["chrom", "pos"] + [c for c in out.columns if c not in ("chrom", "pos")]]
                    for c in ("atlas_gene_id", "atlas_gene_name"):  # fixed schema across parts
                        out[c] = out[c].astype("string")
                    if a.format == "parquet":
                        out.to_parquet(f"{dump_dir}/part-{part:05d}.parquet", index=False)
                    else:
                        out.to_csv(f"{dump_dir}/part-{part:05d}.tsv.gz", sep="\t", index=False)
                    part += 1
                    n_dumped += len(out)
                    df = df[df["pos_key"].isin(var_pos)]  # only your positions go on to the merge
                if not df.empty:
                    buf.append(df)
                    buf_rows += len(df)
                    if buf_rows >= CACHE_PART_ROWS:
                        flush()
                if n % 25 == 0 or n == len(intervals):
                    log.info("fetched %d/%d intervals", n, len(intervals))
        flush()
        tracks = track_table(track_meta)
        tracks.to_csv(os.path.join(cache_dir, "tracks.tsv"), sep="\t", index=False)
        with open(marker, "w") as fh:
            fh.write(fingerprint)
        if dump_dir:
            tracks.to_csv(f"{dump_dir}.tracks.tsv", sep="\t", index=False)
            log.info("Wrote %s rows in %d part file(s) to %s/ (+ %s.tracks.tsv)",
                     f"{n_dumped:,}", part, dump_dir, dump_dir)
        log.info("Saved fetched scores to %s/ (rerun with --resume to skip fetching)", cache_dir)

    # --- process the saved scores one part at a time ---------------------------
    # Each part holds whole positions (intervals never overlap), so per-part
    # results are exact. Memory stays at ~one part, however many variants.
    tracks = pd.read_csv(os.path.join(cache_dir, "tracks.tsv"), sep="\t", dtype={"scorer": str})
    parts = sorted(f for f in os.listdir(cache_dir) if f.endswith(".pkl"))
    long_writer = None if a.no_long else LongWriter(a.out_prefix, a.format, vg, tracks, a.all_genes)
    tops, fis, n_rows = [], [], 0
    fi_cols = ["pos_key", "atlas_ref", "atlas_alt", "scorer", "track_idx", "raw_score"]
    for k, f in enumerate(parts, 1):
        with open(os.path.join(cache_dir, f), "rb") as fh:
            chunk = pickle.load(fh)
        tops.append(compute_top(chunk, tracks))
        fis.append(chunk.loc[chunk["scorer"] == AVI_FI_SCORER, fi_cols])
        if long_writer is not None:
            long_writer.write(chunk)
        n_rows += len(chunk)
        del chunk
        log.info("processed part %d/%d (%s score rows so far)", k, len(parts), f"{n_rows:,}")
    if long_writer is not None:
        long_writer.close()
    top = pd.concat(tops, ignore_index=True) if tops else pd.DataFrame(
        columns=["pos_key", "atlas_ref", "atlas_alt", "scorer", "atlas_gene_name", "atlas_gene_id",
                 "track_idx", "raw_score", "quantile_score", "metric", "track_name"])
    fi_rows = pd.concat(fis, ignore_index=True) if fis else pd.DataFrame(columns=fi_cols)
    del tops, fis
    log.info("%s score rows at your variant positions (all alleles)", f"{n_rows:,}")

    # --- small summary tables --------------------------------------------------
    summary = summarize_variants(vg, top, all_genes=a.all_genes)

    # all possible alleles at every queried position (AVI + track scores),
    # then flag which ones are input variants
    pa = build_position_alleles(fi_rows, top, tracks) if len(top) else pd.DataFrame()
    summary = attach_avi(summary, pa)
    if len(pa):
        id_col = "ID" if "ID" in vg.columns else "variant_key"
        inp = vg.assign(atlas_alt=vg["ALT"]).groupby(["pos_key", "atlas_alt"]).agg(
            input_variant_ids=(id_col, lambda x: ",".join(sorted({str(v) for v in x}))),
            input_ref=("REF", lambda x: ",".join(sorted(set(x)))),
            target_genes=("target_gene", lambda x: ",".join(sorted(set(x))))).reset_index()
        pa = pa.merge(inp, on=["pos_key", "atlas_alt"], how="left")
        pa["is_input_variant"] = pa["input_ref"].notna()
        pg = vg.groupby("pos_key")["target_gene"].agg(lambda x: ",".join(sorted(set(x)))).rename("position_target_genes")
        pa = pa.merge(pg, on="pos_key", how="left")
    pa.to_csv(f"{a.out_prefix}.position_alleles.tsv", sep="\t", index=False)
    log.info("Wrote %s.position_alleles.tsv (%d position x allele rows)", a.out_prefix, len(pa))

    outside = variants[~variants["variant_key"].isin(vg["variant_key"])].drop(columns=["_chrom_norm"])
    if len(outside):
        log.warning("%d variant rows lie outside all gene intervals (kept, atlas_available=False)", len(outside))
        summary = pd.concat([summary, outside.assign(atlas_available=False)], ignore_index=True)

    # --- per-variant status (vectorised) -------------------------------------
    queried_pos = set(in_win["pos_key"])
    returned_pos = set(top["pos_key"])
    in_window = set(vg["variant_key"])
    vk, pk = summary["variant_key"], summary["pos_key"]
    summary["atlas_status"] = np.select(
        [summary["atlas_available"].fillna(False).astype(bool),
         ~vk.isin(in_window),
         ~pk.isin(returned_pos) & pk.isin(queried_pos),
         ~pk.isin(returned_pos),
         ~summary["is_snv"].fillna(False).astype(bool)],
        ["scored", "outside_gene_intervals", "queried_no_atlas_result", "not_queried",
         "not_snv_position_alleles_only"],
        default="alt_not_in_atlas_alleles")

    q = pd.DataFrame({"pos_key": sorted(queried_pos)})
    if len(top):
        alle = top.groupby("pos_key").agg(atlas_ref=("atlas_ref", "first"),
                                          atlas_alleles=("atlas_alt", lambda x: ",".join(sorted(set(x)))))
        q = q.merge(alle, on="pos_key", how="left")
    q["returned"] = q["pos_key"].isin(returned_pos)
    q.to_csv(f"{a.out_prefix}.atlas_queried_positions.tsv", sep="\t", index=False)

    counts = summary.drop_duplicates("variant_key")["atlas_status"].value_counts()
    log.info("Variant status: %s", ", ".join(f"{k}={v}" for k, v in counts.items()))
    summary.to_csv(f"{a.out_prefix}.variants_atlas_summary.tsv", sep="\t", index=False)
    log.info("Wrote %s.variants_atlas_summary.tsv", a.out_prefix)

    if long_writer is None:
        log.info("Skipping long table (--no-long)")
    else:
        log.info("Wrote %s (%s rows; join to the summary on variant_key for annotations)",
                 long_writer.path, f"{long_writer.n:,}")
    log.info("%d / %d variants have Atlas scores for their observed ALT",
             int(summary.loc[summary["atlas_status"] == "scored", "variant_key"].nunique()),
             variants["variant_key"].nunique())
    return summary

if __name__ == "__main__":
    main()