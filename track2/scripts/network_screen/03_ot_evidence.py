#!/usr/bin/env python3
"""
Step 3 - clinical / safety / biology evidence for the drugs and targets found in step 2.

Sources
  Open Targets Platform (GraphQL): target expression in affected tissues, target safety liabilities, mouse
      knockout phenotypes, associated diseases matching the patient context, FAERS adverse events per drug.
  ClinicalTrials.gov v2: trials per drug name -> terminated / stopped-for-toxicity / stopped-for-futility /
      trials in related conditions.
Patient context (tissues, disease regexes, safety regexes) comes from context.json (see common.py).

Outputs (in --out):
  target_evidence.csv   one row per target gene
  drug_evidence.csv     one row per drug_key (candidate drugs only: phase >= 1 or approved, capped by --max-drugs)

Every network call fails soft: a failed lookup leaves the evidence columns empty and is printed, never fatal.

Usage:
  python 03_ot_evidence.py --drugs out/drug_target_long.csv --out out
"""
import argparse
import os
import re

import pandas as pd

from common import http, is_symbol, load_context

OT = "https://api.platform.opentargets.org/api/v4/graphql"
CT = "https://clinicaltrials.gov/api/v2/studies"

TARGET_Q = """query($id:String!){target(ensemblId:$id){id approvedSymbol
  safetyLiabilities{event}
  mousePhenotypes{modelPhenotypeLabel}
  baselineExpression(page:{index:0,size:3000}){rows{median unit datasourceId datatypeId
    tissueBiosample{biosampleId biosampleName}
    celltypeBiosample{biosampleId biosampleName}}}
  associatedDiseases(page:{index:0,size:50}){rows{score disease{name}}}}}"""
FAERS_Q = "query($id:String!){drug(chemblId:$id){adverseEvents{rows{name count logLR}}}}"
MAP_Q = 'query($t:[String!]!){mapIds(queryTerms:$t,entityNames:["target"]){mappings{term hits{id name entity}}}}'
SEARCH_Q = ('query($q:String!){search(queryString:$q,entityNames:["target"],page:{index:0,size:5})'
            '{hits{id name entity}}}')


def gql(query, variables=None):
    js = http("POST", OT, json={"query": query, "variables": variables or {}}).json()
    if js.get("errors"):
        # Fail closed: partial GraphQL data must not silently turn missing expression into a negative result.
        raise RuntimeError(str(js["errors"])[:500])
    return js.get("data") or {}


def map_ids(genes):
    out = {}
    for i in range(0, len(genes), 50):
        chunk = genes[i:i + 50]
        try:
            d = gql(MAP_Q, {"t": chunk})
            for m in d["mapIds"]["mappings"]:
                hits = [h for h in m.get("hits") or [] if h.get("entity") == "target"]
                if hits:
                    out[m["term"]] = hits[0]["id"]
        except Exception as e:  # noqa: BLE001
            print("  ! mapIds failed, falling back to search:", e)
    for g in genes:
        if g in out:
            continue
        try:
            for h in gql(SEARCH_Q, {"q": g}).get("search", {}).get("hits", []):
                if h.get("entity") == "target" and str(h.get("name", "")).upper() == g.upper():
                    out[g] = h["id"]
                    break
        except Exception as e:  # noqa: BLE001
            print(f"  ! OT search failed for {g}: {e}")
    return out


def parse_target(gene, data, ctx):
    """Flatten Open Targets data while treating absent/unmatched expression as unknown."""
    t = (data or {}).get("target")
    if not t:
        return {"gene": gene, "ot_status": "no target data", "expression_status": "unknown"}

    affected = [str(x).strip().lower() for x in ctx.get("affected_tissues", []) if str(x).strip()]
    block = t.get("baselineExpression") or {}
    expr_rows = block.get("rows") or []
    matched, detected = [], []
    for row in expr_rows:
        tissue = ((row.get("tissueBiosample") or {}).get("biosampleName") or "").strip()
        cell = ((row.get("celltypeBiosample") or {}).get("biosampleName") or "").strip()
        labels = [x for x in (tissue, cell) if x]
        combined = " | ".join(labels).lower()
        if affected and any(term in combined for term in affected):
            item = {"tissue": tissue, "cell_type": cell, "median": row.get("median"),
                    "unit": row.get("unit") or "", "datasource": row.get("datasourceId") or "",
                    "datatype": row.get("datatypeId") or ""}
            matched.append(item)
            try:
                if row.get("median") is not None and float(row["median"]) > 0:
                    detected.append(item)
            except (TypeError, ValueError):
                pass

    # Median values have source-specific units (TPM, pseudobulk CPM, proteomics units).
    # We use them only as a presence signal, never compare their magnitudes across sources.
    if matched and detected:
        expr_status, expr_bool = "detected_in_matched_tissue", True
    elif matched:
        expr_status, expr_bool = "not_detected_in_matched_tissue", False
    elif len(expr_rows) >= 3000:
        expr_status, expr_bool = "no_match_in_first_3000_rows", None
    else:
        expr_status, expr_bool = "no_matching_tissue_rows", None

    events = [x.get("event") or "" for x in t.get("safetyLiabilities") or []]
    sev = sorted({e for e in events if re.search(ctx["target_safety_regex"], e, re.I)})
    labels = sorted({x.get("modelPhenotypeLabel") or "" for x in t.get("mousePhenotypes") or []})
    lethal = [l for l in labels if re.search(ctx["ko_lethal_regex"], l, re.I)]
    relev = [l for l in labels if re.search(ctx["ko_relevant_regex"], l, re.I)]
    rel = []
    for r in (t.get("associatedDiseases") or {}).get("rows") or []:
        nm = (r.get("disease") or {}).get("name") or ""
        if re.search(ctx["related_disease_regex"], nm, re.I):
            rel.append((r.get("score") or 0, nm))
    rel.sort(reverse=True)
    return {"gene": gene, "ot_status": "ok", "expression_status": expr_status,
            "expressed_affected": expr_bool, "n_affected_tissues_expressed": len({x['tissue'].lower() for x in detected if x['tissue']}),
            "affected_tissues_expressed": "; ".join(sorted({x['tissue'] for x in detected if x['tissue']})[:10]),
            "expression_match_rows": len(matched),
            "expression_sources": "; ".join(sorted({x['datasource'] for x in matched if x['datasource']})),
            "expression_evidence": "; ".join(f"{x['tissue']} / {x['cell_type']} (median={x['median']} {x['unit']}; {x['datasource']})" for x in detected[:10]),
            "n_safety_events": len(events), "safety_events_relevant": "; ".join(sev[:6]),
            "ko_lethal": bool(lethal), "ko_relevant": "; ".join(relev[:6]),
            "related_disease_hits": "; ".join(f"{n} ({sc:.2f})" for sc, n in rel[:5])}


def summarise_trials(studies, ctx):
    """Pure function: ClinicalTrials.gov v2 'studies' list -> counts."""
    c = dict(n_trials=len(studies), n_completed=0, n_terminated=0, n_tox_stop=0, n_futility_stop=0,
             n_related_trials=0)
    for s in studies:
        ps = s.get("protocolSection") or {}
        st = ps.get("statusModule") or {}
        status, why = st.get("overallStatus") or "", st.get("whyStopped") or ""
        conds = " ".join((ps.get("conditionsModule") or {}).get("conditions") or [])
        if status.upper() == "COMPLETED":
            c["n_completed"] += 1
        if re.search(ctx["terminated_status_regex"], status, re.I):
            c["n_terminated"] += 1
            c["n_tox_stop"] += bool(re.search(ctx["toxicity_regex"], why, re.I))
            c["n_futility_stop"] += bool(re.search(ctx["futility_regex"], why, re.I))
        c["n_related_trials"] += bool(re.search(ctx["related_disease_regex"], conds, re.I))
    return c


def trials_for(name, ctx):
    r = http("GET", CT, params={"query.intr": name, "pageSize": 100, "format": "json",
                                "fields": "NCTId|OverallStatus|WhyStopped|Condition|Phase"})
    return summarise_trials(r.json().get("studies", []), ctx)


def faers_for(chembl_id, ctx):
    rows = (gql(FAERS_Q, {"id": chembl_id}).get("drug") or {}).get("adverseEvents", {}) or {}
    rows = sorted(rows.get("rows") or [], key=lambda x: -(x.get("logLR") or 0))[:ctx["faers_top_n"]]
    names = [r.get("name") or "" for r in rows]
    serious = [n for n in names if re.search(ctx["serious_ae_regex"], n, re.I)]
    return {"faers_top": "; ".join(names), "faers_serious": "; ".join(serious), "faers_serious_flag": bool(serious)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drugs", default="out/drug_target_long.csv")
    ap.add_argument("--context", default="context.json")
    ap.add_argument("--max-drugs", type=int, default=150, help="max candidate drugs to look up trials/FAERS for")
    ap.add_argument("--no-trials", action="store_true")
    ap.add_argument("--no-faers", action="store_true")
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ctx = load_context(a.context)

    df = pd.read_csv(a.drugs)
    genes = sorted(g for g in df.gene.dropna().unique() if is_symbol(g))
    print(f"target evidence for {len(genes)} genes ...")
    ids = map_ids(genes)
    rows = []
    for g in genes:
        if g not in ids:
            rows.append({"gene": g, "ot_status": "not found in Open Targets"})
            continue
        try:
            row = parse_target(g, gql(TARGET_Q, {"id": ids[g]}), ctx)
            row["ensembl"] = ids[g]
        except Exception as e:  # noqa: BLE001
            print(f"  ! target {g}: {e}")
            row = {"gene": g, "ensembl": ids[g], "ot_status": "failed"}
        rows.append(row)
    pd.DataFrame(rows).to_csv(f"{a.out}/target_evidence.csv", index=False)
    print(f"  {sum(r.get('ot_status') == 'ok' for r in rows)}/{len(genes)} ok -> {a.out}/target_evidence.csv")

    d = df[df.drug_key.notna() & (df.drug_key.astype(str) != "")].copy()
    d["max_phase"] = pd.to_numeric(d.max_phase, errors="coerce")
    cand = (d.groupby("drug_key").agg(drug=("drug_name", lambda s: next((x for x in s if isinstance(x, str) and x), "")),
                                      max_phase=("max_phase", "max"))
            .reset_index())
    cand = cand[(cand.max_phase.fillna(0) >= 1)].sort_values("max_phase", ascending=False).head(a.max_drugs)
    print(f"drug evidence for {len(cand)} candidate drugs (phase >= 1) ...")
    out = []
    for r in cand.itertuples():
        row = {"drug_key": r.drug_key, "drug": r.drug}
        key = str(r.drug_key)
        if not a.no_faers and key.upper().startswith("CHEMBL"):
            try:
                row.update(faers_for(key.upper(), ctx))
            except Exception as e:  # noqa: BLE001
                print(f"  ! FAERS {key}: {e}")
        if not a.no_trials and r.drug and not str(r.drug).upper().startswith("CHEMBL"):
            try:
                row.update(trials_for(str(r.drug).title(), ctx))
            except Exception as e:  # noqa: BLE001
                print(f"  ! trials {r.drug}: {e}")
        out.append(row)
    pd.DataFrame(out).to_csv(f"{a.out}/drug_evidence.csv", index=False)
    print(f"  -> {a.out}/drug_evidence.csv")


if __name__ == "__main__":
    main()
