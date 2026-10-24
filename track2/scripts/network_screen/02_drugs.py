#!/usr/bin/env python3
"""
Step 2 - pull drug-target evidence for seeds + top network neighbours + rescue targets.

Sources: UniProt (gene->accession), ChEMBL (mechanism of action, approval phase, binders), DGIdb (interactions).
Gene list = seeds (all tiers, incl. the 26-variant genes) + top hop-1 and hop-2 STRING neighbours + manual rescue
targets + top auto-rescue targets from the OmniPath directed layer.
Output: out/drug_target_long.csv  (one row per gene-drug pair, with action class INHIB / ACT / UNK, plus gene_role)
        out/query_genes.csv       (gene, role) - reused by 03_ot_evidence.py

Usage:
  python 02_drugs.py --nodes out/nodes.csv --seeds seeds.csv --rescue rescue_targets.csv --top-neighbours 25
"""
import argparse
import os
import time

import pandas as pd
import requests

from common import is_symbol

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
DGIDB = "https://dgidb.org/api/graphql"

INHIB = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR",
         "DEGRADER", "DISRUPTING AGENT", "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "INVERSE AGONIST",
         "SEQUESTERING AGENT", "INHIBITOR (COMPETITIVE)"}
ACT = {"ACTIVATOR", "AGONIST", "PARTIAL AGONIST", "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR",
       "STABILISER", "STABILIZER", "OPENER", "INDUCER"}


def action_class(txt):
    """Classify mechanism text conservatively; contradictory annotations stay explicit."""
    parts = [p.strip().upper() for p in str(txt or "").replace("|", ";").replace(",", ";").split(";") if p.strip()]
    has_inhib = any(p in INHIB for p in parts)
    has_act = any(p in ACT for p in parts)
    if has_inhib and has_act:
        return "CONFLICT"
    if has_inhib:
        return "INHIB"
    if has_act:
        return "ACT"
    return "UNK"


def get_json(url, params=None, retries=3):
    last = None
    for i in range(retries):
        try:
            r = requests.get(url, params=params, timeout=120, headers={"Accept": "application/json"})
            if r.status_code == 404:
                return {}
            r.raise_for_status()
            return r.json()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    print(f"  ! GET failed {url} {params}: {last}")
    return {}


def uniprot_acc(gene):
    js = get_json(UNIPROT, {"query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
                            "fields": "accession", "format": "json", "size": 1})
    res = js.get("results", [])
    return res[0]["primaryAccession"] if res else None


def chembl_target(acc):
    js = get_json(f"{CHEMBL}/target.json", {"target_components__accession": acc,
                                            "target_type": "SINGLE PROTEIN", "limit": 1})
    t = js.get("targets", [])
    return t[0]["target_chembl_id"] if t else None


def chembl_mechanisms(tid):
    return get_json(f"{CHEMBL}/mechanism.json", {"target_chembl_id": tid, "limit": 1000}).get("mechanisms", [])


def chembl_binders(tid, min_pchembl, top=15):
    js = get_json(f"{CHEMBL}/activity.json", {"target_chembl_id": tid, "pchembl_value__gte": min_pchembl,
                                              "limit": 1000})
    best = {}
    for a in js.get("activities", []):
        m, p = a.get("molecule_chembl_id"), a.get("pchembl_value")
        if m and p:
            best[m] = max(best.get(m, 0), float(p))
    return sorted(best.items(), key=lambda x: -x[1])[:top]


def chembl_molecules(ids):
    out = {}
    ids = sorted(set(i for i in ids if i))
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        js = get_json(f"{CHEMBL}/molecule.json", {"molecule_chembl_id__in": ",".join(chunk), "limit": 1000})
        for m in js.get("molecules", []):
            try:
                ph = float(m.get("max_phase")) if m.get("max_phase") is not None else None
            except (TypeError, ValueError):
                ph = None
            out[m["molecule_chembl_id"]] = {"name": m.get("pref_name"), "max_phase": ph,
                                            "type": m.get("molecule_type"),
                                            "withdrawn": bool(m.get("withdrawn_flag"))}
    return out


def dgidb(genes):
    q = """query($names:[String!]){genes(names:$names){nodes{name interactions{
             drug{name conceptId approved} interactionScore interactionTypes{type directionality}}}}}"""
    rows = []
    for i in range(0, len(genes), 40):
        chunk = genes[i:i + 40]
        try:
            r = requests.post(DGIDB, json={"query": q, "variables": {"names": chunk}}, timeout=120)
            r.raise_for_status()
            nodes = r.json().get("data", {}).get("genes", {}).get("nodes", [])
        except Exception as e:  # noqa: BLE001
            print("  ! DGIdb failed (check schema at dgidb.org/api/graphql):", e)
            continue
        for n in nodes:
            for it in n.get("interactions", []):
                d = it.get("drug") or {}
                types = ";".join(t.get("type") or "" for t in it.get("interactionTypes") or [])
                cid = (d.get("conceptId") or "")
                rows.append({"gene": n["name"], "source": "DGIdb",
                             "drug_id": cid.split(":")[-1].upper() if cid.lower().startswith("chembl:") else "",
                             "drug_name": (d.get("name") or "").upper(), "action_type": types,
                             "moa": "", "max_phase": None, "approved": bool(d.get("approved")),
                             "pchembl": None, "evidence": "dgidb",
                             "dgidb_score": it.get("interactionScore")})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes", default="out/nodes.csv")
    ap.add_argument("--seeds", default="seeds.csv")
    ap.add_argument("--rescue", default="rescue_targets.csv")
    ap.add_argument("--auto-rescue", default="out/auto_rescue_targets.csv")
    ap.add_argument("--top-neighbours", type=int, default=25, help="top hop-1 neighbours")
    ap.add_argument("--top-hop2", type=int, default=25, help="top hop-2 neighbours")
    ap.add_argument("--max-auto", type=int, default=120, help="max auto-rescue genes to look up")
    ap.add_argument("--min-pchembl", type=float, default=6.5)
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    seeds = pd.read_csv(a.seeds)
    nodes = pd.read_csv(a.nodes)
    nb = nodes[nodes.role == "neighbour"]
    if "hop" in nb.columns:
        n1 = nb[nb.hop == 1].head(a.top_neighbours).gene.tolist()
        n2 = nb[nb.hop == 2].head(a.top_hop2).gene.tolist()
    else:
        n1, n2 = nb.head(a.top_neighbours).gene.tolist(), []
    rescue = pd.read_csv(a.rescue).gene.tolist()
    auto = []
    if os.path.exists(a.auto_rescue):
        ar = pd.read_csv(a.auto_rescue)
        auto = ar.sort_values(["level", "support"], ascending=[True, False]).head(a.max_auto).gene.tolist()

    role = {}
    for g, r in ((g, "seed") for g in seeds.gene):
        role.setdefault(g, r)
    for lst, r in ((rescue, "rescue"), (n1, "hop1"), (n2, "hop2"), (auto, "auto_rescue")):
        for g in lst:
            role.setdefault(g, r)
    genes = list(role)
    pd.DataFrame({"gene": genes, "role": [role[g] for g in genes]}).to_csv(f"{a.out}/query_genes.csv", index=False)
    lookup = [g for g in genes if is_symbol(g)]
    binder_genes = set(seeds[seeds.tier.isin(["A", "B"])].gene) | set(n1[:10])
    print(f"{len(lookup)} genes to query ({len(seeds)} seeds, {len(n1)} hop1, {len(n2)} hop2, "
          f"rescue {len(set(rescue))}, auto-rescue {len(auto)}); {len(genes) - len(lookup)} without symbol skipped")

    rows, mol_ids, status = [], [], []
    for g in lookup:
        acc = uniprot_acc(g)
        tid = chembl_target(acc) if acc else None
        status.append({"gene": g, "uniprot": acc, "chembl_target": tid})
        if not tid:
            print(f"  {g}: no ChEMBL single-protein target")
            continue
        mechs = chembl_mechanisms(tid)
        for m in mechs:
            rows.append({"gene": g, "source": "ChEMBL", "drug_id": m.get("molecule_chembl_id"), "drug_name": "",
                         "action_type": m.get("action_type"), "moa": m.get("mechanism_of_action"),
                         "max_phase": None, "approved": None, "pchembl": None, "evidence": "mechanism"})
            mol_ids.append(m.get("molecule_chembl_id"))
        if g in binder_genes:
            for mid, p in chembl_binders(tid, a.min_pchembl):
                rows.append({"gene": g, "source": "ChEMBL", "drug_id": mid, "drug_name": "", "action_type": "",
                             "moa": "bioactivity binder", "max_phase": None, "approved": None,
                             "pchembl": p, "evidence": "binding"})
                mol_ids.append(mid)
        print(f"  {g}: {acc} {tid} mechanisms={len(mechs)}")

    mols = chembl_molecules(mol_ids)
    for r in rows:
        mi = mols.get(r["drug_id"], {})
        r["drug_name"] = (mi.get("name") or r["drug_id"] or "").upper()
        r["max_phase"] = mi.get("max_phase")
        r["approved"] = (mi.get("max_phase") == 4.0) if mi else None
        r["withdrawn"] = mi.get("withdrawn", False)
    rows += dgidb(lookup)

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("No drug rows retrieved - check network access / API errors above")
    df["action_class"] = df.action_type.map(action_class)
    df["drug_key"] = df.apply(lambda r: r.drug_id if r.drug_id else r.drug_name, axis=1)
    df["gene_role"] = df.gene.map(role)
    df.to_csv(f"{a.out}/drug_target_long.csv", index=False)
    pd.DataFrame(status).to_csv(f"{a.out}/gene_target_map.csv", index=False)
    print(f"\n{len(df)} gene-drug rows, {df.drug_key.nunique()} unique drugs -> {a.out}/drug_target_long.csv")


if __name__ == "__main__":
    main()
