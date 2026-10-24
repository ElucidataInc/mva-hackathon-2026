#!/usr/bin/env python3
"""
Step 4 - network proximity + direction-of-effect + OT/trial evidence + tiers.

Needs the STRING human flat files (download once from string-db.org > Download > Homo sapiens 9606):
  9606.protein.links.v12.0.txt.gz   and   9606.protein.info.v12.0.txt.gz

Method (Guney et al. 2016 style): for each drug, d_c = mean over its targets of the shortest-path
distance to the nearest disease-module gene; z-score vs degree-matched random target sets.
Negative z = drug targets sit closer to the module than chance.

v3 additions: conservative mechanism direction, current baseline expression, and separate clinical toxicity
flags from broader safety-review signals. These tiers remain hypotheses, not clinical recommendations.

Usage:
  python 04_rank.py --links 9606.protein.links.v12.0.txt.gz --info 9606.protein.info.v12.0.txt.gz
"""
import argparse
import os

import networkx as nx
import numpy as np
import pandas as pd

MITOTIC = {"TUBB", "TUBB1", "TUBB3", "TUBA1A", "AURKA", "AURKB", "PLK1", "TTK", "KIF11", "CDK1",
           "BUB1", "BUB1B", "MAD2L1", "CDC20", "CENPE", "NEK2"}


def load_string(links, info, cutoff):
    inf = pd.read_csv(info, sep="\t")
    inf.columns = [c.lstrip("#") for c in inf.columns]
    name = dict(zip(inf["string_protein_id"], inf["preferred_name"]))
    L = pd.read_csv(links, sep=" ")
    L = L[L.combined_score >= cutoff]
    L = L.assign(a=L.protein1.map(name), b=L.protein2.map(name)).dropna(subset=["a", "b"])
    G = nx.from_pandas_edgelist(L, "a", "b")
    G.remove_edges_from(nx.selfloop_edges(G))
    lcc = max(nx.connected_components(G), key=len)
    return G.subgraph(lcc).copy()


def proximity(G, module, drug_targets, n_perm=1000, seed=7):
    """Return {drug: (d_c, z, p, n_targets_in_graph)}; degree-matched null."""
    rng = np.random.default_rng(seed)
    mod = [m for m in module if m in G]
    dist = nx.multi_source_dijkstra_path_length(G, mod)          # nearest-module distance for every node
    nodes = np.array(list(G.nodes))
    d_all = np.array([dist[n] for n in nodes], dtype=float)
    deg = np.array([G.degree(n) for n in nodes])
    q = min(100, max(2, len(nodes) // 50))
    bins = pd.qcut(pd.Series(deg).rank(method="first"), q=q, labels=False).to_numpy()
    members = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    idx = {n: i for i, n in enumerate(nodes)}
    out = {}
    for drug, tg in drug_targets.items():
        t = [idx[g] for g in tg if g in idx]
        if not t:
            out[drug] = (np.nan, np.nan, np.nan, 0)
            continue
        obs = d_all[t].mean()
        null = np.empty(n_perm)
        pools = [members[bins[i]] for i in t]
        for k in range(n_perm):
            null[k] = np.mean([d_all[rng.choice(p)] for p in pools])
        sd = null.std()
        z = (obs - null.mean()) / sd if sd > 0 else 0.0
        p = (np.sum(null <= obs) + 1) / (n_perm + 1)
        out[drug] = (obs, z, p, len(t))
    return out


def bh(p):
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty_like(pv)
    tmp[order] = np.minimum(ranked, 1)
    q[ok] = tmp
    return q


def direction_for_row(r, lof_genes, rescue, auto):
    """Mechanism-aware hypothesis flags; rescue/activation is not proof of functional rescue."""
    g, a = r["gene"], str(r.get("action_class", "UNK")).upper()
    flags = []
    relevant_gene = g in lof_genes or g in rescue or g in auto or g in MITOTIC
    if a == "CONFLICT":
        flags.append("conflicting_mechanism")
    elif a not in ("INHIB", "ACT"):
        flags.append("mechanism_unknown_relevant" if relevant_gene else "mechanism_unknown_network")

    # A drug's nominal action on a LoF gene does not establish rescue of the patient's allele.
    if g in lof_genes:
        if a == "INHIB":
            flags.append("lof_inhibition_risk")
        elif a == "ACT":
            flags.append("lof_activation_unvalidated")

    if g in rescue:
        want, ctx = rescue[g]
        if a == want and a in ("INHIB", "ACT") and a != "CONFLICT":
            flags.append(f"manual_rescue_hypothesis:{ctx}")
        elif a in ("INHIB", "ACT") and a != want:
            flags.append("opposes_manual_rescue")
    elif g in auto:
        want, ctx = auto[g]
        if a == want and a in ("INHIB", "ACT"):
            flags.append(f"auto_rescue_hypothesis:{ctx}")

    if g in MITOTIC and a == "INHIB":
        flags.append("mitotic_disruptor_review")
    return ";".join(dict.fromkeys(flags)) if flags else "network_context_only"


def summarise_direction(flag_series):
    f = ";".join(str(x) for x in flag_series)
    bad = "opposes_manual_rescue" in f or "lof_inhibition_risk" in f
    manual = "manual_rescue_hypothesis:" in f
    auto = "auto_rescue_hypothesis:" in f
    uncertain = any(x in f for x in ("mechanism_unknown_relevant", "conflicting_mechanism", "lof_activation_unvalidated"))
    mitotic = "mitotic_disruptor_review" in f
    if bad and (manual or auto):
        return "mixed_review"
    if bad or uncertain:
        return "review"
    if mitotic:
        return "caution_mitotic"
    if manual:
        return "manual_hypothesis"
    if auto:
        return "auto_hypothesis"
    return "network_context"


def tier(r, z_cut):
    ph = r["max_phase"] if pd.notna(r["max_phase"]) else 0
    close = pd.notna(r["z"]) and r["z"] <= z_cut
    d = r["direction"]
    if d in ("avoid", "mixed", "review", "mixed_review"):
        return "X-avoid/review"
    if d == "caution_mitotic":
        return "3-caution"
    if ph >= 3 and (d == "manual_hypothesis" or close):
        t = "1"
    elif d == "auto_hypothesis":
        t = "2" if ph >= 1 else "3"          # inferred network rescue never reaches tier 1 alone
    elif (ph >= 1 and (d == "manual_hypothesis" or close)) or d == "manual_hypothesis":
        t = "2"
    else:
        t = "3"
    # Broad safety evidence triggers review/demotion, not a declaration that a drug is unsafe.
    if t == "1" and bool(r.get("safety_review_flag", False)):
        return "2"
    return t


def tf(x):
    return x is True or str(x).lower() == "true"


def evidence_for_drug(tg, te, de_row):
    """Return notes, hard clinical-toxicity flag, broad safety-review flag, related-trial count."""
    notes, hard_tox, safety_review = [], False, False
    # Only call expression negative when the new API actually matched the requested tissue(s).
    if "expression_status" in te.columns and te["expression_status"].notna().any():
        detected = [g for g in tg if g in te.index and te.at[g, "expression_status"] == "detected_in_matched_tissue"]
        not_detected = [g for g in tg if g in te.index and te.at[g, "expression_status"] == "not_detected_in_matched_tissue"]
        unknown_expr = [g for g in tg if g in te.index and str(te.at[g, "expression_status"]).startswith("no_")]
        if detected:
            notes.append("baseline expression detected in matched tissue(s): " + ",".join(detected))
        if not_detected:
            notes.append("no baseline expression detected in matched tissue rows (not proof of absence): " + ",".join(not_detected))
        if unknown_expr:
            notes.append("expression unknown/no matching tissue rows: " + ",".join(unknown_expr))
    else:
        notes.append("expression evidence unavailable: rerun Step 3 with baselineExpression API")

    lethal = [g for g in tg if g in te.index and "ko_lethal" in te.columns and tf(te.at[g, "ko_lethal"])]
    safety = [g for g in tg if g in te.index and "safety_events_relevant" in te.columns and
              str(te.at[g, "safety_events_relevant"]).strip().lower() not in ("", "nan", "none")]
    if lethal:
        notes.append("mouse knockout lethality signal (biological caution; not equivalent to drug toxicity): " + ",".join(lethal))
        safety_review = True
    if safety:
        notes.append("target safety-liability signal: " + ",".join(safety))
        safety_review = True
    n_rel = 0
    if de_row is not None:
        def num(k):
            try:
                v = de_row.get(k)
                return float(v) if pd.notna(v) else 0.0
            except (TypeError, ValueError):
                return 0.0
        tox = int(num("n_tox_stop"))
        futility = int(num("n_futility_stop"))
        if tox:
            notes.append(f"{tox} trial(s) reported stopped for toxicity/safety (verify trial identity and reason)")
            hard_tox = True
            safety_review = True
        if futility:
            notes.append(f"{futility} trial(s) reported stopped for futility (not a safety finding)")
        if tf(de_row.get("faers_serious_flag")):
            notes.append("FAERS serious-event signal (spontaneous reports; not proof of causality): " + str(de_row.get("faers_serious")))
            safety_review = True
        n_rel = int(num("n_related_trials"))
        if n_rel:
            notes.append(f"{n_rel} trial(s) in related conditions (condition-name match; inspect individual trials)")
    return "; ".join(notes), hard_tox, safety_review, n_rel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--links", required=True)
    ap.add_argument("--info", required=True)
    ap.add_argument("--drugs", default="out/drug_target_long.csv")
    ap.add_argument("--seeds", default="seeds.csv")
    ap.add_argument("--rescue", default="rescue_targets.csv")
    ap.add_argument("--auto-rescue", default="out/auto_rescue_targets.csv")
    ap.add_argument("--target-evidence", default="out/target_evidence.csv")
    ap.add_argument("--drug-evidence", default="out/drug_evidence.csv")
    ap.add_argument("--seed-report", default="out/seed_report.csv")
    ap.add_argument("--module-tiers", default="AB")
    ap.add_argument("--lof-tiers", default="ABC",
                    help="tiers whose LoF genes make their inhibitors 'avoid' (tier D = passenger/context, excluded)")
    ap.add_argument("--string-cutoff", type=int, default=700)
    ap.add_argument("--perm", type=int, default=1000)
    ap.add_argument("--z-cut", type=float, default=-1.64)
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    seeds = pd.read_csv(a.seeds)
    resc = pd.read_csv(a.rescue)
    rescue = {r.gene: (r.desired_action, r.context) for r in resc.itertuples()}
    lof = set(seeds[(seeds.variant_mechanism == "LoF") & seeds.tier.isin(list(a.lof_tiers))].gene)
    variant_genes = set(seeds[seeds.in_variant26 == "yes"].gene) if "in_variant26" in seeds else set()
    auto = {}
    if os.path.exists(a.auto_rescue):
        ar = pd.read_csv(a.auto_rescue)
        auto = {r.gene: (r.desired_action, r.context) for r in ar.itertuples() if r.gene not in rescue}
        print(f"auto-rescue targets: {len(auto)} (manual: {len(rescue)})")
    te = pd.read_csv(a.target_evidence).set_index("gene") if os.path.exists(a.target_evidence) else None
    if te is not None:
        for c in ("ot_status", "expression_status", "expressed_affected", "ko_lethal", "n_safety_events", "safety_events_relevant"):
            if c not in te.columns:
                te[c] = None
    de = pd.read_csv(a.drug_evidence).set_index("drug_key") if os.path.exists(a.drug_evidence) else None
    if te is None or de is None:
        print("note: OT/trial evidence files missing - run 03_ot_evidence.py for safety/expression flags")
    module = seeds[seeds.tier.isin(list(a.module_tiers))].gene.tolist()

    df = pd.read_csv(a.drugs)
    df = df[df.drug_key.notna() & (df.drug_key.astype(str) != "")]
    df["direction_row"] = df.apply(lambda r: direction_for_row(r, lof, rescue, auto), axis=1)

    print("Loading STRING ...")
    G = load_string(a.links, a.info, a.string_cutoff)
    print(f"  interactome: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges; module genes in graph: "
          f"{sum(m in G for m in module)}/{len(module)}")

    targets = df.groupby("drug_key").gene.agg(lambda s: sorted(set(s))).to_dict()
    prox = proximity(G, module, targets, a.perm)

    def agg(g):
        names = g.drug_name[g.drug_name != ""]
        return pd.Series({
            "drug": names.iloc[0] if len(names) else g.name,
            "targets": ", ".join(sorted(set(g.gene))),
            "n_targets": g.gene.nunique(),
            "direct_seed_target": bool(set(g.gene) & set(module)),
            "hits_variant_gene": ", ".join(sorted(set(g.gene) & variant_genes)),
            "gene_roles": ", ".join(sorted(set(g.gene_role.dropna()))) if "gene_role" in g else "",
            "action_types": ", ".join(sorted(set(x for x in g.action_type.dropna().astype(str) if x))),
            "max_phase": g.max_phase.max(),
            "approved": bool((g.approved == True).any()),  # noqa: E712
            "best_pchembl": g.pchembl.max(),
            "sources": ", ".join(sorted(set(g.source))),
            "evidence": ", ".join(sorted(set(g.evidence))),
            "direction_detail": " | ".join(sorted(set(f"{x.gene}:{x.direction_row}" for x in g.itertuples()))),
            "direction": summarise_direction(g.direction_row),
            "contexts": ", ".join(sorted(set(c for f in g.direction_row for c in
                                             [s.split(":")[1] for s in f.split(";")
                                              if s.startswith(("rescue:", "auto:"))]))),
        })

    res = df.groupby("drug_key").apply(agg, include_groups=False).reset_index()
    res["d_c"] = res.drug_key.map(lambda k: prox[k][0])
    res["z"] = res.drug_key.map(lambda k: prox[k][1])
    res["p"] = res.drug_key.map(lambda k: prox[k][2])
    res["fdr"] = bh(res.p.values)
    ev = []
    for r in res.itertuples():
        tg = [x.strip() for x in r.targets.split(",")]
        drow = de.loc[r.drug_key].to_dict() if (de is not None and r.drug_key in de.index) else None
        ev.append(evidence_for_drug(tg, te, drow) if te is not None else ("expression/safety evidence unavailable", False, False, 0))
    res["evidence_notes"] = [e[0] for e in ev]
    res["hard_flag"] = [e[1] for e in ev]  # specifically a clinical trial toxicity-stop signal
    res["safety_review_flag"] = [e[2] for e in ev]
    res["n_related_trials"] = [e[3] for e in ev]
    res["tier"] = res.apply(lambda r: tier(r, a.z_cut), axis=1)
    # setting: tumour-only rescue drugs are NOT germline options - keep them visibly separate
    def setting(c):
        c = str(c)
        if c == "tumour":
            return "tumour-only"
        if "germline" in c:
            return "germline-pathway"
        if "pathway_inferred" in c:
            return "pathway-inferred"
        return "network/other"
    res["setting"] = res.contexts.map(setting)
    res = res.sort_values(["tier", "z"], ascending=[True, True], na_position="last")
    cols = ["tier", "setting", "drug", "drug_key", "targets", "direction", "contexts", "evidence_notes",
            "hits_variant_gene", "n_related_trials", "z", "p", "fdr", "d_c", "max_phase",
            "approved", "action_types", "best_pchembl", "direct_seed_target", "gene_roles", "sources", "evidence",
            "direction_detail", "hard_flag", "safety_review_flag"]
    res = res[cols]

    with pd.ExcelWriter(f"{a.out}/drug_repurposing_ranked.xlsx") as xw:
        res[~res.tier.str.startswith("X")].to_excel(xw, sheet_name="ranked", index=False)
        res[res.tier.str.startswith("X")].to_excel(xw, sheet_name="avoid_or_review", index=False)
        df.to_excel(xw, sheet_name="gene_drug_evidence", index=False)
        if te is not None:
            te.reset_index().to_excel(xw, sheet_name="target_evidence", index=False)
        if de is not None:
            de.reset_index().to_excel(xw, sheet_name="drug_evidence", index=False)
        if os.path.exists(a.seed_report):
            pd.read_csv(a.seed_report).to_excel(xw, sheet_name="seed_report", index=False)
        pd.DataFrame({"parameter": ["module genes", "LoF-avoid tiers", "string cutoff", "permutations", "z cutoff", "note"],
                      "value": [", ".join(module), a.lof_tiers, a.string_cutoff, a.perm, a.z_cut,
                                "computational hypotheses only; not validated"]}).to_excel(xw, sheet_name="params", index=False)
    res.to_csv(f"{a.out}/drug_repurposing_ranked.csv", index=False)
    print(res.tier.value_counts().to_string())
    print(res.head(20)[["tier", "setting", "drug", "targets", "direction", "z", "max_phase"]].to_string(index=False))


if __name__ == "__main__":
    main()
