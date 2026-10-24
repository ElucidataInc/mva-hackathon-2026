#!/usr/bin/env python3
"""
Step 1 - networks around the seed genes.

  A. STRING, 2 hops:  hop 1 = seeds + STRING-added interactors (score >= --score),
                      hop 2 = partners of the best hop-1 nodes, kept only at score >= --score2 (default 0.9).
  B. OmniPath, directed, --levels up and down (default 2): regulators of each seed (upstream) and genes the
     seed regulates (downstream), with activation/inhibition sign. Level 2 needs >= --min-support2 sources and
     skips hub nodes. Signs are propagated to produce AUTO rescue targets (see auto_rescue_targets.csv).

Outputs (in --out):
  edges.csv nodes.csv network.graphml network.sif network.png reactome_enrichment.csv   (STRING layer)
  seed_report.csv           per seed: in network?, links, distance to the A/B module, directed up/down counts
  directed_edges.csv directed_nodes.csv directed.graphml directed.sif   (OmniPath layer)
  auto_rescue_targets.csv   sign-propagated desired action for downstream (LoF seeds) / upstream (unknown seeds)

Usage:
  python 01_network.py --seeds seeds.csv --out out --score 0.7 --add-nodes 40 --score2 0.9
"""
import argparse
import io
import os
import sys

import networkx as nx
import pandas as pd
from networkx.algorithms import community

from common import http, is_symbol

STRING = "https://string-db.org/api"
REACTOME = "https://reactome.org/AnalysisService/identifiers/projection"
OMNIPATH = "https://omnipathdb.org/interactions"
CALLER = "mva_repurposing_hackathon"


# ----------------------------------------------------------------------------- STRING
def post_json(url, data=None, **kw):
    return http("POST", url, data=data, **kw).json()


def string_network(genes, score, add_nodes):
    data = {"identifiers": "\r".join(genes), "species": 9606, "required_score": int(round(score * 1000)),
            "add_nodes": add_nodes, "network_type": "functional", "caller_identity": CALLER}
    return post_json(f"{STRING}/json/network", data=data)


def string_partners(genes, score, limit):
    data = {"identifiers": "\r".join(genes), "species": 9606, "required_score": int(round(score * 1000)),
            "limit": limit, "caller_identity": CALLER}
    return post_json(f"{STRING}/json/interaction_partners", data=data)


def add_rows(G, rows):
    for e in rows:
        a, b, s = e["preferredName_A"], e["preferredName_B"], float(e["score"])
        if a == b:
            continue
        if G.has_edge(a, b):
            G[a][b]["score"] = max(G[a][b]["score"], s)
        else:
            G.add_edge(a, b, score=s)
    return G


def finish_graph(G):
    for a, b, d in G.edges(data=True):
        d["dist"] = 1.0 / max(d["score"], 1e-6)
    return G


def filter_hop2(G, hop, score2):
    """Edges that touch a hop-2 node must reach score2; hop<=1 pairs keep their (lower) cutoff."""
    drop = [(u, v) for u, v, d in G.edges(data=True)
            if (hop.get(u, 2) == 2 or hop.get(v, 2) == 2) and d["score"] < score2]
    G.remove_edges_from(drop)
    return G


def rank_hop1(G1, seeds):
    sc = {}
    for n in G1.nodes:
        if n in seeds:
            continue
        s_nb = [m for m in G1.neighbors(n) if m in seeds]
        sc[n] = (len(s_nb), sum(G1[n][m]["score"] for m in s_nb))
    return sorted(sc, key=lambda n: sc[n], reverse=True)


def node_table(G, seeds_df, hop):
    seed_tier = dict(zip(seeds_df.gene, seeds_df.tier))
    seeds = set(seed_tier)
    for g in seeds:
        if g not in G and is_symbol(g):
            G.add_node(g)
    btw = nx.betweenness_centrality(G, weight="dist", normalized=True)
    comms = community.greedy_modularity_communities(G, weight="score")
    comm_of = {n: i for i, c in enumerate(comms) for n in c}
    recs = []
    for n in G.nodes:
        nb = list(G.neighbors(n))
        s_nb = [m for m in nb if m in seeds]
        h1_nb = [m for m in nb if hop.get(m, 2) <= 1]
        recs.append({"gene": n, "role": "seed" if n in seeds else "neighbour",
                     "hop": 0 if n in seeds else hop.get(n, 2), "tier": seed_tier.get(n, ""),
                     "degree": len(nb), "seed_links": len(s_nb),
                     "seed_score_sum": round(sum(G[n][m]["score"] for m in s_nb), 3),
                     "hop1_links": len(h1_nb), "betweenness": round(btw[n], 5),
                     "community": comm_of.get(n, -1)})
    df = pd.DataFrame(recs)
    s = df[df.role == "seed"].sort_values(["tier", "degree"], ascending=[True, False])
    nb = df[df.role == "neighbour"].sort_values(
        ["hop", "seed_links", "seed_score_sum", "hop1_links", "betweenness"],
        ascending=[True, False, False, False, False])
    return pd.concat([s, nb]).reset_index(drop=True)


def seed_report(G, seeds_df, module_tiers="AB", max_len=6):
    module = [g for g in seeds_df[seeds_df.tier.isin(list(module_tiers))].gene if g in G]
    rows = []
    for r in seeds_df.itertuples():
        g = r.gene
        if not is_symbol(g):
            rows.append({"gene": g, "tier": r.tier, "in_variant26": r.in_variant26, "in_string_network": False,
                         "string_degree": 0, "dist_to_module": None,
                         "status": "no gene symbol: STRING/OmniPath skipped (OT-only)"})
            continue
        if g not in G or G.degree(g) == 0:
            rows.append({"gene": g, "tier": r.tier, "in_variant26": r.in_variant26, "in_string_network": g in G,
                         "string_degree": 0, "dist_to_module": None, "status": "isolated (no STRING edge at cutoff)"})
            continue
        lens = nx.single_source_shortest_path_length(G, g, cutoff=max_len)
        d = min([lens[m] for m in module if m != g and m in lens], default=None)
        st = "in module" if (g in module) else ("connected to module" if d is not None and d <= 3
                                                 else "weakly/not connected to module")
        rows.append({"gene": g, "tier": r.tier, "in_variant26": r.in_variant26, "in_string_network": True,
                     "string_degree": G.degree(g), "dist_to_module": d, "status": st})
    return pd.DataFrame(rows)


def reactome_enrichment(genes, top=30):
    url = (f"{REACTOME}?pageSize={top}&page=1&sortBy=ENTITIES_FDR&order=ASC"
           f"&resource=TOTAL&pValue=1&includeDisease=true")
    js = http("POST", url, data="\n".join(genes), headers={"Content-Type": "text/plain"}).json()
    out = []
    for p in js.get("pathways", []):
        en = p.get("entities", {})
        out.append({"stId": p.get("stId"), "pathway": p.get("name"), "found": en.get("found"),
                    "total": en.get("total"), "pValue": en.get("pValue"), "fdr": en.get("fdr")})
    return pd.DataFrame(out)


def draw(G, nodes_df, path, label_top=15):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"A": "#c0392b", "B": "#e67e22", "C": "#8e44ad", "D": "#16a085"}
    info = nodes_df.set_index("gene")
    pos = nx.spring_layout(G, weight="score", seed=1, k=0.6)
    fig, ax = plt.subplots(figsize=(14, 11))
    nx.draw_networkx_edges(G, pos, alpha=0.12, width=[0.4 + 2 * G[a][b]["score"] for a, b in G.edges], ax=ax)
    for n in G.nodes:
        r = info.loc[n]
        if r.role == "seed":
            c, size, ec = colors.get(r.tier, "#c0392b"), 140 + 40 * r.degree, "black"
        elif r.hop == 1:
            c, size, ec = "#95a5a6", 80 + 25 * r.degree, "none"
        else:
            c, size, ec = "#d5dbdb", 40 + 10 * r.degree, "none"
        nx.draw_networkx_nodes(G, pos, nodelist=[n], node_color=c, node_size=size, ax=ax,
                               edgecolors=ec, linewidths=0.8)
    lab = set(nodes_df[nodes_df.role == "seed"].gene) | set(nodes_df[nodes_df.role == "neighbour"].head(label_top).gene)
    nx.draw_networkx_labels(G, pos, labels={n: n for n in G.nodes if n in lab}, font_size=8, ax=ax)
    ax.set_title("STRING 2-hop network: seeds (red A, orange B, purple C, teal D), hop1 grey, hop2 light grey")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


# ----------------------------------------------------------------------------- OmniPath (directed)
def load_omnipath(datasets, tf, cache):
    if os.path.exists(cache):
        return pd.read_csv(cache, sep="\t", low_memory=False)
    params = {"genesymbols": "1", "fields": "sources,references,curation_effort,type",
              "datasets": datasets + (",tf_target" if tf else ""), "license": "academic", "format": "tsv"}
    if tf:
        params["dorothea_levels"] = "A,B"
    r = http("GET", OMNIPATH, params=params, timeout=600)
    df = pd.read_csv(io.StringIO(r.text), sep="\t", low_memory=False)
    df.to_csv(cache, sep="\t", index=False)
    return df


def build_digraph(df):
    """Directed, signed gene-level DiGraph. sign: +1 stimulation, -1 inhibition, 0 unknown/both."""
    D = nx.DiGraph()
    col = {c.lower(): c for c in df.columns}
    need = ["source_genesymbol", "target_genesymbol", "is_directed", "is_stimulation", "is_inhibition"]
    if any(n not in col for n in need):
        raise RuntimeError(f"unexpected OmniPath columns: {list(df.columns)}")
    recs = df.to_dict("records")
    for r in recs:
        if int(r.get(col["is_directed"], 0) or 0) != 1:
            continue
        u, v = str(r[col["source_genesymbol"]]), str(r[col["target_genesymbol"]])
        if not (is_symbol(u) and is_symbol(v)) or u == v:
            continue
        stim = int(r.get(col["is_stimulation"], 0) or 0) == 1
        inh = int(r.get(col["is_inhibition"], 0) or 0) == 1
        srcs = str(r.get(col.get("sources", ""), "") or "")
        n_src = len([s for s in srcs.split(";") if s]) if srcs and srcs != "nan" else 0
        ce = r.get(col.get("curation_effort", ""), 0)
        try:
            ce = float(ce)
        except (TypeError, ValueError):
            ce = 0.0
        sup = max(n_src, ce)
        if D.has_edge(u, v):
            e = D[u][v]
            e["stim"] |= stim
            e["inh"] |= inh
            e["support"] = max(e["support"], sup)
        else:
            D.add_edge(u, v, stim=stim, inh=inh, support=sup)
    for u, v, e in D.edges(data=True):
        e["sign"] = 1 if (e["stim"] and not e["inh"]) else (-1 if (e["inh"] and not e["stim"]) else 0)
    return D


ARROW = {1: "->", -1: "-|", 0: "--"}


def expand_directed(D, seeds, levels=2, hub_deg=400, min_support2=2):
    """Per seed, BFS downstream (seed regulates X) and upstream (X regulates seed) up to `levels`."""
    rows = []
    for s in seeds:
        if s not in D:
            continue
        for direction in ("down", "up"):
            seen = {s}
            frontier = [(s, 1, s)]
            for level in range(1, levels + 1):
                nxt = []
                for node, sg, path in frontier:
                    if level > 1 and D.degree(node) > hub_deg:
                        continue
                    nbrs = D.successors(node) if direction == "down" else D.predecessors(node)
                    for m in nbrs:
                        e = D[node][m] if direction == "down" else D[m][node]
                        if level > 1 and e["support"] < min_support2:
                            continue
                        cs = sg * e["sign"] if (sg != 0 and e["sign"] != 0) else 0
                        a = ARROW[e["sign"]]
                        p = f"{path} {a} {m}" if direction == "down" else f"{m} {a} {path}"
                        if m in seen:
                            continue
                        seen.add(m)
                        rows.append({"seed": s, "gene": m, "direction": direction, "level": level,
                                     "sign": cs, "support": e["support"], "path": p})
                        nxt.append((m, cs, p))
                frontier = nxt
    return pd.DataFrame(rows, columns=["seed", "gene", "direction", "level", "sign", "support", "path"])


def auto_rescue(rows, seed_mech, exclude):
    """
    Desired pharmacological action implied by the signed path.
      LoF seed, downstream X: seed -> X (+)  => X falls  => ACT X ;  seed -| X => X rises => INHIB X
      'unknown' seed (possible hypomorph), upstream X: X -> seed => ACT X (boost residual) ; X -| seed => INHIB X
    Genes with equal-level conflicting actions are dropped as ambiguous.
    """
    cand = []
    for r in rows.itertuples():
        if r.gene in exclude or r.sign == 0:
            continue
        mech = seed_mech.get(r.seed, "none")
        if r.direction == "down" and mech == "LoF":
            ctx = "pathway_inferred"
        elif r.direction == "up" and mech == "unknown":
            ctx = "pathway_inferred_hypomorph"
        else:
            continue
        cand.append({"gene": r.gene, "desired_action": "ACT" if r.sign > 0 else "INHIB", "context": ctx,
                     "level": r.level, "support": r.support, "seed": r.seed,
                     "because": f"{r.path} (seed {mech})", "source": f"omnipath_L{r.level}"})
    if not cand:
        return pd.DataFrame(columns=["gene", "desired_action", "context", "level", "support", "seeds", "because", "source"])
    c = pd.DataFrame(cand)
    out = []
    for g, grp in c.groupby("gene"):
        best_level = grp.level.min()
        top = grp[grp.level == best_level]
        if top.desired_action.nunique() > 1:
            continue
        b = top.sort_values("support", ascending=False).iloc[0]
        out.append({"gene": g, "desired_action": b.desired_action, "context": b.context, "level": int(best_level),
                    "support": b.support, "seeds": ";".join(sorted(set(grp.seed))),
                    "because": b.because, "source": b.source})
    return pd.DataFrame(out).sort_values(["level", "support"], ascending=[True, False]).reset_index(drop=True)


def directed_files(D, rows, out):
    rows.to_csv(f"{out}/directed_nodes.csv", index=False)
    keep = set(rows.gene) | set(rows.seed)
    H = D.subgraph(keep).copy()
    edges = [{"source": u, "target": v, "sign": e["sign"], "support": e["support"]} for u, v, e in H.edges(data=True)]
    pd.DataFrame(edges).to_csv(f"{out}/directed_edges.csv", index=False)
    G2 = nx.DiGraph()
    for u, v, e in H.edges(data=True):
        G2.add_edge(u, v, sign=int(e["sign"]), support=float(e["support"]))
    nx.write_graphml(G2, f"{out}/directed.graphml")
    rel = {1: "activates", -1: "inhibits", 0: "regulates"}
    with open(f"{out}/directed.sif", "w") as fh:
        for u, v, e in H.edges(data=True):
            fh.write(f"{u}\t{rel[e['sign']]}\t{v}\n")


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="seeds.csv")
    ap.add_argument("--out", default="out")
    ap.add_argument("--score", type=float, default=0.7, help="STRING cutoff for hop 1 (0-1)")
    ap.add_argument("--add-nodes", type=int, default=40, help="extra interactors STRING adds (hop 1)")
    ap.add_argument("--score2", type=float, default=0.9, help="STRING cutoff for any edge touching a hop-2 node")
    ap.add_argument("--hop2-from", type=int, default=50, help="expand the best N hop-1 nodes to hop 2")
    ap.add_argument("--hop2-limit", type=int, default=10, help="max hop-2 partners per expanded node")
    ap.add_argument("--no-hop2", action="store_true")
    ap.add_argument("--levels", type=int, default=2, help="directed levels up and down (OmniPath)")
    ap.add_argument("--hub-degree", type=int, default=400)
    ap.add_argument("--min-support2", type=float, default=2, help="min sources for a level-2 directed edge")
    ap.add_argument("--op-datasets", default="omnipath,pathwayextra,kinaseextra")
    ap.add_argument("--tf", action="store_true", help="add TF->target edges (tf_target, DoRothEA A/B)")
    ap.add_argument("--no-omnipath", action="store_true")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    seeds = pd.read_csv(a.seeds)
    sym = [g for g in seeds.gene if is_symbol(g)]
    skipped = [g for g in seeds.gene if not is_symbol(g)]
    print(f"{len(seeds)} seed rows; {len(sym)} symbols queried; no-symbol (OT-only): {', '.join(skipped) or '-'}")

    # ---- STRING hop 1
    rows1 = string_network(sym, a.score, a.add_nodes)
    if not rows1:
        sys.exit("STRING returned no edges - lower --score")
    G = finish_graph(add_rows(nx.Graph(), rows1))
    seedset = set(sym)
    hop = {n: (0 if n in seedset else 1) for n in G.nodes}
    print(f"hop1: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # ---- STRING hop 2 (higher cutoff)
    if not a.no_hop2:
        expand = rank_hop1(G, seedset)[:a.hop2_from]
        rows2 = string_partners(expand, a.score2, a.hop2_limit) if expand else []
        add_rows(G, rows2)
        for n in G.nodes:
            hop.setdefault(n, 2)
        allnodes = list(G.nodes)
        rows3 = string_network(allnodes, a.score2, 0)          # edges among all nodes at >= score2
        add_rows(G, rows3)
        for n in G.nodes:
            hop.setdefault(n, 2)
        filter_hop2(G, hop, a.score2)
        finish_graph(G)
        print(f"hop2: {sum(v == 2 for v in hop.values())} extra nodes (edge score >= {a.score2}); "
              f"total {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    nodes = node_table(G, seeds, hop)
    edges = pd.DataFrame([{"gene_a": u, "gene_b": v, "score": d["score"]} for u, v, d in G.edges(data=True)])
    edges.to_csv(f"{a.out}/edges.csv", index=False)
    nodes.to_csv(f"{a.out}/nodes.csv", index=False)
    H = G.copy()
    ninfo = nodes.set_index("gene")
    for n in H.nodes:
        r = ninfo.loc[n]
        H.nodes[n].update(role=r.role, tier=str(r.tier), hop=int(r.hop), community=int(r.community),
                          degree=int(r.degree), seed_links=int(r.seed_links))
    nx.write_graphml(H, f"{a.out}/network.graphml")
    with open(f"{a.out}/network.sif", "w") as fh:
        for u, v, d in G.edges(data=True):
            fh.write(f"{u}\tppi\t{v}\n")
    draw(G, nodes, f"{a.out}/network.png")

    rep = seed_report(G, seeds)
    # ---- OmniPath directed
    if not a.no_omnipath:
        try:
            op = load_omnipath(a.op_datasets, a.tf, f"{a.out}/omnipath_raw.tsv")
            D = build_digraph(op)
            print(f"OmniPath: {D.number_of_nodes()} genes, {D.number_of_edges()} directed edges")
            drows = expand_directed(D, sym, a.levels, a.hub_degree, a.min_support2)
            mech = dict(zip(seeds.gene, seeds.variant_mechanism))
            ar = auto_rescue(drows, mech, exclude=seedset)
            ar.to_csv(f"{a.out}/auto_rescue_targets.csv", index=False)
            directed_files(D, drows, a.out)
            cnt = drows.groupby(["seed", "direction", "level"]).size().unstack(["direction", "level"]).fillna(0).astype(int)
            cnt.columns = [f"{d}_L{l}" for d, l in cnt.columns]
            rep = rep.merge(cnt, left_on="gene", right_index=True, how="left")
            print(f"directed: {len(drows)} gene-level rows; auto rescue targets: {len(ar)}")
        except Exception as e:  # noqa: BLE001
            print("OmniPath layer skipped:", e)
    rep.to_csv(f"{a.out}/seed_report.csv", index=False)
    print("\nSeed connectivity:\n", rep[["gene", "tier", "string_degree", "dist_to_module", "status"]].to_string(index=False))

    try:
        rx = reactome_enrichment(list(G.nodes))
        rx.to_csv(f"{a.out}/reactome_enrichment.csv", index=False)
        print("Top Reactome pathways:\n", rx.head(10).to_string(index=False))
    except Exception as e:  # noqa: BLE001
        print("Reactome step skipped:", e)

    print(f"\nNetwork: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print("\nTop NEW genes to search (hop1 first, ranked by links to your seeds):")
    print(nodes[nodes.role == "neighbour"].head(25)[["gene", "hop", "seed_links", "seed_score_sum", "betweenness", "community"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
