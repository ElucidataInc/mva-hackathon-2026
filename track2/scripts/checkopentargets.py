
import argparse, csv, json, re, sys, time, urllib.request, urllib.error

API = "https://api.platform.opentargets.org/api/v4/graphql"
TERMS = re.compile(r"aneuploid|mosaic variegated|chromosom|microcephal|premature chromatid|"
                   r"cancer predisposition|wilms|mitotic|segregation", re.I)
VAR_RE = re.compile(r"\(chr([0-9XYM]+):(\d+) ([ACGT]+)>([ACGT]+)\)")


def gql(query, variables=None, retries=3):
    body = json.dumps({"query": query, "variables": variables or {}}).encode()
    req = urllib.request.Request(API, body, {"Content-Type": "application/json",
                                             "Accept": "application/json",
                                             "User-Agent": "top10-helper-check/1.0"})
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                out = json.load(r)
            return out
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code}: {e.read()[:300].decode(errors='replace')}"
            if e.code < 500:
                break
        except Exception as e:  # network error, timeout
            last = str(e)
        time.sleep(2 * (i + 1))
    return {"errors": [{"message": last}]}


def read_table(path):
    if path.lower().endswith((".xlsx", ".xlsm")):
        from openpyxl import load_workbook
        ws = load_workbook(path, read_only=True, data_only=True).active
        rows = [["" if c is None else str(c) for c in r] for r in ws.iter_rows(values_only=True)]
    else:
        with open(path, newline="", encoding="utf-8") as fh:
            delim = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
            rows = list(csv.reader(fh, delimiter=delim))
    hdr = rows[0]
    for need in ("Evidence count (0-6)", "Variant", "Gene(s)"):
        if need not in hdr:
            sys.exit(f"Column '{need}' not found. Columns are: {hdr}")
    ix = {h: hdr.index(h) for h in ("Evidence count (0-6)", "Variant", "Gene(s)")}
    out = []
    for n, r in enumerate(rows[1:], start=2):
        if not r or not r[ix["Evidence count (0-6)"]].strip().isdigit():
            continue  # footnotes / blank lines
        out.append(dict(row=n, ev=int(r[ix["Evidence count (0-6)"]]),
                        variant=r[ix["Variant"]].strip(), genes=r[ix["Gene(s)"]].strip()))
    return out


def gene_symbols(s):
    return [g.strip() for g in s.split(",") if g.strip() and not g.strip().startswith("ENSG")]


_target_cache = {}


def check_gene(sym):
    if sym in _target_cache:
        return _target_cache[sym]
    res = {"symbol": sym, "ensembl": "", "n_assoc": "", "mva_hits": "", "top5": "", "error": ""}
    s = gql("""query($q:String!){search(queryString:$q, entityNames:["target"], page:{index:0,size:5}){hits{id name entity}}}""",
            {"q": sym})
    if "errors" in s and not s.get("data"):
        res["error"] = "search: " + s["errors"][0]["message"]; _target_cache[sym] = res; return res
    hits = [h for h in (s["data"]["search"]["hits"] or []) if h["entity"] == "target"]
    exact = [h for h in hits if (h["name"] or "").upper() == sym.upper()]
    pick = (exact or hits or [None])[0]
    if not pick:
        res["error"] = "target not found"; _target_cache[sym] = res; return res
    res["ensembl"] = pick["id"]
    d = gql("""query($id:String!){target(ensemblId:$id){approvedSymbol
                 associatedDiseases(page:{index:0,size:300}){count rows{score disease{name}}}}}""",
            {"id": pick["id"]})
    if "errors" in d and not d.get("data"):
        res["error"] = "associations: " + d["errors"][0]["message"]
    else:
        t = d["data"]["target"]
        rows = t["associatedDiseases"]["rows"]
        res["n_assoc"] = t["associatedDiseases"]["count"]
        rows = sorted(rows, key=lambda r: -r["score"])
        res["top5"] = "; ".join(f'{r["disease"]["name"]} ({r["score"]:.2f})' for r in rows[:5])
        res["mva_hits"] = "; ".join(f'{r["disease"]["name"]} ({r["score"]:.2f})'
                                    for r in rows if TERMS.search(r["disease"]["name"])) or "none"
    time.sleep(0.2)
    _target_cache[sym] = res
    return res


_schema = {}


def type_fields(tname):
    """{field_name: base_type_name} for a GraphQL type, via introspection (cached)."""
    if tname in _schema:
        return _schema[tname]
    out = gql("""query($t:String!){__type(name:$t){fields{name args{name} type{name kind ofType{name kind ofType{name kind ofType{name}}}}}}}""",
              {"t": tname})
    f = {}
    t = (out.get("data") or {}).get("__type")
    for fld in (t or {}).get("fields") or []:
        ty = fld["type"]
        while ty and not ty.get("name") and ty.get("ofType"):
            ty = ty["ofType"]
        f[fld["name"]] = {"type": (ty or {}).get("name"), "args": [a["name"] for a in fld.get("args") or []]}
    _schema[tname] = f
    return f


# what we would like from a Variant: name -> list of sub-fields (None = scalar)
WISH = [("id", None), ("variantId", None), ("rsIds", None),
        ("mostSevereConsequence", ["label"]),
        ("alleleFrequencies", ["populationName", "alleleFrequency"]),
        ("variantEffect", ["method", "assessment", "score"])]


def build_variant_query():
    """Build the query from fields the API really has; returns (query, arg_name, notes)."""
    vf = type_fields("Variant")
    q = type_fields("Query").get("variant")
    if not vf or not q:
        return None, None, "could not read API schema (Variant type or variant query missing)"
    arg = q["args"][0] if q["args"] else "variantId"
    parts = []
    for name, subs in WISH:
        if name not in vf:
            continue
        if subs is None:
            parts.append(name)
        else:
            sf = type_fields(vf[name]["type"])
            ok = [x for x in subs if x in sf] if sf else []
            if ok:
                parts.append(f"{name}{{{' '.join(ok)}}}")
    if not parts:
        return None, None, "no usable Variant fields; Variant has: " + ",".join(vf)
    return "query($v:String!){variant(%s:$v){%s}}" % (arg, " ".join(parts)), arg, ""


_vq = None


def check_variant(vid):
    global _vq
    res = {"found": "no", "rsids": "", "consequence": "", "gnomad": "", "effects": "", "error": ""}
    if _vq is None:
        _vq = build_variant_query()
        print("  variant query:", _vq[0] or _vq[2])
    q, _, why = _vq
    if not q:
        res["found"] = "error"; res["error"] = why; return res
    out = gql(q, {"v": vid})
    if "errors" in out and not (out.get("data") or {}).get("variant"):
        res["found"] = "error"; res["error"] = out["errors"][0]["message"][:200]; return res
    v = (out.get("data") or {}).get("variant")
    if not v:
        return res  # API answered, variant genuinely not in Open Targets
    res["found"] = "yes"
    res["rsids"] = ",".join(v.get("rsIds") or [])
    res["consequence"] = (v.get("mostSevereConsequence") or {}).get("label", "")
    af = v.get("alleleFrequencies") or []
    res["gnomad"] = "; ".join(f'{a["populationName"]}={a["alleleFrequency"]:.2e}' for a in af[:8])
    res["effects"] = "; ".join(f'{e.get("method")}:{e.get("assessment") or ""}:{e.get("score")}'
                               for e in (v.get("variantEffect") or []))
    time.sleep(0.2)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--min-evidence", type=int, default=3)
    ap.add_argument("--out", default="opentargets_check.tsv")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    sel = [r for r in read_table(a.table) if r["ev"] >= a.min_evidence]
    print(f"{len(sel)} variants with evidence >= {a.min_evidence}")
    for r in sel:
        m = VAR_RE.search(r["variant"])
        r["vid"] = f"{m.group(1)}_{m.group(2)}_{m.group(3)}_{m.group(4)}" if m and "<" not in r["variant"] else ""
        r["syms"] = gene_symbols(r["genes"])
        kind = "variant+gene" if r["vid"] else "gene only"
        print(f"  row {r['row']:>3}  ev={r['ev']}  {kind:12}  {r['variant'][:70]}")
    if a.dry_run:
        print("\n--dry-run: nothing queried."); return

    cols = ["sheet_row", "evidence", "variant", "ot_variant_id", "ot_variant_found", "rsIDs",
            "most_severe_consequence", "gnomAD_AF", "variant_effect_predictions", "gene",
            "ensembl_id", "n_disease_associations", "MVA_related_associations",
            "top5_associations", "notes"]
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(cols)
        for r in sel:
            vr = check_variant(r["vid"]) if r["vid"] else None
            note = []
            if not r["vid"]:
                note.append("no variant-level lookup (SV/CNV or alleles not given)")
            if vr and vr["error"]:
                note.append(vr["error"])
            genes = r["syms"] or [""]
            for g in genes:
                gr = check_gene(g) if g else {"ensembl": "", "n_assoc": "", "mva_hits": "", "top5": "", "error": "no gene symbol"}
                n = note + ([gr["error"]] if gr["error"] else [])
                w.writerow([r["row"], r["ev"], r["variant"], r["vid"], vr["found"] if vr else "n/a",
                            vr["rsids"] if vr else "", vr["consequence"] if vr else "",
                            vr["gnomad"] if vr else "", vr["effects"] if vr else "", g, gr["ensembl"],
                            gr["n_assoc"], gr["mva_hits"], gr["top5"], " | ".join(n)])
            print(f"done row {r['row']}: {r['variant'][:55]}")
    print(f"\nWrote {a.out}")
    errs = [s for s, v in _target_cache.items() if v["error"]]
    if errs:
        print("Gene lookups with problems:", ", ".join(errs))


if __name__ == "__main__":
    main()
