#!/usr/bin/env python3
"""08a: list AlphaGenome track metadata (public, no patient data) to choose ontology terms for predict_variant."""
import os, re, pandas as pd
from alphagenome.models import dna_client
OUT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/alphagenome/v2_team11/alphagenome_model/track_metadata"))
model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
md = model.output_metadata(dna_client.Organism.HOMO_SAPIENS)
parts = []
for ot in dna_client.OutputType:
    df = md.get(ot) if hasattr(md, "get") else getattr(md, ot.name.lower(), None)
    if df is None or len(df) == 0: continue
    df = df.copy(); df["output_type"] = ot.name; parts.append(df)
allmd = pd.concat(parts, ignore_index=True)
allmd.to_csv(f"{OUT}/track_metadata_all.tsv", sep="\t", index=False)
print("== tracks per output type =="); print(allmd.groupby("output_type").size().to_string())
rna = allmd[allmd.output_type == "RNA_SEQ"]
print("\n== RNA-seq biosample types =="); print(rna.biosample_type.value_counts().to_string())
pat = re.compile(r"stem|neural|neuro|progenitor|fetal|embryo|brain|fibroblast|lymph|T-cell|B cell|blood|K562|H1|H9", re.I)
hits = rna[rna.biosample_name.astype(str).str.contains(pat)][["ontology_curie", "biosample_name", "biosample_type"]].drop_duplicates()
hits.to_csv(f"{OUT}/candidate_ontology_terms.tsv", sep="\t", index=False)
print(f"\n== candidate RNA-seq biosamples for MVA-relevant tissues ({len(hits)}) ==")
print(hits.head(60).to_string(index=False))
