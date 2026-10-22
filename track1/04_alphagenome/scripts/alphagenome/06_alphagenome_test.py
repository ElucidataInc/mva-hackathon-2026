#!/usr/bin/env python3
"""Connectivity test: score the public AlphaGenome documentation example variant (NOT patient data)."""
import os
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

OUT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/alphagenome/v2_team11/alphagenome_model/test"))
model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])

variant = genome.Variant(chromosome="chr22", position=36201698,
                         reference_bases="A", alternate_bases="C")   # public docs example
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
scores = model.score_variant(interval=interval, variant=variant, variant_scorers=scorers)
df = variant_scorers.tidy_scores([scores])
df.to_csv(f"{OUT}/test_example_variant_scores.tsv", sep="\t", index=False)

print("rows:", len(df))
print(df["output_type"].value_counts().to_string())
