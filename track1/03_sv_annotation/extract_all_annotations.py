import pandas as pd

input_file = "combined.annotsv.tsv"
output_file = "combined.MANE_gnomAD_AnnotSV.tsv"

df = pd.read_csv(input_file, sep="\t", dtype=str).fillna("")

if "INFO" not in df.columns:
    raise SystemExit("ERROR: INFO column not found.")

fields = [
    "MANE_GENE",
    "Max_AF",
    "Max_Het",
    "Max_HomAlt",
    "Max_PopMax_AF",
    "gnomAD_Count",
    "Best_gnomAD_ID",
    "Best_gnomAD_OFP",
    "Best_gnomAD_AF",
    "Best_gnomAD_Het",
    "Best_gnomAD_HomAlt",
    "Best_gnomAD_PopMax_AF",
]

def extract_info(info, field):
    for item in info.split(";"):
        if item.startswith(field + "="):
            return item.split("=", 1)[1]
    return ""

for field in fields:
    df[field] = df["INFO"].apply(lambda x: extract_info(x, field))

df.to_csv(output_file, sep="\t", index=False)

print("Output:", output_file)
print("Rows:", len(df))
print("Columns:", len(df.columns))
print("MANE annotated rows:", (df["MANE_GENE"] != "").sum())
print("Rows with gnomAD matches:",
      (pd.to_numeric(df["gnomAD_Count"], errors="coerce").fillna(0) > 0).sum())
