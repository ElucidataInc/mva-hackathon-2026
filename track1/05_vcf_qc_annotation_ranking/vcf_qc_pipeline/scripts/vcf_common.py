"""Shared helpers for the VCF QC pipeline: opening VCFs, header parsing,
caller / build / annotation detection, allele classification, FASTA access."""

import gzip
import io
import json
import os
import re

# ----------------------------------------------------------------------------
# File access
# ----------------------------------------------------------------------------

def open_vcf(path):
    """Open plain / gzip / bgzip VCF as text."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def write_json(obj, path):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=_json_default)


def _json_default(o):
    try:
        import numpy as np
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
    except ImportError:
        pass
    if isinstance(o, set):
        return sorted(o)
    raise TypeError(str(type(o)))


# ----------------------------------------------------------------------------
# Header parsing
# ----------------------------------------------------------------------------

_META_RE = re.compile(r'^##(\w+)=<(.*)>$')
_KV_RE = re.compile(r'(\w+)=("(?:[^"\\]|\\.)*"|[^,]*)')


def _parse_struct(body):
    d = {}
    for k, v in _KV_RE.findall(body):
        d[k] = v.strip('"')
    return d


def parse_header(path):
    """Return dict with header metadata and the sample list."""
    hdr = {
        "fileformat": None, "info": {}, "format": {}, "filter": {}, "contigs": [],
        "alt": {}, "other": [], "samples": [], "n_header_lines": 0, "sources": [],
        "commands": [], "reference": None, "tumor_sample": None, "normal_sample": None,
    }
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("##"):
                hdr["n_header_lines"] += 1
                line = line.rstrip("\n")
                m = _META_RE.match(line)
                if m and m.group(1) in ("INFO", "FORMAT", "FILTER", "contig", "ALT"):
                    st = _parse_struct(m.group(2))
                    key = m.group(1)
                    if key == "contig":
                        hdr["contigs"].append({"ID": st.get("ID"),
                                               "length": int(st["length"]) if st.get("length", "").isdigit() else None,
                                               "assembly": st.get("assembly")})
                    else:
                        hdr[key.lower()][st.get("ID")] = st
                    continue
                k, _, v = line[2:].partition("=")
                if k == "fileformat":
                    hdr["fileformat"] = v
                elif k == "source":
                    hdr["sources"].append(v)
                elif k == "reference":
                    hdr["reference"] = v
                elif k == "tumor_sample":
                    hdr["tumor_sample"] = v
                elif k == "normal_sample":
                    hdr["normal_sample"] = v
                elif (k.endswith("Command") or k.endswith("CommandLine") or k in ("commandline", "DeepVariant_version",
                      "GATKCommandLine", "DRAGENCommandLine", "fileDate", "FILTER_NOTE")
                      or "version" in k.lower()):
                    hdr["commands"].append({"key": k, "value": v[:600]})
                else:
                    hdr["other"].append({"key": k, "value": v[:300]})
            elif line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                hdr["samples"] = cols[9:]
                break
            else:
                break
    return hdr


# ----------------------------------------------------------------------------
# Provenance detection
# ----------------------------------------------------------------------------

CALLER_PATTERNS = [
    ("GATK Mutect2", [r"Mutect2", r"\bTLOD\b"]),
    ("GATK HaplotypeCaller", [r"HaplotypeCaller"]),
    ("GATK GenotypeGVCFs", [r"GenotypeGVCFs"]),
    ("GATK FilterMutectCalls", [r"FilterMutectCalls"]),
    ("GATK VQSR/ApplyVQSR", [r"ApplyVQSR|VariantRecalibrator|VQSLOD"]),
    ("DeepVariant", [r"DeepVariant"]),
    ("GLnexus", [r"GLnexus"]),
    ("Illumina DRAGEN", [r"DRAGEN"]),
    ("Strelka2", [r"strelka"]),
    ("Manta (SV)", [r"[Mm]anta"]),
    ("FreeBayes", [r"freeBayes|freebayes"]),
    ("bcftools call", [r"bcftools_call"]),
    ("bcftools mpileup", [r"bcftools_mpileup"]),
    ("VarScan", [r"VarScan"]),
    ("Octopus", [r"octopus"]),
    ("Platypus", [r"Platypus"]),
    ("Clair3", [r"Clair3|clair3"]),
    ("Sentieon", [r"Sentieon|sentieon"]),
    ("Parabricks", [r"[Pp]arabricks"]),
    ("Ensembl VEP", [r"\bVEP\b|VEP=|\bCSQ\b"]),
    ("SnpEff", [r"SnpEff|\bANN\b"]),
    ("ANNOVAR", [r"ANNOVAR"]),
]


def detect_callers(hdr, header_text):
    found = []
    blob = header_text
    for name, pats in CALLER_PATTERNS:
        if any(re.search(p, blob) for p in pats):
            found.append(name)
    return found


def read_header_text(path, limit_chars=2_000_000):
    out = []
    n = 0
    with open_vcf(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            out.append(line)
            n += len(line)
            if n > limit_chars:
                break
    return "".join(out)


BUILD_LENGTHS = {
    # chr1 length -> build
    248956422: "GRCh38/hg38", 249250621: "GRCh37/hg19", 247249719: "NCBI36/hg18",
    195471971: "GRCm39/mm39 (mouse, chr1)", 197195432: "GRCm38/mm10 (mouse, chr1)",
    248387328: "T2T-CHM13v2.0",
}


def guess_build(hdr):
    for c in hdr["contigs"]:
        if c["ID"] in ("1", "chr1") and c["length"]:
            return BUILD_LENGTHS.get(c["length"], f"unknown (chr1 length {c['length']})")
    ref = (hdr.get("reference") or "").lower()
    for tag, b in (("38", "GRCh38/hg38 (from ##reference)"), ("37", "GRCh37 (from ##reference)"),
                   ("hg19", "GRCh37/hg19 (from ##reference)"), ("chm13", "T2T-CHM13")):
        if tag in ref:
            return b
    return "unknown (no contig lengths in header)"


def detect_mode(hdr, header_text):
    """germline / somatic / gvcf / sites-only"""
    if not hdr["samples"]:
        return "sites-only"
    if "TLOD" in hdr["info"] or hdr["tumor_sample"] or re.search(r"Mutect2|Strelka.*somatic|somatic", header_text):
        return "somatic"
    if "<NON_REF>" in header_text or "<*>" in header_text and "END" in hdr["info"]:
        return "gvcf"
    return "germline"


def detect_annotations(hdr):
    info = hdr["info"]
    ann = {}
    if "CSQ" in info:
        desc = info["CSQ"].get("Description", "")
        m = re.search(r"Format:\s*([^\"]+)", desc)
        ann["VEP_CSQ"] = m.group(1).strip().split("|") if m else []
    if "ANN" in info:
        ann["SnpEff_ANN"] = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID",
                             "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c",
                             "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance", "Errors"]
    if "BCSQ" in info:
        ann["bcftools_csq"] = ["Consequence", "gene", "transcript", "biotype", "strand", "aa_change", "dna_change"]
    gnomad = [k for k in info if re.search(r"gnomad", k, re.I)]
    if gnomad:
        ann["gnomAD_fields"] = gnomad
    clin = [k for k in info if k.startswith("CLN")]
    if clin:
        ann["ClinVar_fields"] = clin
    if "DB" in info:
        ann["dbSNP_flag"] = ["DB"]
    if "RS" in info:
        ann["dbSNP_RS"] = ["RS"]
    func = [k for k in info if k.startswith(("Func.", "ExonicFunc.", "Gene."))]
    if func:
        ann["ANNOVAR_fields"] = func
    return ann


# ----------------------------------------------------------------------------
# Allele classification
# ----------------------------------------------------------------------------

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def classify_allele(ref, alt):
    """Return (class, indel_length). indel_length>0 insertion, <0 deletion."""
    if alt in ("*",):
        return "spanning_deletion", 0
    if alt in (".", "<NON_REF>", "<*>"):
        return "ref_block", 0
    if alt.startswith("<"):
        return "symbolic_SV:" + alt.strip("<>").split(":")[0], 0
    if "[" in alt or "]" in alt:
        return "breakend", 0
    lr, la = len(ref), len(alt)
    if lr == la:
        if lr == 1:
            return "SNV", 0
        # trim common prefix / suffix to see whether it is really an SNV
        diff = sum(1 for a, b in zip(ref, alt) if a != b)
        if diff == 1:
            return "SNV", 0
        return "MNV", 0
    d = la - lr
    if ref[0] == alt[0] and (alt.startswith(ref) or ref.startswith(alt)):
        return ("insertion" if d > 0 else "deletion"), d
    # left-aligned but not anchored -> complex indel
    return "complex_indel", d


def snv_bases(ref, alt):
    for a, b in zip(ref, alt):
        if a != b:
            return a.upper(), b.upper()
    return ref[0].upper(), alt[0].upper()


def chrom_sort_key(c):
    c = str(c)
    s = c[3:] if c.lower().startswith("chr") else c
    if s.isdigit():
        return (0, int(s), "")
    order = {"X": 1, "Y": 2, "M": 3, "MT": 3}
    return (1, order.get(s.upper(), 9), s)


def norm_chrom(c):
    c = str(c)
    s = c[3:] if c.lower().startswith("chr") else c
    return "MT" if s.upper() == "M" else s.upper()


def is_autosome(c):
    return norm_chrom(c).isdigit()


# GRCh38 / GRCh37 pseudo-autosomal regions on X (1-based inclusive)
PAR = {
    "GRCh38": [(10001, 2781479), (155701383, 156030895)],
    "GRCh37": [(60001, 2699520), (154931044, 155260560)],
}


def in_par(pos, build):
    key = "GRCh38" if "38" in build else ("GRCh37" if "37" in build or "19" in build else "GRCh38")
    return any(a <= pos <= b for a, b in PAR[key])


# ----------------------------------------------------------------------------
# Minimal indexed FASTA reader (no external deps)
# ----------------------------------------------------------------------------

class Fasta:
    def __init__(self, path):
        self.path = path
        fai = path + ".fai"
        if not os.path.exists(fai):
            self._build_fai(fai)
        self.idx = {}
        with open(fai) as fh:
            for line in fh:
                name, length, offset, lb, lw = line.split("\t")[:5]
                self.idx[name] = (int(length), int(offset), int(lb), int(lw))
        self.fh = open(path, "rb")
        self._alias = {}
        for n in self.idx:
            self._alias[n] = n
            alt = n[3:] if n.startswith("chr") else "chr" + n
            self._alias.setdefault(alt, n)
        self._alias.setdefault("chrM", self._alias.get("MT"))
        self._alias.setdefault("MT", self._alias.get("chrM"))

    def _build_fai(self, fai):
        if self.path.endswith(".gz"):
            raise RuntimeError("Gzipped FASTA without .fai; run `samtools faidx` first")
        entries = []
        with open(self.path, "rb") as fh:
            name = None
            offset = 0
            length = lb = lw = 0
            pos = 0
            for line in fh:
                if line.startswith(b">"):
                    if name:
                        entries.append((name, length, offset, lb, lw))
                    name = line[1:].split()[0].decode()
                    length = lb = lw = 0
                    offset = pos + len(line)
                else:
                    if lb == 0:
                        lb = len(line.rstrip(b"\r\n"))
                        lw = len(line)
                    length += len(line.rstrip(b"\r\n"))
                pos += len(line)
            if name:
                entries.append((name, length, offset, lb, lw))
        with open(fai, "w") as out:
            for e in entries:
                out.write("\t".join(map(str, e)) + "\n")

    def has(self, chrom):
        return self._alias.get(chrom) in self.idx

    def fetch(self, chrom, start1, end1):
        """1-based inclusive."""
        name = self._alias.get(chrom)
        if name is None:
            return None
        length, offset, lb, lw = self.idx[name]
        if start1 < 1 or end1 > length:
            return None
        s0 = start1 - 1
        e0 = end1
        fs = offset + (s0 // lb) * lw + s0 % lb
        fe = offset + ((e0 - 1) // lb) * lw + (e0 - 1) % lb + 1
        self.fh.seek(fs)
        raw = self.fh.read(fe - fs)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()


# ----------------------------------------------------------------------------
# Annotation field helpers
# ----------------------------------------------------------------------------

def parse_info(s):
    d = {}
    if s == "." or not s:
        return d
    for kv in s.split(";"):
        k, eq, v = kv.partition("=")
        d[k] = v if eq else True
    return d


def to_float(x):
    try:
        if x in (None, ".", "", True):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def first_float(x):
    if x is None or x is True:
        return None
    return to_float(str(x).split(",")[0])


def read_sex_file(path):
    """TSV/CSV: sample<TAB>sex with sex in {M,F,male,female,1,2}. Also accepts a PED."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = re.split(r"[\t, ]+", line.strip())
            if len(parts) >= 6:          # PED: fam iid pat mat sex pheno
                sid, sex = parts[1], parts[4]
            elif len(parts) >= 2:
                sid, sex = parts[0], parts[1]
            else:
                continue
            sex = sex.strip().lower()
            out[sid] = {"1": "M", "m": "M", "male": "M", "2": "F", "f": "F", "female": "F"}.get(sex, "U")
    return out


def read_ped(path):
    """Return list of trios (child, father, mother) present in a PED file."""
    trios = []
    if not path or not os.path.exists(path):
        return trios
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = re.split(r"\s+", line.strip())
            if len(p) < 4:
                continue
            fam, iid, pat, mat = p[:4]
            if pat not in ("0", ".") and mat not in ("0", "."):
                trios.append((iid, pat, mat))
    return trios
