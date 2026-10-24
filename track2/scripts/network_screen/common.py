"""Shared helpers: HTTP with retries, context config, clinical-stage parsing."""
import json
import os
import re
import time

import requests

DEFAULT_CONTEXT = {
    "affected_tissues": ["brain", "cerebral cortex", "cerebellum", "hippocampus", "bone marrow", "blood",
                         "lymph node", "spleen", "kidney", "fibroblast"],
    "related_disease_regex": ("microcephal|aneuploid|premature chromatid|leuk|lymphoma|wilms|nephroblastoma|"
                              "rhabdomyosarcoma|neuroblastoma|cancer predisposition|noonan|rasopath|"
                              "primordial dwarfism|seckel|fanconi|developmental delay|intellectual disability|"
                              "growth retardation|atelis|mulibrey"),
    "toxicity_regex": "toxic|safety|adverse|death|fatal|hepatotox|cardiotox|side.?effect|tolerab|SAE|dose.limiting",
    "futility_regex": "futil|lack of (clinical )?efficacy|no efficacy|ineffective|not effective|did not meet|"
                      "insufficient efficacy|failed to (meet|show)|interim analysis",
    "terminated_status_regex": "terminat|withdrawn|suspend",
    "active_status_regex": "complet|recruit|active|enrol|approved|not yet",
    "serious_ae_regex": ("death|cardiac|arrhythm|QT|hepat|liver|neutropen|pancytopen|agranulocyt|leuk|marrow|"
                         "seizure|encephal|suicid|renal failure|embol|thrombo|sepsis|myelodysplas|neoplasm"),
    "target_safety_regex": ("cardi|arrhythm|QT|heart|neuro|seizure|cerebel|hepat|liver|marrow|haemat|hemat|"
                            "leuk|neutropen|renal|kidney|death|lethal|carcino|tumou?r|embryo|teratogen|"
                            "genotox|aneuploid|chromosom"),
    "ko_lethal_regex": "lethal|embryonic death|prenatal death|perinatal|postnatal death|premature death",
    "ko_relevant_regex": "neuro|brain|microceph|seizure|behavio|growth|body size|body weight|anemia|blood|"
                         "bone marrow|thymus|leuk|neoplasm|tumou?r|cardi|chromosom|aneuploid|cell death",
    "expressed_min_rna_level": 1,
    "faers_top_n": 5,
}


def load_context(path="context.json"):
    ctx = dict(DEFAULT_CONTEXT)
    if path and os.path.exists(path):
        with open(path) as fh:
            ctx.update(json.load(fh))
    return ctx


def http(method, url, retries=3, timeout=180, **kw):
    last = None
    for i in range(retries):
        try:
            r = requests.request(method, url, timeout=timeout, **kw)
            r.raise_for_status()
            return r
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"{method} failed {url}: {last}")


def stage_to_num(x):
    """Clinical stage / phase -> float. Handles 4, 'PHASE_3', 'PHASE_2_3', 'APPROVAL', 'EARLY_PHASE_1'."""
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).upper()
    if "APPROV" in s:
        return 4.0
    if "EARLY" in s:
        return 0.5
    nums = [int(n) for n in re.findall(r"\d", s)]
    if nums:
        return float(max(nums))
    return None


def is_symbol(g):
    """Gene strings we can send to STRING / UniProt / OmniPath."""
    g = str(g)
    return bool(g) and not g.upper().startswith("ENSG") and "_" not in g
