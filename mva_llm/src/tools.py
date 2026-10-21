"""Local evidence tools used by the automatic workflow."""

import functools
import json
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

GENE_MOL = str(Path(__file__).resolve().parents[2] / "gene_mol")
sys.path.insert(0, str(Path(GENE_MOL) / "src"))
import model_store
from fingerprints import smiles_to_csr
SERVING_DIR = os.path.join(GENE_MOL, "models/serving_v2")
KG = str(Path(__file__).resolve().parents[1] / "data/kg")
CLINVAR_FULL = os.path.join(KG, "clinvar_full.parquet")

@functools.lru_cache(maxsize=1)
def _ext_calibration():
    import json as _json
    with open(os.path.join(SERVING_DIR, "unified.json")) as fh:
        d = _json.load(fh)
    bands = {int(k): (np.asarray(v["z"]), np.asarray(v["p"]))
             for k, v in d["bands"].items()}
    return bands, d["band_edges_log10_actives"], d.get("ext_gate_p", 0.05)

def _ext_precision(n_active, z_value):
    """Calibrated P(this target is a measured active) for one extended hit."""
    bands, edges, _ = _ext_calibration()
    b = int(np.digitize(np.log10(max(float(n_active), 1.0)), edges))
    xs, ys = bands.get(b, bands[max(bands)])
    return float(np.interp(float(z_value), xs, ys))

@functools.lru_cache(maxsize=1)
def _serving():
    man = pd.read_csv(f"{SERVING_DIR}/manifest.csv")
    man = man[man.list.isin(["core", "extended"])].reset_index(drop=True)
    models = [model_store.LaplacianModel.load(os.path.join(GENE_MOL, f))
              for f in man.model_file]
    return man, models

def _none_if_nan(v):
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else v

LIABILITY = {
    "KCNH2": "hERG potassium channel — QT prolongation, arrhythmia",
    "CYP3A4": "major drug-metabolising enzyme — interaction risk",
    "CYP2D6": "drug-metabolising enzyme — interaction risk",
    "CYP2C9": "drug-metabolising enzyme — interaction risk",
    "HTR2B": "5-HT2B — valvular heart disease with chronic agonism",
    "AR": "androgen receptor — endocrine effects",
    "NR3C1": "glucocorticoid receptor — endocrine and growth effects",
    "PGR": "progesterone receptor — endocrine effects",
    "ADRB1": "beta-1 adrenoceptor — cardiac effects",
    "SCN5A": "cardiac sodium channel — conduction effects",
}

@functools.lru_cache(maxsize=1)
def _safety_tables():
    base = os.path.join(KG, "drug_safety")
    return (pd.read_parquet(f"{base}.parquet"),
            pd.read_parquet(f"{base}_warnings.parquet"),
            pd.read_parquet(f"{base}_indications.parquet"),
            pd.read_parquet(f"{base}_atc.parquet"))

def drug_safety(drug_name: str) -> dict:
    """Withdrawals, black-box warnings, route and what the drug was trialled for.

    Development history is evidence: a compound trialled repeatedly and never
    approved usually failed for a reason, and the indications it was tried in
    say what exposure pattern it was designed for.
    """
    mol, warn, ind, atc = _safety_tables()
    q = str(drug_name).strip().upper()
    hit = mol[mol.pref_name.fillna("").str.upper() == q]
    if hit.empty:
        return {"ok": False, "error": "drug not found in ChEMBL by name",
                "drug": drug_name}
    r = hit.sort_values("max_phase", ascending=False).iloc[0]
    cid = r.chembl_id
    w = warn[warn.chembl_id == cid]
    i = ind[ind.chembl_id == cid].sort_values("max_phase_for_ind", ascending=False)
    a = atc[atc.chembl_id == cid]
    routes = [k for k, v in (("oral", r.oral), ("parenteral", r.parenteral),
                             ("topical", r.topical)) if v == 1]
    return {
        "ok": True, "drug": r.pref_name, "chembl_id": cid,
        "max_phase": _none_if_nan(r.max_phase),
        "first_approval_year": _none_if_nan(r.first_approval),
        "withdrawn": bool(r.withdrawn_flag == 1),
        "black_box_warning": bool(r.black_box_warning == 1),
        "routes": routes or None,
        "warnings": [{"type": x.warning_type, "class": x.warning_class,
                      "description": x.warning_description,
                      "country": x.warning_country, "year": _none_if_nan(x.warning_year)}
                     for x in w.itertuples()][:12],
        "atc_class": sorted({x for x in a.level2_description.dropna()}) or None,
        "n_indications_tried": int(len(i)),
        "indications": [{"indication": x.mesh_heading or x.efo_term,
                         "max_phase_reached": float(x.max_phase_for_ind)}
                        for x in i.head(10).itertuples()],
        "note": ("No warning recorded is not evidence of safety: ChEMBL records "
                 "withdrawals and black-box warnings for approved drugs, so an "
                 "investigational compound will usually have none simply because "
                 "it never reached a label."),
        "source": "ChEMBL 37 drug_warning / drug_indication / ATC",
    }

def safety_panel(smiles: str, top: int = 6) -> dict:
    """Predicted off-target liabilities for a structure (hERG, CYPs, and so on).

    Runs the NB panel against a fixed list of targets whose engagement carries
    a known clinical liability. Same caveat as every structure-based call here:
    it recognises chemotype, so a hit is a reason to measure, not a finding.
    """
    man, models = _serving()
    X = smiles_to_csr([smiles], radius=2, n_jobs=1, verbose=False)
    if X.nnz == 0:
        return {"ok": False, "error": "could not parse SMILES", "smiles": smiles}
    out = []
    for gene, why in LIABILITY.items():
        row = man.index[man.gene == gene]
        if len(row) == 0:
            continue
        i = int(row[0])
        z = float(models[i].zscore(X)[0])
        p = (_ext_precision(man.iloc[i].n_train_active, z)
             if pd.notna(man.iloc[i].n_train_active) else None)
        out.append({"gene": gene, "liability": why, "z": round(z, 2),
                    "estimated_precision": None if p is None else round(p, 3)})
    out.sort(key=lambda d: -d["z"])
    flagged = [d for d in out if d["z"] >= 4.0]
    return {"ok": True, "smiles": smiles,
            "n_liability_targets_screened": len(out),
            "flagged": flagged,
            "all_scores": out[:top],
            "verdict": ("flags" if flagged else "no liability target above z=4"),
            "note": ("Structure-based prediction over liability targets. A flag "
                     "warrants a measured assay (hERG patch clamp, CYP "
                     "inhibition panel); absence of a flag is weak evidence of "
                     "safety, not a clearance."),
            "source": "gene_mol NB panel"}

def clean(obj):
    """Replace pandas/numpy NaN with None, recursively.

    NaN is truthy, so it slips past "if not value" guards, and json.dumps emits
    a bare NaN token that is not valid JSON - which would put unparseable tool
    output into the training data.
    """
    if isinstance(obj, dict):
        return {k: clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, float) and (obj != obj):          # NaN
        return None
    if obj is pd.NaT or (obj is not None and obj is pd.NA):
        return None
    try:
        if pd.isna(obj):
            return None
    except (TypeError, ValueError):
        pass
    return obj

@functools.lru_cache(maxsize=8)
def _kg(name: str) -> pd.DataFrame:
    return pd.read_parquet(f"{KG}/{name}.parquet")

def lookup_gene(symbol: str) -> dict:
    """Function, location, constraint and pathways for a gene symbol."""
    g = _kg("genes")
    hit = g[g.symbol.str.upper() == symbol.upper()]
    if hit.empty:
        return {"ok": False, "error": "gene symbol not found", "symbol": symbol,
                "note": "Not in Open Targets; check the symbol or use an alias."}
    r = hit.iloc[0]
    paths = _kg("gene_pathway")
    p = paths[paths.gene_id == r.gene_id]
    return {
        "ok": True, "symbol": r.symbol, "gene_id": r.gene_id, "name": r["name"],
        "uniprot": r.uniprot, "function": r.function,
        "subcellular_location": r.subcellular, "target_class": r.target_class,
        "lof_constraint_score": None if pd.isna(r.pli) else float(r.pli),
        "pathways": p[["pathway_id", "pathway", "top_level"]].head(25).to_dict("records"),
        "source": "Open Targets 26.09 (UniProt, Reactome)",
    }

@functools.lru_cache(maxsize=1)
def _clinvar() -> pd.DataFrame:
    """The complete ClinVar table, VUS and conflicting records included.

    rare_real's clinvar_set.tsv.gz drops Uncertain_significance and
    Conflicting_classifications because it is a classifier training set. Using
    it here made the tool answer "not in ClinVar" for ~1.3M VUS records - a
    false negative, and the most dangerous kind for this use.
    """
    return pd.read_parquet(CLINVAR_FULL)

def lookup_variant(variant: str) -> dict:
    """ClinVar record for a variant given as chrom:pos:ref:alt."""
    cv = _clinvar()
    hit = cv[cv.var_key == variant]
    if hit.empty:
        return {"ok": False, "error": "variant not in ClinVar", "variant": variant,
                "note": ("Absence from ClinVar means it has not been submitted, "
                         "not that it is benign. Judging it needs the gene, the "
                         "predicted consequence and population frequency.")}
    r = hit.iloc[0]
    sig = str(r.clnsig)
    stars = int(r.stars)
    out = {"ok": True, "variant": variant, "gene": r.genes,
           "clinical_significance": sig, "review_stars": stars,
           "review_status": r.get("review_status"),
           "molecular_consequence": r.get("mc"),
           "variant_type": r.get("clnvc", r.get("var_type")),
           "clinvar_id": r.clinvar_id,
           "reported_conditions": r.get("conditions"),
           "source": "ClinVar 2026-10-04"}
    # Say plainly what the classification does and does not support, so the
    # model is not left to infer confidence from the label alone.
    if "Conflicting" in sig:
        out["interpretation_note"] = (
            "Submitters disagree about this variant. There is no consensus "
            "classification; the disagreement itself is the finding.")
    elif "Uncertain" in sig or sig == "VUS-mid":
        out["interpretation_note"] = (
            "A VUS is not evidence of pathogenicity and not evidence against "
            "it. It should not be used to support or exclude a diagnosis.")
    elif stars <= 1:
        out["interpretation_note"] = (
            f"Only {stars} review star(s): a single submitter or no assertion "
            "criteria, so this classification is weakly supported.")
    return out

def gene_diseases(symbol: str, limit: int = 15) -> dict:
    """Diseases associated with a gene, with Open Targets association scores."""
    g = _kg("genes")
    hit = g[g.symbol.str.upper() == symbol.upper()]
    if hit.empty:
        return {"ok": False, "error": "gene symbol not found", "symbol": symbol}
    gid = hit.iloc[0].gene_id
    gd = _kg("gene_disease")
    d = _kg("diseases").set_index("disease_id")
    rows = gd[gd.gene_id == gid].sort_values("assoc_score", ascending=False).head(limit)
    out = []
    for _, r in rows.iterrows():
        name = d.loc[r.disease_id]["name"] if r.disease_id in d.index else None
        out.append({"disease_id": r.disease_id, "disease": name,
                    "association_score": round(float(r.assoc_score), 3),
                    "n_evidence": int(r.n_evidence) if "n_evidence" in r else None})
    return {"ok": True, "gene": symbol, "diseases": out,
            "note": ("Association score reflects weight of evidence, not causality "
                     "or effect size."),
            "source": "Open Targets 26.09"}

def drugs_for_target(symbol: str, approved_only: bool = True, limit: int = 20) -> dict:
    """Drugs acting on a gene's protein product, with mechanism and phase."""
    g = _kg("genes")
    hit = g[g.symbol.str.upper() == symbol.upper()]
    if hit.empty:
        return {"ok": False, "error": "gene symbol not found", "symbol": symbol}
    gid = hit.iloc[0].gene_id
    dt = _kg("drug_target")
    dr = _kg("drugs").set_index("drug_id")
    rows = dt[dt.gene_id == gid]
    out = []
    for _, r in rows.iterrows():
        if r.drug_id not in dr.index:
            continue
        d = dr.loc[r.drug_id]
        if approved_only and not bool(d.is_approved):
            continue
        out.append({"drug_id": r.drug_id, "drug": d.drug_name,
                    "mechanism": r.mechanism, "action_type": r.action_type,
                    "max_clinical_stage": d.max_stage,
                    "approved": bool(d.is_approved),
                    "smiles": d.smiles})
    return {"ok": True, "gene": symbol, "n_found": len(out), "drugs": out[:limit],
            "source": "Open Targets 26.09 / ChEMBL"}

def _wrap(fn):
    """Every tool result goes through clean() before it reaches the model."""
    @functools.wraps(fn)
    def inner(*a, **k):
        return clean(fn(*a, **k))
    return inner

TOOL_SCHEMAS = [{'type': 'function',
  'function': {'name': 'lookup_gene',
               'description': 'Function, subcellular location, loss-of-function constraint and '
                              'Reactome pathways for a human gene.',
               'parameters': {'type': 'object',
                              'properties': {'symbol': {'type': 'string',
                                                        'description': 'HGNC gene symbol, e.g. '
                                                                       'BUB1B'}},
                              'required': ['symbol']}}},
 {'type': 'function',
  'function': {'name': 'lookup_variant',
               'description': 'ClinVar record for one variant given as chrom:pos:ref:alt on '
                              'GRCh38, with clinical significance and review-star confidence.',
               'parameters': {'type': 'object',
                              'properties': {'variant': {'type': 'string',
                                                         'description': 'e.g. 15:40170081:C:T'}},
                              'required': ['variant']}}},
 {'type': 'function',
  'function': {'name': 'gene_diseases',
               'description': 'Diseases associated with a gene and the strength of the evidence, '
                              'from Open Targets.',
               'parameters': {'type': 'object',
                              'properties': {'symbol': {'type': 'string'},
                                             'limit': {'type': 'integer'}},
                              'required': ['symbol']}}},
 {'type': 'function',
  'function': {'name': 'drugs_for_target',
               'description': "Drugs acting on a gene's protein product, with mechanism of action, "
                              'action type and clinical phase. Use to find repurposing candidates '
                              'for a target.',
               'parameters': {'type': 'object',
                              'properties': {'symbol': {'type': 'string'},
                                             'approved_only': {'type': 'boolean'},
                                             'limit': {'type': 'integer'}},
                              'required': ['symbol']}}},
 {'type': 'function',
  'function': {'name': 'drug_safety',
               'description': 'Withdrawals, black-box warnings, administration route and the full '
                              'list of indications a drug has been trialled for with the phase it '
                              'reached. Use on any drug being proposed, to check what its '
                              'development history and label already say about it.',
               'parameters': {'type': 'object',
                              'properties': {'drug_name': {'type': 'string'}},
                              'required': ['drug_name']}}},
 {'type': 'function',
  'function': {'name': 'safety_panel',
               'description': 'Screen a structure against targets whose engagement carries a known '
                              'clinical liability (hERG/QT, CYP interactions, 5-HT2B, endocrine '
                              "receptors). Use on a candidate's SMILES as a toxicity pre-check.",
               'parameters': {'type': 'object',
                              'properties': {'smiles': {'type': 'string'},
                                             'top': {'type': 'integer'}},
                              'required': ['smiles']}}}]

TOOLS = {name: _wrap(globals()[name]) for name in ['lookup_gene', 'lookup_variant', 'gene_diseases', 'drugs_for_target', 'drug_safety', 'safety_panel']}

from tool_codependency import codependency as _fn, SCHEMA as _schema
TOOLS["codependency"] = _wrap(_fn)
TOOL_SCHEMAS.append(_schema)

from tool_sibling_diseases import sibling_diseases as _fn, SCHEMA as _schema
TOOLS["sibling_diseases"] = _wrap(_fn)
TOOL_SCHEMAS.append(_schema)

from tool_essentiality import target_essentiality as _fn, SCHEMA as _schema
TOOLS["target_essentiality"] = _wrap(_fn)
TOOL_SCHEMAS.append(_schema)

from tool_literature import literature_evidence as _fn, SCHEMA as _schema
TOOLS["literature_evidence"] = _wrap(_fn)
TOOL_SCHEMAS.append(_schema)

from tool_protein_regulation import protein_regulation as _fn, SCHEMA as _schema
TOOLS["protein_regulation"] = _wrap(_fn)
TOOL_SCHEMAS.append(_schema)

# Preserve the exact final-run schema order used in model prompts.
_SCHEMA_ORDER = ["lookup_gene", "lookup_variant", "gene_diseases", "drugs_for_target", "target_essentiality", "drug_safety", "safety_panel", "codependency", "sibling_diseases", "literature_evidence", "protein_regulation"]
TOOL_SCHEMAS.sort(key=lambda item: _SCHEMA_ORDER.index(item["function"]["name"]))
