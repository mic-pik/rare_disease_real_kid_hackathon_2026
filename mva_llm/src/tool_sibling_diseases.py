"""What is actually given to patients with the *other* genome-instability syndromes.

Why this tool exists
--------------------
For mosaic variegated aneuploidy (biallelic BUB1B) the target axis is
exhausted: everything reachable from the mitotic-checkpoint pathway is
pan-essential, so `predict_targets` -> `target_essentiality` correctly ends at
"this is chemotherapy, not therapy". This tool opens the second axis -
mechanism-class transfer. Fanconi anaemia, ataxia-telangiectasia, Bloom,
Werner, Nijmegen breakage, dyskeratosis congenita, Rothmund-Thomson and their
relatives share the *failure mode* - constitutional genome or chromosome
instability with growth failure, marrow failure and cancer predisposition -
and differ in the molecular lesion. What has been trialled in those children
is a hypothesis here that already carries paediatric safety and dosing data.

This is explicitly not phenotype similarity. HPO overlap was tested on this
patient and returned noise (pemphigoid gestationis, Legg-Calve-Perthes),
because the patient's HPO annotation is generic. The sibling set is a curated,
enumerated list of disorders selected on mechanism; see
`gene_mol/src/build_sibling_diseases.py` for the inclusion rule and provenance.

Signature reasoning
-------------------
The useful entry points are a *gene* and a *disease*, and both resolve to the
same object - a set of sibling disorders - so one function takes either:

  sibling_diseases()                     the whole mechanism class
  sibling_diseases(gene="BUB1B")         resolve the gene to its disorder, then
                                         return the *other* members of its class
  sibling_diseases(disease="Fanconi Anemia")   one named sibling

The return is deliberately *not* a phase-ranked list. Sorting Fanconi's 23
agents by phase puts busulfan, fludarabine, cyclophosphamide and alemtuzumab
on top - transplant conditioning, which treats the consequence and transfers
nothing. The agents worth transferring (metformin, quercetin, niacinamide,
amifostine) sit at phase 1-2 and would be buried. So the result is bucketed by
a hand-curated `intent` first and ranked by phase only within a bucket.
"""

from __future__ import annotations

import functools
import os

import pandas as pd

KG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                  "data", "kg")
PARQUET = os.path.join(KG, "sibling_diseases.parquet")

MECHANISM_CLASS = "constitutional genome / chromosome instability"

_BUCKETS = [
    ("instability_directed", "agents_aimed_at_the_instability_itself"),
    ("other_mechanism_directed", "agents_aimed_at_a_different_molecular_lesion"),
    ("consequence_directed", "agents_treating_the_consequence"),
    ("symptomatic", "symptomatic_agents"),
    ("gene_specific_correction", "gene_specific_therapies"),
    ("unclassified", "unclassified_agents"),
]

_CAVEAT = (
    "An agent trialled in a sibling syndrome is a HYPOTHESIS TRANSFER, not "
    "evidence in this disease. The only thing shared is the mechanism class "
    "(constitutional genome instability); the molecular lesion differs - "
    "crosslink repair in Fanconi, ATM signalling in ataxia-telangiectasia, "
    "spindle-checkpoint failure in mosaic variegated aneuploidy - so a drug "
    "may transfer, may be irrelevant, or may be actively harmful. Published "
    "exposure history can inform follow-up, but age, dose, duration, endpoint "
    "and adverse events must be checked in the primary study. Tolerability "
    "does not transfer automatically between diseases or populations."
)

_SPLIT = (
    "Read the buckets, not the phase column. 'agents_treating_the_consequence' "
    "are transplant conditioning (busulfan, treosulfan, fludarabine, "
    "melphalan, cyclophosphamide, alemtuzumab), GvHD prophylaxis (methotrexate, "
    "ciclosporin, mycophenolate, steroids) and marrow support (eltrombopag, "
    "filgrastim, danazol, oxandrolone). They say a transplant pathway exists "
    "and nothing about modifying instability, and several are themselves "
    "genotoxic. 'agents_aimed_at_the_instability_itself' - metformin, "
    "quercetin, niacinamide, amifostine - are the transferable hypotheses, and "
    "they sit at phase 1-2, below the conditioning agents, so any phase-ranked "
    "list buries them. 'agents_aimed_at_a_different_molecular_lesion' "
    "(lonafarnib, progerinin) are disease-modifying but for a lesion this "
    "patient does not have. 'gene_specific_therapies' (mozafancogene "
    "autotemcel) are the strongest results in the table and the least "
    "transferable."
)


@functools.lru_cache(maxsize=1)
def _table() -> pd.DataFrame:
    return pd.read_parquet(PARQUET)


def _none_if_nan(v):
    return None if v is None or pd.isna(v) else v


def _agent(r) -> dict:
    return {
        "drug": r.drug,
        "trialled_in": r.disease,
        "max_phase_for_indication": _none_if_nan(r.max_phase_for_ind),
        "molecule_approved_somewhere": bool(_none_if_nan(r.molecule_max_phase) == 4),
        "mechanism_of_action": _none_if_nan(r.moa),
        "atc_class": _none_if_nan(r.atc_class),
        "black_box_warning": bool(r.black_box_warning),
        "withdrawn": bool(r.withdrawn),
        "why_it_was_given": r.intent_rationale,
        "chembl_id": r.chembl_id,
    }


def _resolve_disease(df: pd.DataFrame, disease: str):
    q = str(disease).strip().lower()
    names = df.disease.str.lower()
    mesh = df.mesh_heading.fillna("").str.lower()
    hit = df[(names == q) | (mesh == q)]
    if hit.empty:
        hit = df[names.str.contains(q, regex=False)
                 | mesh.str.contains(q, regex=False)]
    return hit


def sibling_diseases(gene: str | None = None, disease: str | None = None,
                     intent: str | None = None, min_phase: float = 0.0,
                     include_adjacent: bool = False, top: int = 40) -> dict:
    """Agents trialled in the other constitutional genome-instability syndromes.

    Mechanism-class transfer, not phenotype similarity: the sibling set is a
    curated list of disorders whose primary defect is genome or chromosome
    maintenance (DNA repair, replication-fork and helicase function, damage
    signalling, telomere maintenance, sister-chromatid cohesion, mitotic
    segregation). Results are bucketed by what the agent was actually *for* -
    the instability itself versus its downstream consequences - because a
    phase-ranked list is dominated by transplant conditioning.

    Args:
        gene: HGNC symbol. Resolved to its disorder; the siblings returned are
            the other members of the class. Use for the index patient's gene.
        disease: disorder name or ChEMBL MeSH heading, e.g. "Fanconi Anemia".
        intent: restrict to one bucket, e.g. "instability_directed".
        min_phase: drop indications below this max clinical phase.
        include_adjacent: also include the adjacent tier (progeria,
            Shwachman-Diamond, Diamond-Blackfan) - shared clinical problem,
            different mechanism class.
        top: cap on agents returned per bucket.
    """
    df = _table()
    tiers = ["core", "adjacent"] if include_adjacent else ["core"]
    pool = df[df.tier.isin(tiers)]

    query = {"gene": gene, "disease": disease, "intent": intent,
             "min_phase": min_phase, "include_adjacent": include_adjacent}
    notes = []
    home = None

    if disease is not None:
        hit = _resolve_disease(df, disease)
        if hit.empty:
            return {"ok": False,
                    "error": ("disease not in the curated genome-instability "
                              "sibling set"),
                    "query": query,
                    "known_disorders": sorted(df.disease.unique().tolist()),
                    "source": "curated sibling set, see build_sibling_diseases.py"}
        pool = hit
        if not hit.tier.isin(tiers).any():
            notes.append(f"'{hit.disease.iloc[0]}' is tier "
                         f"'{hit.tier.iloc[0]}', returned because it was asked "
                         "for by name.")

    if gene is not None:
        g = str(gene).strip().upper()
        mask = df.genes.fillna("").str.split().apply(lambda xs: g in xs)
        matched = df[mask]
        if matched.empty:
            notes.append(f"{g} is not a causal gene of any disorder in the "
                         "curated set, so the whole mechanism class is "
                         "returned rather than a gene-specific view.")
        else:
            home = matched.disease.iloc[0]
            pool = pool[pool.disease != home]
            notes.append(f"{g} causes '{home}' "
                         f"({matched.mechanism_subclass.iloc[0]}); that "
                         "disorder is excluded and its siblings returned.")

    searched = pool.drop_duplicates("disease")
    agents = pool[pool.drug.notna()].copy()
    if intent is not None:
        agents = agents[agents.intent == intent]
    if min_phase:
        agents = agents[agents.max_phase_for_ind.fillna(-1) >= float(min_phase)]

    buckets = {}
    for key, label in _BUCKETS:
        sub = agents[agents.intent == key].sort_values(
            ["max_phase_for_ind", "drug"], ascending=[False, True])
        if len(sub):
            buckets[label] = [_agent(r) for r in sub.head(top).itertuples()]

    empty = searched[searched.n_agents == 0]
    n_total = int(len(agents))

    if n_total == 0:
        if disease is not None:
            name = searched.disease.iloc[0]
            inchembl = bool(searched.in_chembl.iloc[0])
            err = (f"No agents recorded for '{name}'. "
                   + ("ChEMBL has the MeSH heading but no drug indication "
                      "mapped to it."
                      if inchembl else
                      "ChEMBL has no MeSH indication heading for this disorder "
                      "at all, so nothing could be found even in principle."))
        else:
            err = "No agents matched the filters applied."
        return {
            "ok": True, "query": query, "mechanism_class": MECHANISM_CLASS,
            "n_agents": 0,
            "disorders_searched": searched.disease.tolist(),
            "finding": err,
            "notes": notes or None,
            "interpretation": (
                err + " This is absence of recorded trials, not evidence that "
                "nothing works - these disorders have a few hundred patients "
                "worldwide and most have never had a registered "
                "interventional trial. Do not report it as a negative result. "
                "Query the class as a whole (no disease argument) to see what "
                "has been tried in the better-studied siblings, chiefly "
                "Fanconi anaemia and ataxia-telangiectasia. " + _CAVEAT),
            "source": ("ChEMBL 37 drug_indication (MeSH-coded) over a curated "
                       "genome-instability sibling set; Orphanet gene "
                       "annotation; intent labels hand-curated"),
        }

    return {
        "ok": True,
        "query": query,
        "mechanism_class": MECHANISM_CLASS,
        "index_disorder_excluded": home,
        "n_sibling_disorders_searched": int(len(searched)),
        "disorders_with_agents": [
            {"disease": r.disease, "mechanism_subclass": r.mechanism_subclass,
             "n_agents": int(r.n_agents)}
            for r in searched[searched.n_agents > 0]
            .sort_values("n_agents", ascending=False).itertuples()],
        "disorders_with_no_agents_in_chembl": [
            {"disease": r.disease, "mechanism_subclass": r.mechanism_subclass,
             "reason": ("no drug indication mapped to its MeSH heading"
                        if r.in_chembl else
                        "no MeSH indication heading exists in ChEMBL")}
            for r in empty.itertuples()],
        "n_agents": n_total,
        "n_instability_directed": int((agents.intent == "instability_directed").sum()),
        **buckets,
        "notes": notes or None,
        "interpretation": _CAVEAT + " " + _SPLIT + (
            " Empty buckets mean nothing of that kind was trialled, and the "
            "disorders in 'disorders_with_no_agents_in_chembl' have no "
            "recorded agents at all - that is a gap in the literature for "
            "orphan diseases, not a negative finding. Any candidate taken from "
            "here still needs a stated mechanistic reason why it should act on "
            "*this* lesion, and a check that it is not pan-essential or "
            "genotoxic (`target_essentiality`, `drug_safety`)."),
        "source": ("ChEMBL 37 drug_indication (MeSH-coded) over a curated "
                   "genome-instability sibling set; Orphanet gene annotation; "
                   "intent labels hand-curated, see "
                   "gene_mol/src/build_sibling_diseases.py"),
    }


SCHEMA = {"type": "function", "function": {
    "name": "sibling_diseases",
    "description": (
        "What has actually been given to patients with the OTHER constitutional "
        "genome-instability syndromes (Fanconi anaemia, ataxia-telangiectasia, "
        "Bloom, Werner, Nijmegen breakage, dyskeratosis congenita, "
        "Rothmund-Thomson, xeroderma pigmentosum, cohesinopathies and "
        "relatives). Use when the target-based search for a chromosome-"
        "instability disorder has run out - it transfers hypotheses by "
        "mechanism class, not by phenotype. Agents are bucketed by whether "
        "they were aimed at the instability itself or at its downstream "
        "consequences (transplant conditioning, GvHD prophylaxis, marrow "
        "support); that split, not the clinical phase, is the point. Returns "
        "zero agents honestly for disorders nothing has been trialled in."),
    "parameters": {"type": "object", "properties": {
        "gene": {"type": "string",
                 "description": ("HGNC symbol, e.g. BUB1B. Resolves to its "
                                 "disorder and returns the siblings of that "
                                 "disorder, excluding it.")},
        "disease": {"type": "string",
                    "description": ("A single sibling disorder, by name or "
                                    "ChEMBL MeSH heading, e.g. "
                                    "'Fanconi Anemia'.")},
        "intent": {"type": "string",
                   "enum": ["instability_directed", "other_mechanism_directed",
                            "consequence_directed", "symptomatic",
                            "gene_specific_correction", "unclassified"],
                   "description": "Restrict to one intent bucket."},
        "min_phase": {"type": "number",
                      "description": ("Drop indications below this maximum "
                                      "clinical phase (0-4).")},
        "include_adjacent": {"type": "boolean",
                             "description": ("Include the adjacent tier - "
                                             "progeria, Shwachman-Diamond, "
                                             "Diamond-Blackfan - which share "
                                             "the clinical problem but not the "
                                             "mechanism class. Default false.")},
        "top": {"type": "integer",
                "description": "Cap on agents returned per bucket (default 40)."},
    }, "required": []}}}
