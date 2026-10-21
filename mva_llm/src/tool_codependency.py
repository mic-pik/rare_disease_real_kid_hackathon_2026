"""Gene-gene co-dependency from DepMap CRISPR.

Built by ../gene_mol/src/build_codependency.py; see that docstring for method
and storage choices. This module is standalone so it can be folded into
tools.py without touching anything already there.
"""

from __future__ import annotations

import functools
import os

import pandas as pd

KG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                  "data/kg")
CODEP = os.path.join(KG, "codependency.parquet")

# Anchors for reading the magnitude. These are real values from this table.
_REFERENCE = {
    "CTNNB1 - TCF7L2 (obligate transcriptional partners)": 0.68,
    "BRAF - MAP2K1 (direct kinase-substrate)": 0.47,
    "EGFR - SHC1 (receptor and its adaptor)": 0.40,
    "EIF4E - EIF4G1 (same complex)": 0.22,
    "median gene's single best partner": 0.26,
}

_INTERPRETATION = (
    "These correlations are modest in absolute terms - a strong co-dependency "
    "here is around 0.25, and the whole-genome median for a gene's single best "
    "partner is 0.26. Do NOT dismiss a result for having r=0.25: that is the "
    "normal scale of this measurement, not a weak version of a bigger effect. "
    "What validates a result is the *composition* of the top list - whether the "
    "partners form a coherent complex or pathway - not the magnitude of any one "
    "correlation. Equally, do not over-read it: r=0.25 means these genes are "
    "needed under similar conditions across cell lines, which is a hypothesis "
    "about shared requirement. It is not evidence of a physical interaction, "
    "not a direction of causality, and not a synthetic-lethal relationship, "
    "which would need the dependency measured specifically in cells carrying "
    "the lesion. Partners are also not screened for druggability or for "
    "selectivity: check target_essentiality on anything before proposing it, "
    "because a co-dependency can perfectly well be a pan-essential gene, in "
    "which case inhibiting it is chemotherapy."
)


@functools.lru_cache(maxsize=1)
def _codep() -> pd.DataFrame:
    return pd.read_parquet(CODEP)


@functools.lru_cache(maxsize=1)
def _genes() -> frozenset:
    return frozenset(_codep().gene.unique())


def codependency(gene: str, top: int = 25) -> dict:
    """Genes whose CRISPR dependency profile tracks this gene's across the panel.

    Two genes that cell lines need to the same degree, line by line, are
    needed under the same conditions - same complex, same pathway, or same
    buffering relationship - whether or not anyone has annotated them
    together. This is the data-driven alternative to hand-curating a broken
    gene's pathway neighbours, which tends to return the textbook core of the
    pathway and therefore the genes every dividing cell needs.
    """
    g = str(gene).upper().strip()
    if g not in _genes():
        return {"ok": False, "gene": gene,
                "error": "gene not in the co-dependency table",
                "note": ("Either it is absent from the DepMap CRISPR library or "
                         "it was missing in >10% of cell lines and dropped, in "
                         "which case its profile is too sparse to correlate. "
                         "This is not evidence the gene has no partners.")}
    c = _codep()
    hits = c[c.gene == g].nsmallest(max(1, int(top)), "rank")
    return {
        "ok": True, "gene": g,
        "n_partners_returned": int(len(hits)),
        "n_cell_lines": 1178,
        "partners": [{"gene": h.partner,
                      "correlation": round(float(h.correlation), 3),
                      "rank": int(h.rank)}
                     for h in hits.itertuples()],
        "interpretation": _INTERPRETATION,
        "reference_points": _REFERENCE,
        "method": ("Pearson correlation of Chronos gene-effect profiles across "
                   "1,178 cell lines; genes missing in >10% of lines dropped, "
                   "remaining gaps filled with the gene's panel mean. Only the "
                   "top 50 positive partners per gene are stored."),
        "source": "DepMap CRISPR (Chronos) gene effect, gene-gene correlation",
    }


SCHEMA = {"type": "function", "function": {
    "name": "codependency",
    "description": (
        "Genes whose CRISPR dependency profile correlates with this gene's "
        "across 1,178 DepMap cell lines - i.e. genes needed under the same "
        "conditions, which is the data-driven way to go from a broken or "
        "lost gene to candidate vulnerabilities. Use it instead of listing a "
        "gene's known pathway neighbours from memory, which returns the "
        "textbook pathway core and therefore mostly pan-essential genes. "
        "Correlations are modest by construction (~0.25 is strong here); judge "
        "a result by whether the top partners cohere, not by magnitude. "
        "Complements dependency_targets, which starts from a drug rather than "
        "a gene."),
    "parameters": {"type": "object", "properties": {
        "gene": {"type": "string", "description": "HGNC symbol, e.g. BUB1B"},
        "top": {"type": "integer",
                "description": "partners to return, default 25, max 50"}},
        "required": ["gene"]}}}
