"""target_essentiality v2 - the mean plus the spread it was hiding.

v1 reported a single mean Chronos score per gene and a class derived from the
fraction of dependent lines. That is structurally blind to a dependency that
is selective in the context that matters: KIF18A, the chromosomal-instability
-selective target with oral phase-1 inhibitors (sovilnesib, VLS-1488), has a
mean of -1.30 and lands in "broadly required" next to genes nothing can be
done about, because its mean is unremarkable while its sd (0.65, twice
BUB1B's 0.31) and its range (-3.84 to +0.14) are not.

This module keeps every v1 field and class value and adds the distributional
ones, plus interpretation text that tells the model a high-variance "broadly
required" gene is not disqualified. Table built by
../gene_mol/src/build_essentiality.py, which justifies the statistic.
"""

from __future__ import annotations

from pathlib import Path
import functools
import os

import pandas as pd

GENE_MOL = str(Path(__file__).resolve().parents[2] / "gene_mol")
ESSENTIALITY = os.path.join(GENE_MOL, "work/essentiality.parquet")

# Whole-genome percentiles of selectivity_score, for turning the scalar into a
# verdict. Score 0 means no line depends on the gene at all OR every line does.
_SELECTIVE = 0.45      # ~92nd percentile
_SOME = 0.20           # ~82nd percentile

_CLASS_NOTE = {
    "pan-essential": ("Required by nearly every cell line tested. Systemic "
                      "inhibition is antiproliferative chemotherapy rather "
                      "than a targeted intervention; expect dose-limiting "
                      "marrow, gut and growth toxicity, and do not assume a "
                      "genotype-specific window exists."),
    "broadly required": ("Needed by a large fraction of lines. A window may "
                         "exist but has to be demonstrated, not assumed."),
    "selectively essential": ("Needed by a minority of lines, which is the "
                              "profile of a tractable target - the question "
                              "becomes whether this patient's cells are among "
                              "them."),
    "rarely essential": ("Few lines depend on it. Loss is generally "
                         "tolerated, so efficacy rather than toxicity is "
                         "the harder question."),
}

_SPREAD_NOTE = (
    "Read the spread, not only the mean. The mean answers 'how much does the "
    "average cell need this', which is the wrong question when the point is "
    "to kill one kind of cell and spare another. A gene can be classed "
    "'broadly required' on its mean and still be a real target if its "
    "distribution is wide and a substantial population of lines tolerates its "
    "loss: that spread is where a therapeutic window lives. KIF18A is the "
    "worked example - mean -1.30, sd 0.65, range -3.84 to +0.14, classed "
    "'broadly required', and nonetheless the chromosomal-instability-selective "
    "target now in oral phase-1 trials. The converse also holds: a deep mean "
    "with a narrow spread (AURKB -2.30, sd 0.34, every line dependent) is "
    "pan-essential and no amount of genotype rationale creates a window. "
    "So: do not reject a high-variance gene because its mean looks pan-"
    "essential-ish, and do not accept a low-variance deep one because the "
    "biology story is good. What the spread does NOT tell you is which "
    "lines sit in the resistant mode or why - this table has no lineage or "
    "genotype breakdown, so 'selective' here is a hypothesis that the "
    "patient's cells are on the sensitive side, to be checked against the "
    "specific lesion, not a demonstration of it."
)

_REFERENCE_POINTS = {
    "RPL5 (ribosomal, pan-essential)": {"mean": -2.55, "sd": 0.35,
                                        "selectivity_score": 0.0},
    "AURKB (pan-essential)": {"mean": -2.30, "sd": 0.34,
                              "selectivity_score": 0.0},
    "PLK1 (pan-essential)": {"mean": -2.69, "sd": 0.46,
                             "selectivity_score": 0.16},
    "BUB1B (deep mean, narrow spread)": {"mean": -1.05, "sd": 0.31,
                                         "selectivity_score": 0.35},
    "BRAF (selective, lineage/genotype)": {"mean": -0.11, "sd": 0.33,
                                           "selectivity_score": 0.57},
    "EGFR (selective, 19% of lines)": {"mean": -0.29, "sd": 0.34,
                                       "selectivity_score": 0.67},
    "KIF18A (selective despite a deep mean)": {"mean": -1.30, "sd": 0.65,
                                               "selectivity_score": 0.91},
}


@functools.lru_cache(maxsize=1)
def _essentiality() -> pd.DataFrame:
    return pd.read_parquet(ESSENTIALITY).set_index("gene")


def _verdict(sel: float, cls: str, frac_dep: float) -> str:
    # A window needs somewhere to hide: if essentially every line is dependent
    # the score can still be non-trivial (BUB1B, 0.35, because the dependent
    # lines die hard) and that must not read as a window.
    resistant = 1.0 - frac_dep
    if resistant >= 0.05:
        if sel >= _SELECTIVE:
            return "selective dependency - a therapeutic window is plausible"
        if sel >= _SOME:
            return "partially selective - a window is possible but unproven"
    if cls in ("pan-essential", "broadly required"):
        return ("uniform dependency - almost every line is dependent, so there "
                "is no window; inhibition is chemotherapy")
    return "no appreciable dependency in any line - efficacy is the open question"


def target_essentiality(gene: str) -> dict:
    """How much every cell needs this gene, and how much that varies by line.

    DepMap Chronos gene effect across ~1,200 cell lines: 0 means knockout does
    nothing, -1 is the median common essential gene. A target that is
    *selectively* essential can be inhibited with a therapeutic window; one
    that is pan-essential cannot, however good the genotype rationale is, and
    that is the difference between a targeted agent and chemotherapy. v2 adds
    the spread, because selectivity lives in the variance and the mean hides
    it.
    """
    e = _essentiality()
    g = str(gene).upper().strip()
    if g not in e.index:
        return {"ok": False, "error": "gene not in the CRISPR screen", "gene": gene}
    r = e.loc[g]
    cls = str(r.essentiality_class)
    sel = float(r.selectivity_score)
    return {
        "ok": True, "gene": g,
        # --- v1 fields, unchanged ---
        "mean_chronos_gene_effect": float(r.mean_chronos),
        "fraction_of_lines_dependent": float(r.frac_dependent),
        "fraction_strongly_dependent": float(r.frac_strongly_dependent),
        "n_cell_lines": int(r.n_lines),
        "essentiality_class": cls,
        # --- v2: the distribution behind the mean ---
        "sd_chronos_gene_effect": float(r.sd_chronos),
        "chronos_range": [float(r.min_chronos), float(r.max_chronos)],
        "chronos_percentiles": {"p10": float(r.p10_chronos),
                                "median": float(r.median_chronos),
                                "p90": float(r.p90_chronos)},
        "mean_chronos_in_dependent_lines": float(r.mean_chronos_in_dependent_lines),
        "fraction_of_lines_resistant": float(r.frac_resistant),
        "skewness": float(r.skew_chronos),
        "bimodality_coefficient": float(r.bimodality_coef),
        "selectivity_score": round(sel, 3),
        "selectivity_percentile_genome_wide": float(r.selectivity_pct),
        "selectivity_verdict": _verdict(sel, cls, float(r.frac_dependent)),
        "interpretation": _CLASS_NOTE[cls] + " " + _SPREAD_NOTE,
        "selectivity_score_definition": (
            "2*sqrt(f*(1-f)) * |mean Chronos of the dependent lines|, where f "
            "is the fraction of lines below -0.5. It is the killing depth "
            "available *selectively*: maximal when the panel splits into "
            "dependent and resistant halves and the dependent half dies hard, "
            "and 0 when every line depends on the gene (pan-essential) or none "
            "does. A score of 0 is therefore ambiguous on its own - read it "
            "with essentiality_class. Roughly: >=0.45 is selective (~92nd "
            "percentile genome-wide), 0.20-0.45 partially selective - but only if at least 5% of lines are non-dependent, since a window needs a resistant population."),
        "reference_points": _REFERENCE_POINTS,
        "source": "DepMap CRISPR (Chronos) gene effect, 1,178 cell lines",
    }


SCHEMA = {"type": "function", "function": {
    "name": "target_essentiality",
    "description": (
        "How broadly cells depend on a gene, from DepMap CRISPR - the mean "
        "Chronos effect AND the spread across cell lines (sd, range, "
        "percentiles, a selectivity score). Use before proposing that a target "
        "be inhibited in a patient, and especially before proposing chronic "
        "dosing: a pan-essential target cannot be inhibited selectively "
        "whatever the genotype rationale, while a selectively essential one "
        "can. Check the spread and not just the mean - a gene with a deep mean "
        "but a wide, bimodal distribution (KIF18A) can still be a tractable "
        "target in the right context, whereas a deep mean with a narrow spread "
        "(AURKB, PLK1, RPL5) cannot."),
    "parameters": {"type": "object", "properties": {
        "gene": {"type": "string", "description": "HGNC symbol, e.g. KIF18A"}},
        "required": ["gene"]}}}
