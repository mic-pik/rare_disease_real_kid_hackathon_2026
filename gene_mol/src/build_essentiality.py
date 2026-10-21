"""Per-gene CRISPR essentiality, for the question "can this be inhibited chronically?".

DepMap Chronos gene effect: 0 = knockout does nothing, -1 = the median common
essential gene. `prism_tool/crispr_z.npy` is z-scored per gene, which is right
for correlating a drug's killing profile against a dependency profile but
destroys the absolute signal - a ribosomal protein looks like any other gene.
This keeps the absolute scores.

The distinction that matters for a therapeutic hypothesis:
  selective  a dependency in some lineages only (EGFR, BRAF) - a drug target
  pan-essential  required by nearly every dividing cell (AURKB, PLK1, RPL5) -
                 inhibiting it systemically is antiproliferative chemotherapy,
                 whatever the genotype-specific rationale says

Why the spread, not just the mean (added 2026-10-15)
----------------------------------------------------
A mean alone structurally cannot see a dependency that is selective in the
context that matters, and it mislabels the clearest example we have. KIF18A is
the chromosomal-instability-selective target - CIN/aneuploid cells need it to
cluster supernumerary centrosomes, near-diploid cells do not, and it has oral
phase-1 inhibitors (sovilnesib, VLS-1488). Its mean is -1.30, unremarkable
next to AURKB at -2.30, so a mean-only table files it under "broadly required"
and the hypothesis dies there. Its sd is 0.65 - twice BUB1B's 0.31 - and its
range runs -3.84 to +0.14. The variance is the whole signal.

So this build also emits, per gene: sd, skewness, min/max, the 10th/50th/90th
percentiles, the resistant fraction, Sarle's bimodality coefficient, and a
`selectivity_score`.

Choice of selectivity statistic
-------------------------------
`selectivity_score = 2 * sqrt(f * (1 - f)) * depth`
  f     = fraction of lines with Chronos < -0.5 (the same threshold the
          existing frac_dependent column uses, so the two agree by
          construction)
  depth = |mean Chronos of those dependent lines|

Read it as the Chronos killing depth that is actually available *selectively*:
the Bernoulli-spread term 2*sqrt(f(1-f)) is 1 when the panel splits evenly
into dependent and resistant lines and goes to 0 as every line (pan-essential)
or no line (irrelevant gene) depends on it, and it is multiplied by how hard
the dependent lines actually die. A gene scores high only if it both splits
the panel and kills the half it splits off.

Rejected alternatives:
  * sd alone - ranks KIF18A (0.65) over BUB1B (0.31) correctly but calls EGFR
    (0.335) no more selective than AURKB (0.336). It cannot tell a wide
    distribution from one with a genuinely resistant population.
  * Sarle's bimodality coefficient - computed and kept as a descriptive
    column, but it is driven by skew, so it scores strongly left-skewed
    rarely-essential genes (WRN 0.60, BRAF 0.58) above KIF18A (0.20) and only
    barely separates KIF18A from the pan-essentials (AURKB 0.14, RPL5 0.16).
  * skewed-t likelihood ratio (DepMap's normLRT) - also skew-driven, and
    KIF18A and BUB1B have essentially the same skew (-0.35 vs -0.34), so it
    cannot make the one separation this fix exists to make.
  * Hartigan's dip - no dependable implementation in this venv, and its
    p-value needs tabulated critical values.
The chosen measure is threshold-based and therefore sensitive to the -0.5
cutoff; that cutoff is inherited from the existing columns rather than tuned,
and the raw percentiles are emitted alongside so the number can be checked.

Every pre-existing column and every essentiality_class value is unchanged;
the new columns are additive.

Usage:
    python -I src/build_essentiality.py --out work/essentiality.parquet
"""
import argparse
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _bimodality_coefficient(v: np.ndarray) -> float:
    """Sarle's BC = (skew^2 + 1) / (kurtosis + 3(n-1)^2/((n-2)(n-3))).

    >0.555 (the value for a uniform distribution) is the conventional hint of
    bimodality. Descriptive only here - see the module docstring for why it is
    not the primary statistic.
    """
    n = len(v)
    if n < 4:
        return float("nan")
    sk = float(stats.skew(v))
    ku = float(stats.kurtosis(v))  # excess
    denom = (ku + 3.0) + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    return (sk ** 2 + 1.0) / denom if denom else float("nan")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--crispr", default="data/depmap_raw/crispr_gene_effect.csv")
    p.add_argument("--out", default="work/essentiality.parquet")
    args = p.parse_args()
    os.chdir(ROOT)

    d = pd.read_csv(args.crispr, index_col=0)
    print(f"{d.shape[0]:,} cell lines x {d.shape[1]:,} genes", flush=True)
    rows = []
    for col in d.columns:
        gene = re.sub(r"\s*\(\d+\)$", "", col).strip()
        v = d[col].to_numpy(dtype=float)
        v = v[~np.isnan(v)]
        if len(v) < 100:
            continue
        mean = float(v.mean())
        frac = float((v < -0.5).mean())
        strong = float((v < -1.0).mean())
        if frac >= 0.90:
            cls = "pan-essential"
        elif frac >= 0.25:
            cls = "broadly required"
        elif frac >= 0.03:
            cls = "selectively essential"
        else:
            cls = "rarely essential"

        dep = v[v < -0.5]
        depth = float(-dep.mean()) if len(dep) else 0.0
        sel = 2.0 * float(np.sqrt(frac * (1.0 - frac))) * depth
        q10, q50, q90 = (float(x) for x in np.percentile(v, [10, 50, 90]))

        rows.append({"gene": gene, "mean_chronos": round(mean, 3),
                     "frac_dependent": round(frac, 4),
                     "frac_strongly_dependent": round(strong, 4),
                     "n_lines": len(v), "essentiality_class": cls,
                     # --- added 2026-10-15: the spread the mean was hiding ---
                     "sd_chronos": round(float(v.std(ddof=1)), 3),
                     "min_chronos": round(float(v.min()), 3),
                     "max_chronos": round(float(v.max()), 3),
                     "p10_chronos": round(q10, 3),
                     "median_chronos": round(q50, 3),
                     "p90_chronos": round(q90, 3),
                     "skew_chronos": round(float(stats.skew(v)), 3),
                     "frac_resistant": round(float((v > -0.2).mean()), 4),
                     "mean_chronos_in_dependent_lines": round(-depth, 3),
                     "bimodality_coef": round(_bimodality_coefficient(v), 4),
                     "selectivity_score": round(sel, 4)})
    out = pd.DataFrame(rows).drop_duplicates("gene")
    # Percentile rank of the selectivity score across the genome, so a tool can
    # say "more selective than X% of genes" without re-deriving a threshold.
    out["selectivity_pct"] = (out.selectivity_score.rank(pct=True) * 100).round(1)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_parquet(args.out, index=False)
    print(out.essentiality_class.value_counts().to_string())

    print("\nsanity check (selectivity_score, higher = more selective):")
    chk = out[out.gene.isin(["KIF18A", "BUB1B", "EGFR", "BRAF", "AURKB",
                             "PLK1", "RPL5", "CDK1", "MAD2L1", "TP53"])]
    print(chk.sort_values("selectivity_score", ascending=False)[
        ["gene", "mean_chronos", "sd_chronos", "frac_dependent",
         "bimodality_coef", "selectivity_score", "selectivity_pct",
         "essentiality_class"]].to_string(index=False))
    print(f"\nwrote {args.out} ({len(out):,} genes)")


if __name__ == "__main__":
    main()
