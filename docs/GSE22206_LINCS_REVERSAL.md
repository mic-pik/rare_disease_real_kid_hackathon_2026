# GSE22206 → LINCS: approved-ingredient expression reversal

Completed 2026-10-19. **544 currently eligible FDA ingredients scored; 29 pass the primary computational screen, of which 20 also pass a post-hoc matched-drug-background check. These are expression-reversal hypotheses, not validated rescue drugs.** No existing therapeutic nomination or evidence-review decision was changed.

## Results

The input is the [single-case GSE22206 signature](GSE22206_ANALYSIS.md), restricted to 74 directly measured LINCS landmark genes: 59 increased and 15 decreased. The comparison uses the existing GSE92742 Level 5 matrix. The donor's cells were not treated with these drugs; LINCS profiles come from other cell models.

| Stage | Count |
|---|---:|
| Exact InChIKey matches before FDA product-eligibility filtering | 781 ingredients |
| Exact matches with at least one eligible FDA product | 607 ingredients / 19,204 profiles |
| Quality-filtered, nonoverlapping profiles | 2,630 |
| Ingredients with qualifying profiles | 544 |
| Primary computational follow-up | 29 |
| Also supported by post-hoc matched-drug background | 20 |

The initial metadata inspection counted 3,199 quality profiles before excluding ingredients without eligible products. The final analysis applies the existing FDA catalogue's product eligibility; it does not treat every historical approval as currently eligible. Full InChIKey matching is deliberately conservative and can miss salts or alternative structural representations. No fuzzy names or connectivity-only matches were used.

The highest primary scores among compounds passing both checks are:

| Ingredient | Reversal score | Random-gene BH q | Cell models | Observed conditions |
|---|---:|---:|---:|---|
| Chlorhexidine | 0.289 | 0.0136 | 2 | MCF7 and PC3; 10 µM, 24 h |
| Nisoldipine | 0.220 | 0.0136 | 2 | MCF7 6 h and PC3 24 h; 10 µM |
| Methimazole | 0.195 | 0.0733 | 2 | MCF7 and PC3; 10 µM, 6 h |
| Nitazoxanide | 0.181 | 0.0733 | 2 | A549 and PC3; 10 µM, 6 h |

These are numerical screening leaders, not a suitability ranking. A high reversal score supplies neither a relevant exposure window nor evidence of preserved division, survival or checkpoint function. All four are supported by only two qualifying profiles in different cell models. No compound-specific exposure or clinical suitability review was performed here.

Several previous hypotheses do not receive support from this screen:

- **Metformin:** primary score 0.0333, q 0.549; fibroblast-weighted and larger-effect queries are negative. Fails the primary screen.
- **Sirolimus:** primary score 0.0844, q 0.204. Fails the primary screen.
- **Amifostine:** primary score 0.1157, q 0.295. Fails the primary screen.
- **Acetylcysteine:** passes the primary screen but has sufficient matched-drug background in only one cell model; its adjusted score there is approximately zero. It does not pass the background check. Missing coverage elsewhere is not evidence of no effect.
- Ceritinib and nicotinamide have no scored row under these identity, eligibility and quality rules; absence is not a negative biological result.

## Scoring and controls

1. Retain `trt_cp` signatures with `is_exemplar = 1`, `distil_cc_q75 >= 0.2`, and `distil_nsample >= 3`. Resolve identities to the eligible FDA ingredient catalogue. Check constituent `distil_id` overlap so shared experiments cannot be counted twice. None of the retained profiles required removal for overlap.
2. Restrict the ranking/null universe to 802 landmark genes represented by eligible GSE22206 probes with unambiguous Entrez mapping. Rank each drug's z-scores within that universe, handling ties by average rank, and rescale ranks to [−1, 1]. Only directly measured landmarks are used.
3. Reversal score = half the mean rank of disease-decreased genes minus half the mean rank of disease-increased genes. Positive means opposite transcriptional direction. Both directions have equal total weight despite the 59/15 imbalance. This is a custom descriptive score, **not CMap tau**.
4. Average profiles within exact ingredient/cell/time/dose/unit conditions, then average conditions within each cell, then cells within each ingredient. Do not pick the best dose or best cell. Save every profile and condition score.
5. Evaluate sensitivity queries: 65 genes supported against relatives; the same 74 genes weighted by fibroblast or lymphoblastoid effect magnitude; and 26 genes whose minimum absolute submitted-scale difference is at least 0.25. Tissue weighting is **not** an independent tissue-specific signature analysis.
6. Draw 10,000 random disjoint 59-up/15-down queries from the same 802-gene universe, seed 20261019. Apply identical aggregation, calculate one-sided empirical p-values with a +1 correction, and BH-adjust across all 544 scored ingredients. The null matches gene count and direction balance; it does not preserve gene correlation, pathways, expression magnitude or disease causality.
7. The primary follow-up rule, coded before scoring, requires q ≤ 0.10, at least two cell models, positive scores in at least 75% of cells, and positive aggregate scores for all sensitivity queries. This yields 29 compounds. It is a computational screening rule, not a therapeutic advancement gate.

## Post-hoc matched-drug check

After inspecting the primary results, we added a separate sensitivity check for common transcriptional responses. For each exact cell/time/dose/unit context, compare a drug with the median of other eligible drugs, requiring at least 20 distinct other ingredients and excluding the queried ingredient. Aggregate the residual score and background percentile within cells and then across covered cells.

Among the original 29, require at least two covered cells, positive mean residual and mean background percentile ≥ 75%. Twenty satisfy this additional rule. The original 29 flags and all background coverage counts remain in the output; the post-hoc check does not replace or recalibrate the primary q-values. Coverage differs between compounds and sparse contexts can fail for lack of comparators.

This is a common-response sensitivity check, **not a cytotoxicity or cytostasis adjustment**. No matched growth/viability measurements were integrated. A drug can reverse this signature by changing proliferation, cell state or survival without rescuing BUB1B dysfunction. The random-gene null and background check cannot resolve that ambiguity. The next useful analysis is exposure/biological-mechanism triage of these measured hypotheses, not clinical use or a claim that 20 good drug candidates were found.

## Artifacts and reproduction

Outputs: `runs/gse22206_lincs_20261019/analysis/`.

- `compound_scores.csv`: all 544 compounds, five query scores, q-values, coverage and both follow-up flags.
- `computational_followup.csv`: 29 primary hits, including the flag distinguishing the 20 background-supported results.
- `signature_scores.csv`, `condition_scores.csv`, `cell_scores.csv`: retained experimental context and scores.
- `matched_drug_background.csv`: exact-condition comparator counts, residuals and percentiles.
- `identity_mapping.csv`: exact FDA/LINCS identity links.
- `selected_landmark_matrix.npz`, `null_scores.npz`: extracted expression block and compound null distributions.
- `summary.json`, `manifest.json`: counts, decisions, hashes and provenance. The selected matrix and identifiers are archived and hashed; the full 23 GB GCTX was not rehashed, and its size/mtime are recorded explicitly.

```bash
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 python pipeline/score_gse22206_lincs.py \
  --query runs/gse22206_20261019/analysis \
  --lincs ../gene_mol/data/lincs_raw \
  --fda runs/fda_repurposing_20261019/final/fda_ingredients.parquet \
  --output runs/gse22206_lincs_20261019/analysis
python -m unittest discover -s tests -p test_lincs_reversal.py -v
```

Requires numpy, pandas, pyarrow and h5py. Tests cover reversal direction, balanced query weights, rejecting one-sided queries, multiple-testing order and excluding the queried drug/unmatched contexts from background calculations. All four tests passed. No inference API, new model training or laboratory work was required.
