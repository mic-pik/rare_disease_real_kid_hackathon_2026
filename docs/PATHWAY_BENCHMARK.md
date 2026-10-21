# Uncapped genetics and pathway retrieval: frozen comparison

## Outcome

Implemented and evaluated on 2026-10-19. **Neither change met the prespecified promotion rule.** The existing retrieval baseline and MVA candidate decisions remain unchanged. The new methods are available as experimental alternatives, not validated improvements.

The original six diseases remain development cases. Five new families—glaucoma, osteoporosis, peptic ulcer disease, acute myeloid leukemia and chronic kidney disease—were frozen before scoring. They contain 108 unique approval-metadata pairs in the same 1,443-ingredient FDA universe.

| Split | Method | Macro recall@10 | Macro hit@10 | MRR |
|---|---|---:|---:|---:|
| Development | Top-50 genetic target overlap | 21.58% | 85.77% | 0.629 |
| Development | Uncapped genetic target overlap | 21.41% | 83.33% | 0.622 |
| Development | Uncapped plus pathways | 21.41% | 83.33% | 0.623 |
| Held out | Top-50 genetic target overlap | 4.99% | 45.02% | 0.273 |
| Held out | Uncapped genetic target overlap | 4.79% | 39.84% | 0.257 |
| Held out | Uncapped plus pathways | 4.79% | 39.84% | 0.253 |
| Held out | Pathways alone | 2.81% | 19.62% | 0.071 |
| Held out | Random expectation | 0.69% | 13.95% | 0.063 |

Held-out recall differences versus the capped baseline were −0.194 percentage points for both changes (paired disease-bootstrap interval −0.472 to 0.000 percentage points). Pathway augmentation versus uncapped genetics had zero recall difference in every disease. These are exploratory intervals from only five diseases, not proof of equivalence or general inferiority.

Both changes broaden coverage. Held-out macro positive-score coverage increases from 16.7% for the capped baseline to 17.3% for uncapped targets and 46.0% with pathways. This means more known treatments receive some numerical support; it does not establish correct mechanism or improve the top-ten shortlist. The small recall decline partly reflects removal of chance credit from zero-score ties, rather than losing previously supported top-ten treatments.

## Implemented changes

1. `benchmark_recovery.select_genetics(..., maximum=None)` retains every qualifying genetic association instead of the first 50. The 0.1 minimum genetic score remains fixed. The old capped behavior remains available and unchanged.
2. `pipeline/benchmark_pathways.py` tests the capped target baseline, uncapped target overlap, pathway membership alone, pathway augmentation, target count and random expectation. No weights are fitted.
3. Each pathway connection records a disease-associated seed gene, a different drug-target gene, pathway identifier/name, pathway size, annotated drug actions and direction status. Self-connections are excluded. The score is the seed genetic score divided by `sqrt(pathway_size - 1)`. Augmentation takes `max(direct_score, 0.5 * pathway_score)`.
4. Only pathways with 2–100 distinct genes are used. If any annotation places a pathway in `Disease` or `Drug ADME`, the whole pathway is excluded to reduce direct disease/treatment annotation shortcuts. Duplicate memberships cannot inflate scores. These controls do not eliminate publication overlap or all annotation bias.
5. Direction constraints require a disease, target gene, required increase/decrease, source and source locator. Opposing or mixed actions remain visible in the audit and cannot supply a pathway-expansion score. Unknown requirements remain explicitly unknown. A compatible action is still not functional rescue.

## Direction is an unresolved data requirement

The local table contains unsigned gene–pathway membership, not signed causal edges. A shared pathway does not establish that the drug target is downstream, whether it increases or decreases the relevant function, or whether changing it compensates for the lesion. Reactome distinguishes membership/enrichment analysis from pathway topology in its [pathway browser documentation](https://reactome.org/userguide/pathway-browser).

No independently sourced required-direction constraints were available for this frozen test. **All 17,141 recorded drug-target pathway connections are `unknown_required_direction`.** The action-checking code is exercised with synthetic tests, but this run does not claim biologically validated direction checks. No route is labelled therapeutic compensation, and no MVA candidate is promoted.

The compatible/opposing check applies to pathway expansion. Direct-target scores remain the unsigned benchmark baselines and are not clinical nominations. Broad disease labels can combine causal, symptomatic, supportive and combination uses, so successful indication recovery would still not establish correction of the underlying defect.

## Split and leakage controls

`config/pathway_benchmark.json` was copied to `runs/pathway_benchmark_20261019/frozen/` with a timestamp and SHA-256 before implementation scoring or outcome inspection. Only ontology structure and raw approval-label availability counts were inspected during case selection. No query was selected based on ranking performance.

The family roots are explicit in the configuration. Preparation checks membership in each chosen family and rejects a query lying within an opposite-split family or disease branch. These are selected ontology families, not disjoint therapeutic areas, drugs, pathways or genetic mechanisms. For example, the peptic-ulcer family is `ulcer disease`; the broader digestive-system root would overlap the diabetes development case and was not used. Family audit records preserve ancestors and root names.

Scoring reads an allowlisted feature directory and never reads approval labels. Tests change external labels and require identical scores. Disease genetics excludes known-drug and literature association types. The same candidate and mechanism snapshots are reused from the earlier benchmark. All 8,658 development capped-baseline scores exactly reproduce the previous run.

Executable source snapshots are captured at preparation, and evaluation rejects a changed main implementation. Frozen configuration, features, truth tables, scores, full rankings and reports have manifests. This is not a temporal or publication holdout: present-day genetic and mechanism annotations may reflect established treatment research. The held-out cases have now been evaluated and should not be reused as an unseen test set after tuning.

The primary endpoint was macro recall@10. Pathway promotion required a positive mean paired difference and positive lower 95% paired-bootstrap bound against both capped and uncapped target baselines. Uncapped promotion similarly required a positive lower bound against capped targets. Neither passed. Results were not used to retune size limits, score thresholds or pathway weight.

## Reproduction

Run from the repository root with the chemical environment. To start a new experiment, copy the configuration into a fresh frozen directory and create `freeze.json` containing `created_utc`, `config_sha256` and the freeze rationale **before scoring**. The recorded frozen directory can be reused for an exact rerun; output directories must be fresh.

```bash
../gene_mol/.venv/bin/python pipeline/benchmark_pathways.py prepare \
  --frozen runs/pathway_benchmark_20261019/frozen \
  --base runs/recovery_benchmark_20261019/prepared \
  --kg ../mva_llm/data/kg \
  --ontology ../mva_llm/data/raw/opentargets/disease \
  --chembl ../gene_mol/data/chembl37_db/chembl_37/chembl_37_sqlite/chembl_37.db \
  --out runs/pathway-example/prepared

../gene_mol/.venv/bin/python pipeline/benchmark_pathways.py score \
  --features runs/pathway-example/prepared/features \
  --out runs/pathway-example/scores

../gene_mol/.venv/bin/python pipeline/benchmark_pathways.py evaluate \
  --prepared runs/pathway-example/prepared --scores runs/pathway-example/scores \
  --out runs/pathway-example/evaluation

../gene_mol/.venv/bin/python -m unittest discover -s tests -p test_pathway_benchmark.py -v
```

Recorded results are in `runs/pathway_benchmark_20261019/evaluation/`. Full pathway explanations are in `scores/pathway_routes.parquet`; family checks and executable snapshots are in `prepared/`. `decision.json` explicitly rejects promotion of both experimental methods.

## What this resolves

Removing the cap and adding unsigned pathway neighbours are insufficient to improve prioritisation in this test. These methods can generate a wider review queue, but are not a justified replacement for the baseline. The next evidence requirement is a signed, context-specific mechanism or functional response, with corrective and symptomatic treatment objectives evaluated separately. Increasing pathway weight after seeing these held-out outcomes would be development tuning and would require a new untouched evaluation set.
