# FDA drug recovery benchmark

## Scope and result

This is a fixed retrospective **development pilot** of evidence retrieval and ranking. It does not evaluate the local LLM, patient-specific efficacy, mechanism direction, exposure feasibility or MVA rescue. No weights were fitted. All six disease queries and scoring rules were fixed before the first completed scoring run.

The 2026-10-19 run ranks **1,443 currently eligible mapped FDA ingredients** for cystic fibrosis, rheumatoid arthritis, type 2 diabetes, asthma, Parkinson disease and familial hypercholesterolemia. The evaluation contains **165 unique drug–disease pairs** from Open Targets approval-indication metadata. Other candidates are unlabelled, not confirmed negatives. This includes component-level and symptomatic uses; recovering a component does not validate a combination regimen or disease correction.

| Method | Macro recall@10 | Macro hit@10 | MRR |
|---|---:|---:|---:|
| Target overlap | 21.6% | 85.8% | 0.629 |
| Measured assay overlap | 17.8% | 99.1% | 0.574 |
| Chemical similarity | 13.5% | 83.3% | 0.354 |
| Fixed equal-weight combination | 19.5% | 100% | 0.792 |
| Target-count baseline | 0.5% | 16.7% | 0.179 |
| Random expectation | 0.7% | 17.2% | 0.075 |

Recall measures the fraction of known positives retrieved, whereas hit rate asks whether at least one is retrieved. Metrics average diseases equally and integrate uniformly over tied rankings. A higher hit rate alongside lower recall means finding one answer more consistently while recovering fewer of the known answers. The combined score did not beat simple target overlap on recall@10; it must not be described as an overall improvement.

Two source-checked sentinels make failures concrete. Ivacaftor ranks **3/1,443** for cystic fibrosis under the combined score, while metformin ranks **1,396/1,443** for diabetes. Ivacaftor's indication requires a responsive CFTR mutation; broad disease recovery does not test that eligibility. Metformin is a glycaemic-control treatment, not an example of correcting a constitutional genetic lesion. Sources: [Kalydeco label](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=0ab0c9f8-3eee-4e0f-9f3f-c1e16aaffe25&version=36), [metformin label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=53bde438-4277-f418-e063-6394a90a5ea2).

## Leakage controls and remaining limits

- `prepare` places indication labels and sentinel identities in `truth/`. `score` reads only an allowlisted `features/` directory. Unexpected files, columns, query fields and changed hashes are rejected. Tests mutate external truth and require identical scores.
- Disease inputs use only the direct `genetic_association` datatype. Overall Open Targets scores, `known_drug`, literature, trial phases, drug descriptions and indication text are excluded from scoring. Open Targets distinguishes these evidence types in its [documentation](https://platform-docs.opentargets.org/evidence).
- The full eligible catalogue remains in each ranking, including drugs with no mechanism or assay coverage. Labels do not select candidates or target genes. Parent mappings use the existing ChEMBL hierarchy, not fuzzy name matching.
- This is **not a publication, temporal, disease-family or chemical-scaffold holdout**. Present-day genetic and mechanism annotations may reflect research on already established treatments. It tests recovery with direct indication links hidden, not prospective discovery. These six diseases are now development cases, not a reusable unseen test set.
- FDA ingredient/product eligibility does not prove every metadata indication is an FDA-approved use. Apart from the two sentinel checks, the label set is a metadata reference, not a fully adjudicated clinical gold standard. There were 333 raw approval rows, of which 123 were outside the mapped eligible universe; the remaining rows collapse to 165 parent-level pairs. Those exclusions must not be silently counted as successes or ignored in coverage reporting.
- Biologics and other unmapped ingredients can be missing from the catalogue. No specificity, clinical precision or AUROC is reported because unlabelled candidates are not established failures.
- A nonzero fingerprint similarity is easy to obtain. The reported positive-score coverage of similarity/combined methods is numerical coverage, not biological support.

## Methods

The fixed pilot selects up to 50 genes per disease with genetic-association score at least 0.1. These are operational limits, not biological truth. The target baseline takes the maximum genetic score over a drug's annotated mechanism targets. The assay baseline uses the same aggregation over qualifying unflagged original human single-protein ChEMBL measurements, reusing the FDA review's censoring and quality rules. Activity is not interpreted as beneficial direction.

Chemical similarity takes the maximum Morgan radius-2, 2,048-bit Tanimoto similarity to genetically linked mechanism-target anchors, weighted by the anchor's genetic score. Anchors are selected without indication labels. Self and identical InChIKey connectivity blocks are excluded. This baseline depends on mechanism annotations and is not an independent modality. The combined score is the fixed mean of target, assay and similarity scores. A disease-agnostic mechanism-target count and exact random expectation provide additional baselines.

## What failed and what to implement next

For the 165 metadata-positive pairs, the target route has:

| Target-route status | Pairs |
|---|---:|
| Selected target support | 20 |
| No direct genetic support for annotated targets | 105 |
| No mechanism annotation | 19 |
| Support excluded by the 50-gene limit | 19 |
| Support below the genetic threshold | 2 |

These categories describe database/selection coverage, not biological ineffectiveness. The dominant failure is relying on direct genetic-target overlap for treatments that can act downstream or symptomatically. Even in diseases with many approved treatments, few receive direct support through this route.

The next development change should compare an uncapped genetic-target route and a separately labelled downstream/pathway route, with explicit direction and symptomatic-versus-corrective intent. Retain this pilot as development data. Before tuning, freeze a larger, disease-family-separated evaluation set, adjudicate indication/subtype labels, and check incremental recovery against these baselines. Do not tune a special rule to force metformin into the top ten.

## Reproduce

Run from the repository root using the chemical environment (pandas, pyarrow, numpy and RDKit). Replace paths if the sibling repositories are located elsewhere. Outputs must be fresh directories.

```bash
../gene_mol/.venv/bin/python pipeline/benchmark_recovery.py prepare \
  --config config/recovery_benchmark.json \
  --kg ../mva_llm/data/kg \
  --fda runs/fda_repurposing_20261019/final \
  --chembl ../gene_mol/data/chembl37_db/chembl_37/chembl_37_sqlite/chembl_37.db \
  --out runs/recovery-example/prepared

../gene_mol/.venv/bin/python pipeline/benchmark_recovery.py score \
  --features runs/recovery-example/prepared/features \
  --out runs/recovery-example/scores

../gene_mol/.venv/bin/python pipeline/benchmark_recovery.py evaluate \
  --prepared runs/recovery-example/prepared \
  --scores runs/recovery-example/scores \
  --out runs/recovery-example/evaluation

../gene_mol/.venv/bin/python pipeline/benchmark_recovery.py diagnose \
  --prepared runs/recovery-example/prepared --kg ../mva_llm/data/kg \
  --out runs/recovery-example/diagnostics

../gene_mol/.venv/bin/python -m unittest discover -s tests -p test_recovery_benchmark.py -v
```

The recorded run is `runs/recovery_benchmark_20261019/`. Use `prepared/`, `scores/`, `evaluation/` and `diagnostics/`. Original measurements, selected feature tables, label mapping exclusions, per-disease metrics, positive rank intervals, sentinel ranks and SHA-256 manifests are retained. The full 30-GB ChEMBL database is identified by path/size/mtime in this run; extracted original records are hashed. No private patient data are used.

The proposed uncapped and pathway extensions have now been implemented and
evaluated on five frozen additional disease families. Neither improved the
primary endpoint enough to pass the prespecified promotion rule. See the
[pathway comparison and limitations](PATHWAY_BENCHMARK.md).
