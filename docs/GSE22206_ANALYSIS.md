# Public BUB1B expression data: GSE22206

Analysis date: 2026-10-19. Completed public-data retrieval and descriptive analysis; no laboratory access required.

## Result

We now have an actual patient-derived expression dataset and an exploratory signature that can be compared with the existing LINCS drug signatures. This is a different patient from the hackathon case. No drug has been nominated from this analysis.

| Item | Result |
|---|---:|
| Arrays / probes | 14 / 22,184 |
| Affected donors | 1, represented in fibroblasts and lymphoblastoid cells |
| Unaffected control donors | 1 fibroblast donor; 2 lymphoblastoid donors |
| Additional relatives | 2 lymphoblastoid donors, sensitivity comparisons only |
| Duplicate-culture Pearson correlations | 0.9921–0.9982 |
| Detection-filtered, unambiguously annotated genes | 8,235 |
| Cross-cell-type direction-consistent signature | 685 genes: 469 higher, 216 lower |
| Also direction-consistent against both relatives | 542 genes |
| Signature genes mapped to existing LINCS gene metadata | 566 |
| Directly measured LINCS landmark genes | 74: 59 higher, 15 lower |

The BUB1B probe is lower in the affected donor: fibroblast difference −0.09185 and mean lymphoblastoid difference −0.60498. All case/control culture comparisons for this probe have the same direction. These are differences on the submitted variance-stabilized scale, **not log2 fold changes or percentages of remaining protein**. The small fibroblast difference is particularly important when interpreting effect strength.

The GEO record describes a homozygous c.2386-11A>G BUB1B splice variant and two biological replicates of seven donor/cell-type samples. Duplicate cultures do not create independent affected donors. The study contains neither an isogenic corrected comparator nor compound-treated deficient cells. [GEO GSE22206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22206), [original study](https://www.nejm.org/doi/full/10.1056/NEJMoa1006565).

## Files acquired and analysis rules

Downloaded the GEO series metadata, submitted normalized matrix, non-normalized signal/detection table and GPL6104 probe annotations into `runs/gse22206_20261019/sources/`. Exact source URLs, input/output SHA-256 hashes, script hash and library versions are recorded in `analysis/manifest.json`. The raw image archive was not needed for this bounded analysis and was not downloaded.

The submitted matrix already underwent lumi variance-stabilizing transformation and robust spline normalization. We use it without a second log transform or normalization. This does not independently reproduce the authors' preprocessing. The non-normalized table supplies detection p-values, which are assay-detection statistics, not case-control differential-expression p-values.

1. Verify sample identities and order against descriptions in the raw detection table. Require complete finite expression data and valid detection values. Export sample metadata, distributions and correlations.
2. Compare fibroblasts only with fibroblasts. Compare affected lymphoblastoid cells separately with each of the two unaffected donors, then average their contrasts with equal donor weights. Relatives never enter the primary control mean.
3. A probe must be detected at p < 0.01 in both cultures of at least one donor in each cell type. Requiring detection in either condition retains transcripts reduced or induced in the case. Use unambiguous historical platform symbols; preserve all probe-level results and stable gene IDs.
4. Collapse eligible probes per gene by the median per array. Require every eligible probe and all case/control culture pairings to agree in direction across both cell types and both lymphoblastoid control donors. Require the mean lymphoblastoid case difference to exceed the difference between the two controls in absolute value.
5. Rank the retained genes by their smallest absolute mean difference across fibroblasts and the two separate lymphoblastoid comparisons. This prioritizes consistent effects. It is not an efficacy score or statistical significance measure.
6. Export comparisons against the relatives as sensitivity results. Agreement with relatives cannot establish causality; disagreement may reflect shared biology rather than invalidating a gene.
7. Map single, unambiguous Entrez IDs into the existing GSE92742 LINCS metadata. Exclude duplicate gene-ID collisions. Export measured landmark genes separately from the larger overlap containing inferred measurements.

No disease-associated p-values or FDR are reported: the one affected donor cannot support population-level inference. The consistency filters are transparent exploratory choices, not calibrated false-discovery controls. An effect-size sensitivity check retains 253, 38 and 5 genes at minimum absolute differences of 0.25, 0.5 and 1 submitted units, respectively. Many of the 685 changes are small.

## Interpretation and next computational use

The top-ranked increases include GSTM1, ANXA1, SCHIP1 and MGST2; decreases include CHURC1, GSTT1, ATF5 and SLC7A7. These are observed expression contrasts, not new drug targets. In particular, the large CHURC1 difference against unrelated controls is not direction-consistent against both relatives, illustrating why shared donor/family background matters. Historical probe annotation and possible genotype-dependent hybridization remain unresolved; probes were not remapped to a current genome.

The data support testing a **public single-case expression-reversal hypothesis** using the existing LINCS matrix. The directly measured query has only 74 genes and is imbalanced (59 up / 15 down), so a future comparison should balance directions, report gene coverage and use size/direction-matched null signatures. It should test separate tissue signatures and the relative-supported subset, stratify drug results by cell line/time/dose, and distinguish independent experiments from consensus signatures. This analysis has prepared the query but has not run drug reversal or pathway enrichment.

A reversal score would still be a comparison of separate experiments, not drug treatment of this donor's cells. Expression differences can reflect donor background, immortalization/culture, cell state or consequences of disease. Reversing a protective response could be undesirable. Neither this signature nor a future reversal score establishes a stress phenotype, rescue of chromosome segregation or suitability of an approved drug. Existing candidate decisions remain unchanged.

## Reproduce

Use a Python environment containing numpy and pandas. Download the four public files at the URLs recorded in `pipeline/analyze_gse22206.py` under the corresponding filenames, then run:

```bash
python pipeline/analyze_gse22206.py \
  --sources runs/gse22206_20261019/sources \
  --output runs/gse22206_20261019/analysis \
  --lincs-gene-info ../gene_mol/data/lincs_raw/GSE92742_Broad_LINCS_gene_info.txt.gz
python -m unittest discover -s tests -p test_gse22206.py -v
```

The LINCS metadata argument is optional. Outputs include all probe/gene contrasts, per-sample BUB1B values and detection scores, the full signature, top-100 lists in each direction, the relative-supported signature and the measured-landmark query. Raw data and generated outputs remain ignored by Git; analysis code, tests and this report are distributable.
