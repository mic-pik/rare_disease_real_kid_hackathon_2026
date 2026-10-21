# Secondary stress phenotype: evidence and next experiment

Assessment date: 2026-10-19. Assistant-curated research assessment, not independent expert validation.

## Decision

The inspected data do not establish a metabolic or oxidative-stress phenotype in the relevant BUB1B-deficient model. This is missing evidence, not evidence that stress is absent. Metformin moves from conditional experiment priority to `defer_transfer`. It remains a retrieved FDA-approved ingredient, but no candidate currently passes the curated experiment-priority gate. Earlier review artifacts remain historical snapshots.

## What the local data can establish

| Source inspected | Available evidence | Limitation for this question |
|---|---|---|
| `rare_real/data`, patient work reports | WGS reads, variants, clinical phenotype document, computational variant/discovery reports | No matched cell-culture ROS, respiration, protein, RNA-seq or segregation-rescue measurements identified in the inspected material. Variant prioritisation does not establish a metabolic phenotype. |
| `gene_mol/work/rescue/rescue.json` | BUB1B knockout: five wells, signature strength 1.58, stronger than fraction 0.048 of genes; explicitly `usable: false` | Weak morphology cannot establish absence or presence of oxidative stress. |
| `gene_mol/src/rescue_candidates.py` | Separate compound morphology compared with BUB1/MAD2L1 knockout proxies; pooled PRISM viability filter | Neither a BUB1B-specific stress measurement nor drug treatment of deficient cells. Pooled cancer-line viability does not demonstrate preserved growth in the intended model. |
| `gene_mol/data/lincs_raw/GSE92742…sig_info.txt.gz` | 121 BUB1B rows: 89 individual shRNA, 14 consensus shRNA, 18 overexpression; 36 metformin rows | These are signature rows, not 121 independent loss-of-function experiments. Consensus and constituent shRNAs must not be counted independently. Overexpression is not deficiency. |
| `mva_llm/data` and pipeline reports | Literature, knowledge graphs, training examples and generated reasoning | Retrieved or generated claims are not new patient/model measurements. |

The LINCS metadata share ten cell identifiers across all BUB1B perturbation types and metformin. Metformin doses recorded in `pert_dose` are 0.15625, 0.625, 2.5 and 10 micromolar, at 6 or 24 hours. BUB1B perturbations span 24–168 hours. Shared identifiers do not establish matching perturbation type, exposure time or biological state. Rounded `pert_idose` labels differ from exact `pert_dose`; preserve the latter. No expression-matrix analysis or stress-pathway enrichment was performed in this audit.

The selected metadata, structured counts and input hashes are saved in `runs/stress_phenotype_20261019/`. This inventory is bounded to the repositories and artifacts inspected; external laboratory records may exist.

## Primary literature changes the interpretation

Experimental aneuploidy can trigger a ROS-dependent ATM–p53 response that constrains proliferation. That supports measuring stress but does not establish it in this model. It also means suppressing a stress signal is not automatically beneficial. [Li et al., 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2922543/).

Metformin impaired accumulation of trisomy 13 and 16 mouse cells in a primary study. This raises a growth/selection concern for a rescue objective. It does not establish clinical toxicity or prove the same effect in BUB1B deficiency; supplementary dose details were not independently extracted here. Do not attribute the paper's stronger AICAR/17-AAG results to metformin. [Tang et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3532042/).

The Fanconi study supports a transfer hypothesis, but its single-arm design and lack of significant cohort-level chromosome-breakage or micronucleus improvement do not demonstrate correction of spindle-checkpoint failure. [Pollard et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9631552/).

## Implementable measurement plan

1. **Establish the model.** Identify available patient-derived or engineered deficient cells, matched corrected controls, and independent cultures/clones. Confirm genotype, residual BubR1 abundance and baseline segregation phenotype. Keep uncertain allele phase/function explicit; a knockout surrogate is not automatically equivalent to the case. If these cells or records are unavailable, stop biological advancement at this gate.
2. **Measure baseline stress before adding a drug.** Compare deficient and corrected cells under matched passage, density and media conditions. Use an oxidative-stress readout plus an independent damage/redox measurement; add respiration/energy measurements if the proposed mechanism is metabolic. Record live-cell counts, division rate and death alongside each measurement. A single fluorescent signal or inferred transcriptional signature is insufficient.
3. **Define the analysis before screening.** Use independent cultures as replication units, with cells nested within culture; do not treat thousands of imaged cells as thousands of independent replicates. Randomise wells, blind scoring, specify exclusions and estimate variance in a pilot before fixing sample size. Set biologically meaningful improvement and allowable growth-loss margins before viewing drug effects; no universal margin is justified by the current data.
4. **Only after baseline relevance is supported, test the interaction.** Compare vehicle and metformin in both deficient and corrected cells over a justified exposure range, with measured exposure where feasible. Include a growth-slowing comparator matched to any observed proliferation reduction. Track completed divisions, abnormal segregation events, surviving daughters and cell deaths. Short-term viability alone cannot detect cytostasis or selective depletion.
5. **Advance on functional rescue.** Require a reproducible reduction in segregation errors per completed division, preserved absolute survival and growth within the prespecified margins, and supporting stress improvement. Assess treated deficient versus vehicle deficient cells and the corresponding change in corrected cells. Confirm across independent cultures and a second relevant model if available. Stress-only improvement, fewer divisions, selective death or effects only at unjustified exposure remain deferred or rejected for this objective.

This is a study specification, not a patient dosing protocol. Plasma exposure is a reference, not proof of intracellular exposure. Do not choose a culture dose solely because it appears in LINCS or a different disease's paper.

## Optional computational work while cells are unavailable

A carefully controlled LINCS analysis could test whether BUB1B knockdown consistently induces stress-associated transcription in public models. Restrict the primary analysis to measured landmark genes, keep inferred genes separate, verify shRNA agreement and signature quality, and stratify cell line/time. Use matched random-gene and cytostatic perturbations as nulls. Any subsequent metformin reversal must be labelled a comparison of separate perturbations. Even a positive result would support a surrogate hypothesis only and would not satisfy the deficient-cell rescue gate.

The next dependency is therefore an inventory of available cells and baseline functional measurements, not another unrestricted compound ranking.
