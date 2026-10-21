# BUBR1 functional characterization: experiment handoff

Prepared 2026-10-19. Status: proposed study; no cells acquired, experiments performed or laboratory work scheduled. Model availability is unconfirmed. This is a study design for a collaborating laboratory, whose validated assay and delivery procedures must supply the experimental protocols.

## Question and deliverable

Can restoration of functional BUBR1 improve chromosome segregation in the relevant deficient model without suppressing division or selecting against abnormal cells? If so, is the actionable defect primarily protein abundance, localization/attachment, checkpoint function, or a combination?

The deliverable is a model-specific phenotype and positive-control report. It determines whether a subsequent compound screen has a measurable rescue endpoint. It does not establish a drug candidate or a patient treatment.

The [rescue audit](BUBR1_RESCUE_AUDIT.md) supports this sequence: genetic complementation can restore checkpoint activity, whereas localization rescue has sometimes been demonstrated only as arrested-cell alignment. Alignment can coexist with defective attachments. [Suijkerbuijk et al., 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2887387/), [Xu et al., 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3654266/), [Shrestha et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5532248/).

## 1. Establish a usable model — laboratory and analysis leads

- Inventory available deficient lines, matched corrected lines, independent clones and unrelated healthy controls. Record cell type, provenance, genotype evidence, phase status, passage, authentication and culture quality. Keep donor identities in a private crosswalk.
- Prefer an isogenic corrected comparator; use multiple clones where feasible to assess clone effects. Unrelated healthy lines provide context but cannot isolate the effect of the variant. Engineered models can establish assay feasibility but do not validate rescue in the patient model.
- Resolve the genotype evidence separately. The current computational handoff does not establish the second allele, phase or residual protein function. Do not copy an inferred genotype into a confirmed field.
- Establish whether the cells divide often enough for segregation measurements in a practical observation window. Record passage, density and medium consistently. If the model lacks a reproducible defect, stop and revisit the model or endpoint before screening drugs.

Output: completed `models.csv`, source records and a model-selection note. Real measurements belong under ignored `inputs/functional_characterization/`; the repository contains only blank templates in `examples/functional_characterization/`.

## 2. Run a baseline and genetic-rescue pilot — laboratory lead

| Arm | Purpose |
|---|---|
| Deficient + matched delivery/empty-vector control | Disease-model baseline and delivery control |
| Corrected + matched delivery/empty-vector control | Reference phenotype under the same handling |
| Deficient + wild-type BUBR1 restoration | Test whether the phenotype is rescuable |
| Corrected + wild-type BUBR1 restoration, if feasible | Detect expression/delivery effects in the reference model |

Confirm restoration abundance and localization; target a range comparable to corrected cells using the laboratory's established method. Record achieved expression rather than calling every transduced cell rescued. Supraphysiologic expression is a confounder. Include an unmanipulated baseline if delivery itself alters division or survival. Analyze all assigned cells as the primary population; any expression-selected analysis is secondary and must report exclusions.

Pilot measurements, with matched acquisition windows:

1. **Protein and localization:** BubR1 abundance with an appropriate loading/reference measure, cell-cycle context, and kinetochore localization. Measure BUB1B transcript if useful to distinguish transcript depletion from protein instability. Compare mitotic cells at matched stages; differences in mitotic fraction can distort bulk abundance.
2. **Division fidelity:** blinded live-cell scoring of completed divisions, divisions with at least one prespecified segregation error, mitotic duration, aborted divisions and cells lost to follow-up. Prespecify error categories such as lagging chromosomes and bridges and retain event-level records. A division with several errors counts once in the primary numerator.
3. **Survival and proliferation:** initial/final live-cell counts, completed divisions per initial cell per observation time, deaths, and daughter-cell outcomes at a fixed follow-up. Report incomplete daughter follow-up separately. These measurements distinguish improved fidelity from cytostasis or selective loss.
4. **Mechanism:** a laboratory-validated checkpoint challenge with matched untreated controls, and attachment maturation/localization measurements such as Astrin recruitment or cold-stable fibers. Challenge/arrest assays are separate from spontaneous live-cell fidelity assays; they cannot substitute for them.

The pilot establishes feasibility and variance. It is not a powered efficacy claim. The laboratory must specify independent culture batches and acquisition capacity before execution; no sample-size claim can be justified from the current files.

## 3. Lock the confirmatory analysis — analysis and laboratory leads

Use the pilot to select one primary fidelity endpoint and a fixed observation window. Prespecify the minimum meaningful improvement and acceptable survival/division margins before confirmatory data collection. Use between-culture variance and the intended comparison to choose the number of independent cultures. Do not treat hundreds of imaged cells as hundreds of biological replicates.

- Randomize arm positions and acquisition order within batch, blind image scoring, and distribute arms across independent culture batches. Record donor/clone/culture/well identifiers so technical wells remain nested within cultures.
- Primary descriptive rate: `error_divisions / completed_divisions` for each independent culture and arm. A zero-division culture has an undefined rate, never a zero-error success. Always show the numerator, denominator, division rate, death counts and censoring alongside it.
- Compare deficient restoration against its delivery control, with the corrected comparator showing the reference range. Use paired batch comparisons where the design supports pairing. Preserve donor/clone structure in the statistical model; do not pool all cells into an independence-assuming test.
- Report absolute differences and uncertainty. A hierarchical/binomial analysis or culture-level analysis must reflect the actual experimental unit. Missing observations and exclusion reasons stay visible. Do not manufacture precision when there are too few independent cultures.
- Prespecify how checkpoint, attachment, survival and proliferation findings affect advancement. A lower error fraction caused by fewer divisions, increased death or checkpoint bypass does not meet the intended rescue criterion.

Output: signed-off analysis plan with sample size, endpoint definitions, margins, exclusion rules and the genetic positive-control acceptance criterion. These fields are deliberately unresolved until the pilot and laboratory capacity are known.

## 4. Choose the next experiment from the result

| Observation | Next action |
|---|---|
| Low protein; wild-type restoration improves fidelity and preserves division/survival | Investigate whether endogenous mutant protein retains function before testing abundance-enhancing compounds. Wild-type rescue alone does not establish mutant rescue. |
| Adequate abundance but abnormal localization/attachment; wild-type rescue works | Consider targeted PP2A-B56 recruitment as a mechanistic control, with inactive/matched construct controls and the same fidelity/checkpoint endpoints. |
| Alignment improves but fidelity, checkpoint or daughter survival does not | Record partial endpoint rescue; do not advance on alignment alone. |
| No rescue despite verified wild-type expression/localization | Check assay sensitivity, delivery effects, genotype, irreversible/secondary abnormalities and model relevance before concluding that BUBR1 biology is unrescuable. |
| No reproducible deficient-versus-corrected phenotype | Change the model or endpoint; there is no interpretable rescue screen yet. |
| Genetic rescue reproducibly meets the prespecified criteria | Design a small compound study with exposure justification and matched controls using the validated endpoint. FDA approval and target binding remain insufficient nomination criteria. |

Targeted recruitment is a conditional mechanistic experiment, not a general PP2A-activation recommendation. Excess phosphatase recruitment can impair checkpoint/error correction. [Primary balance study, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10577578/).

## Data handoff and completion criteria

Blank CSVs are provided for models, culture summaries and individual measurements. They contain no synthetic successes and no patient data. Keep raw images/tracks and their scoring provenance; aggregate CSVs alone cannot resolve disputed scoring. The template README defines units and denominators.

The characterization is complete only when the report includes verified model context, independent-culture results, genetic-control behavior, fidelity and survival/division outcomes, missingness, and a justified decision from the table above. At present, the next dependency is access to cells or existing functional measurements and a laboratory owner. Sequencing, separate perturbation signatures and pooled cancer viability data in the workspace cannot supply these measurements.
