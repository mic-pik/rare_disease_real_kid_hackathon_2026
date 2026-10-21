# FDA-approved candidate retrieval

This separate CPU workflow implements the first batch of the repurposing plan:
baseline preservation, a candidate record contract, FDA ingredient mapping,
measured ChEMBL retrieval and mandatory related-disease retrieval. It uses
public resources only and makes no case-model calls. It does not change the
historical hackathon results or claim that retrieved drugs rescue MVA.

## Run

Use the chemical Python environment (pandas and pyarrow). Download a snapshot
from the [official Drugs@FDA data page](https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files)
to `inputs/drugsatfda.zip`. The archive is updated frequently: retain your
download and its run-manifest hash for reproducibility.

Copy `config/fda_repurposing.example.json` to
`config/fda_repurposing.private.json`, then set its three input paths to your
FDA ZIP, original ChEMBL SQLite database and sibling-disease parquet table.
Resource paths resolve relative to the configuration file. Edit the explicit
gene and disease sets for the analysis; these are search scope, not claims
that inhibition of those targets is beneficial.

```bash
gene_mol/.venv/bin/python pipeline/fda_repurposing.py \
  --config config/fda_repurposing.private.json \
  --out runs/fda-repurposing-001
```

Use a fresh output directory. Existing outputs are never overwritten. Inputs,
including the original ChEMBL database, are hashed; this can take time for
large files. SQLite is opened read-only with temporary selection tables.
No resources are downloaded by the runner.

## Identity and eligibility

The FDA parser retains every product and its application, original approval
dates, ingredient text, formulation/route, strength and marketing status.
Eligibility requires original `ORIG/AP` approval evidence and a current
prescription or OTC product status. Tentative, unknown, manufacturing-only,
discontinued-only and conflicting statuses do not enter the active pool.
Discontinuation is not interpreted as a safety withdrawal. Safety-related
withdrawal review is explicitly `not_assessed` in this batch.

Names match ChEMBL preferred names or synonyms exactly after case/whitespace
normalisation. Salt forms collapse only through explicit ChEMBL parent links;
ambiguous parents and missing structures are quarantined. This favours
precision over coverage and is not comprehensive coverage of biologics.
Combination products retain all components, including unresolved ones.
Approval as a combination component does not imply standalone approval:
`eligible_single_ingredient_product_ids` records that distinction.

## Outputs and contract

| File | Meaning |
|---|---|
| `baseline.json` | Repository revision, tracked-file hashes, original candidate hashes/verdicts |
| `baseline_candidates.json` | Historical candidates in the new contract with complete original results preserved |
| `fda_products.parquet` | All product-level evidence, including excluded statuses |
| `fda_ingredients.parquet` | Canonical identities, FDA aliases and product references |
| `identity_mapping_issues.csv` | Unmatched, ambiguous or structureless ingredient names |
| `measured_interactions.parquet` | Original assay rows, including conflicts and censored values |
| `target_coverage.csv` | Human single-protein target mappings for the requested genes |
| `related_disease_candidates.parquet` | FDA-intersected related-disease indication metadata, retaining intent |
| `candidates.json` | Deduplicated union, sorted by identity, with evidence-artifact pointers |
| `candidates.csv` | Compact review table, explicitly distinguishing standalone product availability |
| `summary.json` | Counts, missing target mappings, measurement gaps and disease coverage gaps |
| `manifest.json`, `status.json` | Input/output hashes, configuration, limitations and execution status |

Each candidate has `schema_version`, `identity`, `retrieval_routes`, `evidence`,
`mechanism`, `assessments`, `opposing_evidence`, `unresolved_assumptions`,
`decision`, `decision_reason` and `clinical_suitability`. All newly retrieved
candidates remain `insufficient_evidence`. No ordering is a therapeutic rank.
Historical candidates retain the exact original gate under `original_result`.

Measured retrieval retains binding/functional assays with confidence >=8
against human single proteins. It includes all available activity types,
weak results, validity flags, duplicate flags and source identifiers. It does
not turn a binding result or an IC50 into an inferred mechanism of action.
Parent mapping preserves the actual tested compound identifier, so a salt's
measurement remains attributable. It does not substitute active metabolites
for their parent drugs.

Related-disease rows remain `indication_metadata_only`; their measured endpoint
is `not_retrieved`. Curated intent is not independently validated efficacy.
Primary-paper retrieval and review must precede any biological nomination.
The metformin/Fanconi record is a coverage check, not an MVA recommendation.

## Validation and next batch

```bash
gene_mol/.venv/bin/python -m unittest discover -s tests -v
python3 -m unittest discover -s mva_llm/unified/tests -v
python3 pipeline/check_repository.py
```

Synthetic tests exercise approval/status separation, salt deduplication,
ambiguous-name quarantine, incomplete combinations, censored/conflicting assay
records, unknown action direction, sibling parent joins and non-nomination.
These test data integrity, not clinical performance or held-out retrieval recall.

Remaining work: reviewed identity overrides; primary endpoint extraction;
ExCAPE provenance-preserving retrieval; signed regulator paths; exposure and
safety review; structured ranking; a publication-held-out benchmark; then
qualified LINCS/JUMP branches. Jev is not needed to run this first batch.

## First public-resource run (19 October 2026)

Using the FDA archive downloaded on 19 October, the existing ChEMBL 37
database and the local sibling-disease table, the initial retrieval produced:

| Measure | Count |
|---|---:|
| FDA product rows retained | 51,889 |
| Distinct FDA ingredient names | 2,909 |
| Exact names mapped to structured ChEMBL parents | 2,193 |
| Canonical ingredients with current product evidence | 1,443 |
| Original assay records for the configured target set | 858 |
| Related-disease indication records | 18 |
| Retrieved canonical ingredients | 131 |
| Retrieved ingredients with current standalone products | 128 |

The three combination-only ingredients are niacinamide, defactinib and
avutometinib. Their product evidence must not be described as standalone
approval. Metformin/Fanconi anaemia was recovered through its salt-to-parent
mapping, retaining the `instability_directed` intent and an unverified endpoint.
All candidates remain `insufficient_evidence`.

There are 716 quarantined ingredient names: 566 lack a usable parent structure,
145 are unmatched and five match multiple parents. These are coverage gaps,
not evidence against those ingredients. MAD2L1 has no qualifying target
mapping. Bloom and Werner syndromes have no eligible related-disease records
in this run; this is not a claim that no relevant studies exist.

The automated checks comprise nine FDA contract/integration tests, five
portability tests and the existing 28 workflow tests. The integration test
checks ZIP ingestion, original assay preservation, historical verdicts,
output hashes and refusal to overwrite a completed run. These results do not
measure prospective rescue efficacy, clinical suitability or benchmark recall.

## Evidence-review pass

The next CPU stage verifies upstream hashes, interprets individual assays and
combines those results with explicit source-level curator dossiers:

```bash
gene_mol/.venv/bin/python pipeline/review_fda_candidates.py \
  --input runs/fda-repurposing-001 \
  --review config/fda_review_evidence.json \
  --out runs/fda-review-001
```

The versioned evidence JSON contains five reviewed candidates, source URLs,
source types/locations, opposing findings, exposure limitations and proposed
decisive experiments. These are assistant-curated judgments, not an independent
expert assessment or a benchmark of autonomous reasoning. Remaining candidates
receive assay/route triage only; their clinical exposure and safety are marked
unassessed. Do not interpret an alphabetical position within a tier as an
efficacy rank.

Potency thresholds are operational: <=1 micromolar is a retrieval activity
signal; >=10 micromolar is weak at the screening threshold. Intermediate or
nondecisive censored bounds remain ambiguous. Units and relation operators are
preserved. Zero concentration cannot become infinite potency. Percent assays
require interpretable polarity and an unambiguous test concentration. Generic
fold changes, phosphorylation and other indirect endpoints need manual review.
Binding does not establish functional direction. Flagged duplicates, invalid
measurements and contradictory source comments cannot contribute to positive
support counts. Mixed active/weak measurements remain visible by target.

Outputs include `review_queue.csv` for all ingredients,
`experimental_shortlist.csv`, `assay_triage.parquet`, `target_evidence.csv`,
five source-linked Markdown dossiers, `report.md`, preserved review inputs and
hash manifests. Rejection means the retrieved rationale does not justify this
rescue objective, not that every other use or mechanism is disproved.

On the first 131-ingredient pool, the pass classified 858 records as 62 active,
482 weak/inactive at the screening threshold, 115 ambiguous and 199
uninterpretable. These counts precede exclusion of 87 quality-flagged records
from support counts. Metformin is a conditional experiment priority;
amifostine is deferred; niacinamide, sirolimus and ceritinib have rejected
retrieval rationales documented in their dossiers. No MVA rescue or clinical
lead is established. The shortlist is not padded to an arbitrary size.

Nine additional tests cover censoring, unit conversion, nonpositive potency,
percent polarity, ambiguous concentrations, indirect endpoints, quality flags,
mixed evidence and rejection of source-free curator priority. A full
label/withdrawal review, a publication-held-out benchmark and experiments in a
relevant model remain necessary before stronger conclusions.

The subsequent [stress-phenotype assessment](STRESS_PHENOTYPE_ASSESSMENT.md)
did not establish the required phenotype in the available model data and found
primary aneuploid-cell antiproliferation counterevidence. The current evidence
configuration therefore defers metformin; the regenerated review in
`runs/fda_review_20261019/stress_followup_verified/` has no experiment-priority candidates.
The original `analysis/` output remains a historical snapshot. The assessment
specifies baseline measurements and functional rescue criteria for reopening
this hypothesis.
