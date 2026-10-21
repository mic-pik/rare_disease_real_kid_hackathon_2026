# Functional-characterization intake templates

Use with the [study handoff](../../docs/FUNCTIONAL_CHARACTERIZATION.md). These CSVs are header-only templates, not results. Copy to ignored `inputs/functional_characterization/` before entering real data. Use pseudonymous identifiers and retain source evidence privately.

- `models.csv`: one row per model/clone. `model_role` identifies deficient, corrected, unrelated healthy or engineered comparator. Record unknown genotype/phase explicitly; a prediction is not confirmation. References point to evidence rather than personal identifiers.
- `cultures.csv`: one row per independently established culture and assigned arm. `culture_id` must be unique; arms split from a common preparation share `batch_id`. Record the actual replication design separately. Multiple wells from one culture are technical replicates, not independent cultures. Aggregate only disjoint observed cells, with consistent windows and no double-counting tracks across fields. `arm_id` maps to the laboratory's documented intervention and achieved expression, not merely a presumed successful rescue.
- `measurements.csv`: one row per measured endpoint and sampling unit, linked to `culture_id`. Retain well and assay identifiers, challenge status, units and normalization references. Do not combine arbitrary fluorescence units across acquisition settings without calibration. Preserve single-cell measurements in referenced raw files.

Culture count definitions:

| Field | Definition |
|---|---|
| `initial_live_cells`, `final_live_cells` | Live counts in the same prespecified tracked population at window start/end; document detachment and field entry/exit. |
| `completed_divisions` | Completed division events during `window_hours`, under a prespecified completion definition. Descendant divisions may occur; this is an event count. |
| `error_divisions` | Completed division events with at least one prespecified segregation error; between zero and `completed_divisions`. |
| `aborted_divisions` | Events entering mitosis without successful completion under the scoring definition; keep distinct from completed events. |
| `death_events` | Death events in the tracked population during the window, including daughters when applicable. This is not automatically a fraction of initial cells because descendants can die. |
| `lost_tracks` | Tracks whose fate cannot be scored, with reasons retained in raw data. Loss is not death or successful division. |
| `daughters_eligible` | Daughters from a prespecified division cohort eligible for fixed-duration follow-up. Define how overlapping generations are handled. |
| `daughters_surviving`, `daughters_dead`, `daughters_censored` | Mutually exclusive daughter outcomes at `daughter_followup_hours`; their sum equals `daughters_eligible`. Censored outcomes remain unknown. |

All counts are nonnegative integers when measured. Empty fields mean missing/not measured, never zero. Numeric zero means an observed zero. `window_hours` must be positive. The error fraction is undefined when no divisions complete. Report division events per initial live cell per hour only when both denominators are positive. Report daughter survival among observed outcomes together with censoring and cohort counts; do not silently drop censored daughters.

Retain excluded rows with reasons. Document randomization, blinding, culture independence, raw track identifiers, endpoint definitions and analysis-plan version alongside the CSVs. These templates do not automatically validate data or declare rescue.
