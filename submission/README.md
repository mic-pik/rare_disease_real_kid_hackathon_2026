# Submission reports

This folder contains the Track 1 submission files and the prepared Track 2 methods report. Evidence from the completed automatic run is retained privately and is not distributed here; see [Scope and integrity](#scope-and-integrity). The Track 2 report nominates NMN for preclinical investigation while retaining the recorded `review_required` outcome; it does not turn that outcome into a passed scientific gate.

The [repurposing addendum](working_variant_track2_repurposing_addendum.md) records the later public-expression screen and exposure/mechanism review. None of the four reviewed leaders justifies replacing NMN; NMN's competition eligibility remains unresolved. The addendum is prepared locally and has not been uploaded. The original report and the privately retained run artifacts are unchanged by this review.

The [subsequent fixed-candidate review](../docs/LINCS_FOLLOWUP.md) covers the remaining sixteen hits, cell-cycle exclusion, individual BUB1B shRNAs and contextual PRISM viability. It also finds no justified replacement and identifies poor agreement between the knockdown and patient-expression signatures.

| File | Contents |
|---|---|
| [working_variant_automatic_track1.csv](working_variant_automatic_track1.csv) | Submitted automatic Track 1 predictions |
| [working_variant_automatic_track1_report.md](working_variant_automatic_track1_report.md) | Submitted Track 1 methods report |
| [working_variant_automatic_track1_epcr_audit.json](working_variant_automatic_track1_epcr_audit.json) | Row selection, scores and EPCR calculations |
| [working_variant_track2_report.md](working_variant_track2_report.md) | Prepared Track 2 report and qualified NMN nomination |
| [working_variant_track2_repurposing_addendum.md](working_variant_track2_repurposing_addendum.md) | Subsequent FDA expression-reversal screen and decision against replacing the nomination with its four reviewed leaders |
| [artifact_manifest.json](artifact_manifest.json) | Original source locations relative to the workspace, original hashes, packaged hashes and transformations |

The Track 1 CSV and EPCR audit are byte-identical to the supplied originals. The Track 1 report has its clinical feature narrative generalised: individual proband phenotype terms are removed from the phenotype-ranking and sensitivity sections, leaving scores, methods and conclusions unchanged. The Track 2 report has its generative AI disclosure moved before the abstract and reformatted to match Track 1. The artifact manifest records both changes. Original reports mention development-workspace locations; those references are provenance, not paths required to execute this repository. See the root README for the portable implementation and training recipes.

## Final automatic scientific run

The completed `e2e_enriched_v4` scientific run (workflow 0.4.1) produced a recorded audit trail: the genetic handoff, discovered candidates, per-candidate assessments and gates, the combined final result, and the stage prompts, model generations, tool calls, source ledger and run manifest. Its recorded outcome is `review_required`.

Those artifacts are **not** included in this repository. They embed the proband clinical narrative, HPO terms and variant-level calls carried through every stage prompt, and that content remains subject to the dataset access conditions. They are retained privately and can be supplied to the organisers through an authorised channel on request.

The run's conclusions, gates and limitations are reported in the Track 2 methods report, which is included. The later workflow 0.5.0 added the automatic Track 1 CSV export; the Track 1 audit here comes from that export, rather than being represented as part of the earlier scientific run.

## Scope and integrity

These explicitly included submission artifacts contain case-derived findings and clinical context. The repository still excludes raw case inputs, the complete 52 MB variant ranking, model checkpoints, feature arrays and downloaded training databases. The historical source-code snapshot is also excluded to keep one maintained implementation; its recorded hashes remain in the privately retained run manifest.

The only content change to an included artifact is the generalisation of the clinical feature narrative in the Track 1 report, recorded in `artifact_manifest.json` with both the original and packaged hashes. No scientific claims, scores or review failures were rewritten.

For each artifact, compare SHA-256 of the packaged file with `packaged_sha256` in `artifact_manifest.json`. The original hashes allow comparison with the privately retained source artifacts. This folder is an inspection bundle, not a complete resumable checkpoint: the recorded run artifacts, the full ranking, original absolute paths and other private dependencies are deliberately absent.
