# Follow-up of the fixed 20 LINCS hypotheses

Completed 2026-10-19. **No candidate gains sufficient evidence for a replacement Track 2 nomination.** We reviewed the remaining sixteen ingredients, removed annotated cell-cycle/replication genes, examined individual BUB1B shRNAs, and added contextual PRISM viability measurements. The original candidate list, scores and thresholds remain unchanged.

The most consequential finding is that the BUB1B knockdown signatures do not provide a consistent independent representation of the patient expression signature. This weakens the interpretation of expression reversal as BUB1B rescue. It is not evidence that all twenty drugs are ineffective.

## Exposure-first review of the remaining sixteen

Concentrations below are the actual retained LINCS conditions, not proposed doses. Total plasma measurements are study-specific references, not maximum achievable exposure or intracellular/free-culture measurements. Parent molecular weights were calculated from the frozen catalogue with RDKit; `ng/mL ÷ g/mol = µM`. We do not claim a quantified exposure margin where the available review established only binding or a molecular-species mismatch. No dose escalation is proposed.

| Ingredient | LINCS µM | Source-based exposure/mechanism assessment | Current interpretation |
|---|---:|---|---|
| Pentamidine | 10 | Injection label reports mean total Cmax 612 ng/mL after a 4 mg/kg two-hour infusion in six AIDS patients, approximately 1.80 µM parent equivalent. [Label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=d46aea4f-8f3f-4ce0-93d6-94e147f5a31a) | Exposure gap against this reference; no rescue window established. |
| Midostaurin | 10 | Parent and two active metabolites are >99.8% protein bound; label describes inhibition of proliferation and induction of apoptosis in relevant malignant cells. [Label](https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=11fa3fc9-6776-49a6-b1c1-653f627c3e58) | Parent/metabolite contributions and free exposure unresolved; antiproliferative alternative remains. No quantitative Cmax margin established here. |
| Disulfiram | 10 | Human PK work measures parent plus multiple metabolites. CuET, a copper-containing derivative, has a demonstrated NPL4/protein-homeostasis mechanism in cancer models. [Human study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6379104/), [primary mechanism study](https://www.nature.com/articles/nature25016) | Cannot transfer metabolite/CuET evidence to a nominal-parent rescue claim; copper/species and exposure unresolved. |
| Neratinib | 10 | >99% protein binding; active metabolites and food/interactions affect exposure. [Label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=d1e41dcf-e82a-47c2-a0ad-6c6eef621834) | No verified 10 µM cellular exposure bridge or segregation-rescue mechanism. Binding alone does not prove unreachability. |
| Mebendazole | 10, 40 | Label mean Cmax 56.2 ng/mL in its fed 500-mg study, approximately 0.190 µM; 90–95% binding. Primary cancer-cell work demonstrates abnormal spindles, tubulin depolymerization and mitotic arrest. [Label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=a720c424-517a-461b-8e13-6c37caf8aa3f), [primary study](https://pubmed.ncbi.nlm.nih.gov/12479701/) | Exposure concern plus mechanism directly conflicting with the desired fidelity endpoint. Do not promote. |
| Erlotinib | 0.4, 10 | 93% protein binding. [Label](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=3ecafc1e-979a-4aff-a276-18fc0b110491&version=1) | The only retained 0.4 µM profile has a negative reversal score; high-dose aggregate does not establish low-dose benefit. Absolute PK margin not independently established here. |
| Afatinib | 10 | Approximately 95% binding and EGFR-directed oncology pharmacology. [Label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=ec5b9e00-3eab-4dbd-9bb3-75865b1af485) | No quantified exposure bridge or BUB1B benefit established; growth-related explanation unresolved. |
| Dasatinib anhydrous | 10 | Label steady-state Cmax at 100 mg daily is 82.2 ng/mL, approximately 0.168 µM; binding about 96%. [Label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=76e1833c-e105-5a3d-0f01-1330d8cc7ae1) | Large nominal exposure gap against this reference; PRISM also flags growth loss in several mapped lines. |
| Desipramine | 10 | Primary PK study reports observed peaks of 8–46 ng/mL in its studied single-dose settings, up to approximately 0.173 µM. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2950987/) | Dose- and population-limited exposure concern, not a ceiling across all treatment settings; no benefit mechanism established. |
| Leflunomide | 10 | Parent rarely detectable in plasma; active metabolite is teriflunomide. [Label](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=acf3efe4-0d70-4718-9c50-789d69ab215d) | Parent/metabolite mismatch; do not substitute metabolite exposure for the LINCS parent. |
| Paroxetine | 10 | Label says plasma concentrations normally <400 ng/mL, approximately <1.21 µM, with 93–95% binding in specified assays. [Label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=f39a2a7a-7a7e-4dd8-9704-54d9598504fb) | Exposure concern; no functional-rescue bridge. This is not an assertion that concentrations above the reference never occur. |
| Glipizide | 10 | Mean Cmax 523 ng/mL, approximately 1.17 µM, after the studied immediate-release 5-mg dose; label binding 98–99% and hypoglycemia concern. [Primary PK study](https://pubmed.ncbi.nlm.nih.gov/16509763/), [label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=a59e2aaf-692c-4696-817b-3fed7a9dc52c) | No BUB1B benefit or exposure bridge to balance glucose-lowering effects. |
| Sunitinib | 10 | Label reports combined parent/metabolite plasma concentrations of 63–101 ng/mL by day 14; parent and metabolite binding 95%/90%. [Label](https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=a5d555f5-d66a-4f94-abcf-96fa6d71a32f&type=display) | Combined-species values must not be converted into parent Cmax. No 10 µM parent exposure bridge established. |
| Podofilox | 2.5, 10 | Topical label describes an antimitotic drug and peak serum levels 1–17 ng/mL after the specified topical applications, up to approximately 0.041 µM. [Label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=6af07697-5203-43e0-ba35-83683af851c3) | Route/exposure gap and antimitotic mechanism; PRISM growth-loss signal reinforces concern. Do not promote. |
| Nilotinib | 0.08, 0.4, 2, 10 | Label mean Cmax 1,540 ng/mL at 300 mg twice daily in newly diagnosed adult CML, approximately 2.91 µM total; serum binding ~98%. [Label](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=b34e7fd4-53c5-46d8-bf8a-1fcbe7a16504) | Low-dose profiles merit inspection, but direction is mixed across cell models and no independent BUB1B support is obtained. Total peak overlap alone is insufficient. |
| Sorafenib | 10 | Label protein binding 99.5%; primary work confirms extensive binding in cancer-patient plasma. [Label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=93f9c3b2-d4c5-80d3-5421-8b4c29fbf8a6), [primary binding study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3652000/) | No measured free-culture/exposure bridge established; a positive shRNA reversal in one context is insufficient, particularly when that context does not reproduce the patient signature. |

This is bounded triage, not a complete clinical pharmacology review of all sixteen. Missing quantitative exposure estimates are left unresolved rather than converted to exclusions. None obtains a demonstrated BUB1B-rescue mechanism from the reviewed evidence.

### Lower-dose profiles

These are individual measurements, not a matched dose-response experiment:

| Ingredient | Cell | Time | µM | Original reversal score |
|---|---|---:|---:|---:|
| Erlotinib | BT20 | 24 h | 0.4 | −0.1765 |
| Nilotinib | BT20 | 6 h | 0.4 | −0.0995 |
| Nilotinib | HS578T | 24 h | 0.08 | +0.1341 |
| Nilotinib | MDAMB231 | 24 h | 2 | −0.0463 |
| Nilotinib | SKBR3 | 6 h | 2 | +0.1489 |

## Cell-cycle exclusion

Before calculating the follow-up results, we specified removal of all genes annotated to local Reactome top-level **Cell Cycle** or **DNA Replication**. We did not select exclusions based on which compounds benefited. The exclusions remove ten of the 74 query genes: CCNB1, CHEK1, POLE2, CDC20, CDC45L, RFC5, CCNE2, AKT1, NUP62 and BUB1B. Historical annotation names are retained.

The remaining query has 50 increased and 14 decreased genes. We preserve the original rank scale and rebalance the two directions, then use 10,000 count/direction-matched random queries from 715 nonexcluded genes. BH adjustment still covers all 544 compounds, not just the selected twenty. Nine retain q ≤ 0.10 and positive scores in at least 75% of at least two cells: chlorhexidine, nisoldipine, midostaurin, disulfiram, neratinib, mebendazole, erlotinib, desipramine and nilotinib.

All twenty retain positive aggregate scores. Therefore, it would be incorrect to say eleven signals disappeared: they fail the follow-up statistical/consistency threshold. Several q-values are close to 0.10 and subject to Monte Carlo uncertainty. Removing BUB1B and other genuine mitotic biology also makes this a sensitivity test, not proof that the original signal was a growth artifact.

## Individual BUB1B knockdown check

Use only individual `trt_sh` records, with exemplar=1, replicate correlation ≥0.2 and at least three constituent samples. Consensus signatures and overexpression are excluded. Shared constituent experiments are deduplicated, repeated profiles of the same reagent are averaged, and agreement is evaluated among distinct reagents within cell/time.

Thirty-two qualifying profiles cover nine cell contexts, all at 96 hours. On the 715 retained measured genes, a context requires at least two reagents and median between-reagent centered-rank correlation >0.1. Only ASC (three reagents, correlation 0.151) and HT29 (five, 0.169) pass. Both have negative patient-direction agreement scores (−0.0547 and −0.0503). Most other contexts have near-zero reagent agreement; HCC515 has only one qualifying reagent.

For the two usable contexts, compare each cell-matched drug profile against each reagent using negative centered-rank cosine, then aggregate conditions and cells without choosing the best profile. No original candidate meets the exploratory multi-cell/reagent reversal rule. Of the twenty, thirteen have any usable context coverage and only six have two covered contexts. Missing coverage is not negative efficacy evidence. Sorafenib reverses the measured knockdown in its one covered context, but lacks the required replication and fails the noncycle q threshold.

This check does not establish a trustworthy independent disease signature: the same LINCS platform is used, off-target effects and achieved knockdown remain unresolved, drug exposures are 6/24 h versus 96 h for shRNA, and cell models differ from the patient. The mismatch is a reason to stop interpreting the original ranking as BUB1B-specific rescue, not to label the compounds universally ineffective.

## PRISM contextual viability

The local PRISM matrix maps to eighteen of the twenty ingredients in overlapping cell lines. Mapping uses exact Broad compound-root IDs and unique normalized CCLE cell-name matches, retaining full identifiers. The matching metadata describe 2.5 µM profiles; almost all LINCS conditions differ in dose, and the aggregated PRISM data do not supply a matched LINCS-time endpoint. They are not treated as patient-cell or normal-tissue safety measurements.

Dasatinib and podofilox each have relative viability below one-half (log2 value <−1) in five of six mapped lines. Median log2 values are −1.851 and −2.911, respectively. Mebendazole has two of five below one-half. The complete measurements, including compounds without such reductions, are retained. Absence of a reduction at a different dose/time is not evidence of safety, and these values do not distinguish cytostasis from cell death. Sorafenib and sunitinib have no retained mapped observations in this contextual join.

## Decision and reproducibility

No new nomination is justified. Further reranking of this single-donor signature is unlikely to resolve the central limitation: we have not established which observed transcriptional changes are consequences of BUB1B dysfunction that should be reversed. The completed branch provides hypotheses and counter-evidence, not a validated disease model. The original NMN proposal and its unresolved eligibility are unchanged.

The pre-result specification is `runs/lincs_followup_20261019/plan.json`; generated results are in `runs/lincs_followup_20261019/analysis/`. Important files include `fixed20_checks.csv`, `all_compound_checks.csv`, `query_ablation.csv`, `kd_inventory.csv`, `kd_selected.csv`, `kd_context_agreement.csv`, `drug_vs_individual_shrna.csv`, `kd_condition_scores.csv`, `prism_context.csv` and `prism_summary.csv`. Input and output hashes are recorded in `manifest.json`; the selected shRNA matrix is archived. The original baseline hashes are checked before analysis.

```bash
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 python pipeline/check_lincs_followup.py \
  --baseline runs/gse22206_lincs_20261019/analysis \
  --query runs/gse22206_20261019/analysis \
  --lincs ../gene_mol/data/lincs_raw --kg ../mva_llm/data/kg \
  --prism ../gene_mol/data/depmap_raw \
  --output runs/lincs_followup_20261019/analysis
python -m unittest discover -s tests -p 'test_lincs*.py' -v
```

Seven tests pass, covering scoring direction, multiple testing, condition weighting, conflicting reagents, constant signatures and matched-background exclusion. Sources for exposure review were browser-inspected; supplements and every original pharmacology study were not exhaustively reviewed. No external messages, submissions or laboratory experiments were sent or performed.
