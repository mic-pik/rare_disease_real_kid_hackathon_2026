# BUBR1 rescue evidence: what can actually be transferred to a drug?

Audit date: 2026-10-19. Bounded primary-evidence review by the coding assistant; not a systematic review or independent expert assessment. No private patient data or patient-specific drug recommendation.

## Decision

**There are experimentally supported ways to improve BUBR1-related defects, but this audit does not establish an approved drug that reproduces the relevant rescue at a defensible exposure.** The strongest evidence supports genetic positive controls and identifying the functional defect before further drug ranking. No drug is promoted and historical candidate decisions remain unchanged.

The machine-readable ledger is `config/bubr1_rescue_audit.json`. Ten experiment-level records separate the model, intervention, endpoint, result and transfer limitation. Generated mapping, exposure calculations, source snapshots and manifests are under `runs/bubr1_rescue_audit_20261019/`.

## Findings that change the next step

**Restoring protein can work when useful function remains.** Patient-cell BUBR1 complementation restored checkpoint activity. In engineered replacement cells, raising I909T or L1012P expression restored the nocodazole response. These results support testing abundance and residual function; they do not establish that increasing any mutant protein is beneficial. The same study relates segregation errors to protein depletion, but it does not supply a universal therapeutic protein threshold. [Suijkerbuijk et al., 2010, especially Figure 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC2887387/).

**Spatial PP2A-B56 recruitment is a concrete rescue mechanism.** A MIS12–KARD-3D construct improved chromosome alignment in two MVA patient-derived fibroblast lines. In contrast, the ZM447439 chemical-rescue experiments were performed in engineered HeLa/B56 models. They used 2 µM drug for one hour in MG132-arrested cells. The patient-cell localization experiment and chemical experiment must not be combined into a claim that an Aurora inhibitor rescued patient-cell segregation. [Xu et al., 2013, Figures 3–4](https://pmc.ncbi.nlm.nih.gov/articles/PMC3654266/).

**The endpoint matters.** Later work found defective end-on attachments and Astrin recruitment despite congressed chromosomes in a BubR1 PP2A-docking mutant. Alignment alone can therefore miss an important remaining defect. [Shrestha et al., 2017, Figure 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC5532248/).

**More phosphatase activity is not automatically better.** Engineered phosphatase-dominant conditions weakened checkpoint signalling; excess recruitment in another construct impaired error correction and increased segregation errors. Restoring the right localization and balance is different from globally activating PP2A. [Kinase–phosphatase balance study, 2023, Figures 2 and 7](https://pmc.ncbi.nlm.nih.gov/articles/PMC10577578/).

## Endpoint-specific evidence ledger

| Intervention | Model and measured endpoint | Interpretation |
|---|---|---|
| Wild-type BUBR1 expression | Patient-derived lines; checkpoint restoration | Genetic positive control; not a drug |
| Increased expression of selected unstable mutants | Engineered replacement cells; checkpoint response | Residual-function-dependent abundance route |
| Kinetochore-targeted PP2A-B56 recruitment | Patient-derived fibroblasts; alignment | Spatial mechanism supported; sustained fidelity not established |
| ZM447439 | Engineered HeLa/B56 models; congression/K-fibers | Partial chemical endpoint rescue |
| Hesperadin or ZM447439 | B56-depleted RPE1; cold-stable K-fibers | Related lesion, not direct patient-model rescue |
| BI2536 | B56-depleted cells | Conflicting results across studies; no general PLK1 rescue rule |
| NMN | Mouse tissues and MEF experiments; BubR1 abundance | Abundance evidence, not patient-cell segregation rescue |
| SIRT2 overexpression | BubR1-hypomorphic mice; lifespan | Genetic organismal effect; cannot be assigned to NMN |
| Congression without normal attachment maturation | Engineered BubR1 docking mutant | Counterexample to alignment as a sufficient endpoint |
| Excess phosphatase/recruitment | Engineered mitotic balance | Counterexample to monotonic “more PP2A is better” reasoning |

In the related RPE1 study, Hesperadin 50 nM or ZM447439 1 µM restored cold-stable fibers after short exposure in mitotically arrested cells. BI2536 also restored fibers there, whereas the Xu study reported failure and monopolar spindles under its conditions. That conflict prevents a context-free PLK1-inhibition hypothesis. [Foley et al., 2011, Figures 4–5](https://pmc.ncbi.nlm.nih.gov/articles/PMC3186838/).

The NMN experiment used 1 mM for 24 hours in MEFs. Its abundance result must remain separate from the SIRT2-transgenic lifespan result. Neither establishes an approved pharmacological SIRT2 activator or segregation rescue in the relevant patient model. [North et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4194088/).

## Mapping to approved compounds and exposure

The frozen FDA/ChEMBL branch contains 164 AURKB records across 69 mapped ingredients. Its existing quality rules retain 21 active records across 14 ingredients. Binding, inhibition and conflicting measurements remain separate. No SIRT2 record met the existing active/unflagged threshold; this is a bounded database result, not proof that no relevant compound exists.

Four AURKB candidates received source-level exposure checks: axitinib for strong binding, crizotinib for binding plus explicit inhibition, and pexidartinib/momelotinib for functional inhibition records. The other ten active-signal ingredients have **not** received equivalent exposure reviews and are not ruled out wholesale.

| Ingredient | Local AURKB evidence | Calculated unbound plasma peak reference | Assessment |
|---|---|---:|---|
| Axitinib | Kd 11 and 236 nM | <0.72 nM | Even the lowest binding value exceeds this reference; no cellular rescue evidence |
| Crizotinib | Kd 76 and 1,240 nM; >80% inhibition at 500 nM | ~124 nM, specifically the label's pediatric oncology PK cohort | Binding overlap is insufficient; cellular engagement and rescue unknown |
| Pexidartinib | IC50 643 nM; separate kinobead Kd >30,000 nM | <206 nM | Exposure concern and discordant assays; no rescue window established |
| Momelotinib | <50% inhibition at 100 nM, >50% at 1,000 nM; Kd >30,000 nM in another assay | ~104 nM | Censored potency and assay disagreement do not establish effective cellular exposure |

These calculations use `total Cmax (ng/mL) × 1000 / parent molecular weight × unbound fraction`. Parent molecular weights are calculated from the frozen catalogue's structures. Axitinib and pexidartinib binding is reported as >99%, so the values are **upper bounds**, not point estimates. Concentrations are label population summaries; they are not maxima across patients, intracellular measurements, sustained engagement, rescue doses or safety margins. No direct comparison is made between a substitute drug's plasma concentration and ZM447439's culture concentration.

PK sources: [INLYTA](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=84137882-e000-47da-bd5b-fa76ab3c76f9), [XALKORI](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=2a51b0de-47d6-455e-a94c-d2c737b04ff7), [TURALIO](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=7b0ec4bf-b2ac-49ee-96af-7740ad3f2e16), [OJJAARA](https://dailymed.nlm.nih.gov/dailymed/getFile.cfm?setid=4a672fd1-eb7b-4aa9-9b23-845e4b5dc400). TURALIO also carries a serious liver-injury boxed warning. A full suitability review of all mapped ingredients was outside this audit.

Original assay values, relation operators, descriptions, quality flags and document identifiers are retained in `analysis/aurkb_assays_all.parquet`. The source papers behind all kinase measurements were not independently inspected; those numerical claims come from the original local ChEMBL records. The strongest assays must not erase weaker or conflicting measurements.

## Transfers explicitly rejected or deferred

- **Fingolimod → PP2A-B56 recruitment:** PP2A modulation has been reported in hematologic models, with distinct parent/phosphorylated-species effects. This does not establish kinetochore recruitment or BUB1B rescue. No applicable cellular-rescue exposure margin was established here. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3772499/).
- **MG132 → approved proteasome inhibitor:** MG132 was used to arrest cells for measurement. Its inclusion does not identify it as the rescuing treatment or justify substituting an approved proteasome inhibitor.
- **Research inhibitor → approved inhibitor with a similar target:** requires compound-specific cellular engagement, rescue, exposure and growth data. Shared binding alone is insufficient.
- **Abundance increase → correction of the case:** requires residual functional protein; no allele-specific functional conclusion is made in this audit.

## Concrete next dependency

Before another drug shortlist, establish whether the relevant model is rescuable and which defect dominates:

1. Quantify residual BubR1 abundance/localization and baseline checkpoint, attachment and segregation defects in deficient versus corrected cells.
2. Use wild-type BUBR1 restoration as a positive control. If localization is implicated, assess a targeted PP2A-B56 restoration control. These are mechanistic experiments, not proposed patient therapies.
3. Measure segregation errors per completed division, mature attachments, absolute daughter-cell survival and proliferation. Keep arrested-cell alignment as a secondary endpoint.
4. Only nominate a drug after it reproduces relevant rescue without merely slowing division, killing abnormal cells or bypassing the checkpoint, at an independently justified exposure.

No cell experiment was performed or scheduled. Without access to a relevant model or existing functional measurements, additional drug nominations remain deferred. The audit supports a narrower experimental question, not another unrestricted compound-ranking run.

## Provenance and limits

Four primary full-text XML snapshots and the existing North full-text cache were archived with hashes. The 2010 paper's abstract and indexed Results were inspected, but fullTextXML retrieval failed twice with HTTP 500; supplements were not inspected. Label and fingolimod sources were browser-inspected, with exact URLs recorded rather than claimed as archived snapshots.

The reproducible local analysis script is `runs/bubr1_rescue_audit_20261019/audit_local_mapping.py`. It verifies upstream manifests, exports all AURKB rows, checks exact candidate identities and computes exposure references. It does not score drug efficacy. `analysis/manifest.json` records inputs and output hashes. Research-control absence from the mapped FDA catalogue is explicitly not treated as a comprehensive regulatory determination.
