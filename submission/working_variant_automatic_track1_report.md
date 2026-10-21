# MVA Hackathon 2026, Track 1: automatic phenotype-guided variant prioritisation

| Field | Description |
|---|---|
| Team | Working variant |
| Proband | PROBAND01 |
| Analysis date | 16 October 2026 |
| Submission file | `working_variant_automatic_track1.csv` (10 automatically selected primary hypotheses) |
| Primary result | BUB1B c.2210T>G p.(Leu737Ter) paired with c.3006T>G p.(Asn1002Lys) |
| Scope | Local variant scoring, automatic phenotype–gene ranking, partner prioritisation and structured Track 2 handoff |

## Generative AI use

Evo 2 and ESM-1v provide locally computed sequence representations for a supervised variant classifier. Track 1 gene ranking, partner selection and submission EPCR assignment are deterministic Python calculations; no conversational language model selects the lead gene or pair. The connected Track 2 workflow uses a locally hosted Muse-Glimmer-30B model with local evidence tools. Patient inputs, model inference and case-run evidence retrieval remain on the local cluster; the automatic case workflow implements no external API retrieval.

AI coding assistance was used to develop the workflow and prepare this methods description. The automatic Track 1 inference path uses the inputs and calculations documented below. Expected answers are read only by a separate post-run evaluator.
Models used: Claude Opus 5.5 (Anthropic API), GPT-6 Astra (OpenAI API), commercial terms, no training on customer content.


## Method abstract (≤ 500 words)

We implemented a local workflow that connects genomic foundation-model features to automatic phenotype-guided gene ranking, second-allele prioritisation, submission generation and a structured mechanistic-research handoff.

A supervised variant classifier uses Evo 2 DNA likelihood changes and embeddings, ESM-1v protein scores and embeddings, and VEP annotations. It was trained on public ClinVar labels with gene-grouped train/validation/test splits; BUB1B, CEP57 and TRIP13 were held out from training. The final classifier averages LightGBM and MLP predictions. Population frequency and gene identity are excluded from its learned features. Held-out variant AUROC was 0.978 overall and 0.969 for missense variants.

For the proband, the workflow scored 40,629 prepared candidates using precomputed sequence features. Seven observed HPO terms were matched to local gene annotations through ontology ancestors and information content. A transparent heuristic combines the highest eligible variant score with phenotype similarity. Partner hypotheses are ordered by selected coding/splice consequences, absent annotated population frequency and pathogenicity score.

BUB1B ranked first. Its lead variant, c.2210T>G p.(Leu737Ter), scored 0.951; the first partner, c.3006T>G p.(Asn1002Lys), scored 0.075. The pipeline selected this pair without manual gene promotion or partner substitution during execution. The submission contains 10 automatically selected primary hypotheses. EPCR is calculated as the gene ranking score divided by the within-gene hypothesis rank; the leading pair receives 0.848188180. These are uncalibrated ranking estimates, not posterior causal probabilities.

Repeated scoring reproduced all variant scores exactly. Input-order permutation preserved the ranking, while removing the single most influential observed phenotype term or the lead variant's strong score changed the leading gene. The workflow passed 28 regression tests, including submission export checks.

Unknown phase and unresolved residual function are preserved in the Track 2 handoff. Local candidate discovery and evidence auditing can return review required when scientific claims remain unsupported. A prioritised variant pair does not establish compound heterozygosity or a treatment indication.

The core scoring and HPO-ranking components are reusable for other rare-disease cases with compatible variant and phenotype inputs. Cross-disease performance has not been established. General use requires disease-independent candidate preparation, appropriate inheritance interpretation and disease-relevant Track 2 evidence. Workflow development used knowledge of this case's outcome, so this case is not a blinded diagnostic evaluation. Raw-read processing, automated structural-variant follow-up and calibrated pair probabilities are outside the demonstrated automatic path. From prepared inputs, gene and partner ranking, EPCR assignment and competition-CSV export execute automatically with an inspectable audit trail.


## 1. Approach in detail

The workflow combines learned sequence-based pathogenicity scores with the proband's observed clinical features. It ranks genes, selects a lead variant within each gene, orders candidate partner alleles and creates a structured handoff for mechanistic research. These steps, including submission-row selection and CSV export, execute automatically from prepared, annotated candidate variants and manually transcribed observed HPO terms.

### Applicability to other cases and diseases

**The core Track 1 scoring and phenotype-ranking components can be applied to other rare-disease cases; their diagnostic performance beyond this case has not yet been established.** The classifier scores sequence variants without using gene identity as a learned feature, and the phenotype-ranking algorithm uses supplied HPO terms and public gene annotations rather than a fixed MVA diagnosis or a hard-coded BUB1B preference.

For another case, the required inputs are that patient's prepared, annotated candidate variants and observed HPO terms. The scoring path also requires compatible Evo 2 and ESM-1v features, either generated for that case or supplied from an existing cache. Alternatively, the ranking stage can consume a compatible scored-variant table. Reference assembly, transcript annotations and feature formats must match the trained model's requirements. The model need not be retrained simply to execute on a different patient, although performance in a new disease or variant distribution requires evaluation.

The complete demonstrated workflow is **not yet a validated general-purpose diagnostic pipeline**. Before using it outside MVA, candidate preparation should remove the MVA-gene exceptions and apply an appropriate disease-independent policy. The partner-ordering rule is a simple heuristic, not a general model of dominant, recessive, X-linked or de novo disease; inheritance interpretation requires further development and review. HPO annotation coverage also limits phenotype matching. Raw-read processing, phase determination and structural-variant follow-up are not automated by this workflow.

The Track 1-to-Track 2 interface is reusable, but the current Track 2 evidence collection and some tools are MVA-focused. A different disease requires suitable evidence coverage and evaluation of the downstream mechanism and safety assessments. Thus, the implementation provides a reusable research framework for other cases, while the demonstrated result remains an analysis of one case rather than evidence of diagnostic or therapeutic generalisation.

### Training data and sequence models

1. **Training labels.** ClinVar (GRCh38, 2026-10-04) germline variants with at least 1 review star:
   - Pathogenic/Likely pathogenic vs Benign/Likely benign, conflicting and VUS removed.
   - SNVs and indels up to 50 bp on chr1–22, X and Y.
   - 1,704,763 variants in total.
2. **Gene-grouped split, 80/10/10.** Genes that share variants were merged into connected components, and whole components were assigned to train, validation or test. **BUB1B, CEP57 and TRIP13 were forced into the test split**, so a ClinVar entry for the proband's causal variant could not leak into training.
3. **Subsampling for compute.** Within each split and consequence group we kept all P/LP and up to 2× as many B/LB variants. Truncating and canonical-splice P/LP were capped at about 20k each. This gave 271,170 variants, with a `sample_weight` to restore the ClinVar mix.
4. **Annotation.** Ensembl VEP 115 (offline cache, `--pick` MANE-first, `--distance 0`). VEP consequence, not ClinVar's, is used as a feature, because that is all that is available for a new genome.
5. **Features.**
   - **Evo 2 7B** (Apache-2.0) on 4 kb reference and alternate windows centred on the variant:
     - change in log-likelihood (alternate minus reference);
     - `blocks.28.mlp.l3` embeddings at the variant and averaged over the window (alternate minus reference, plus reference).

     We benchmarked 4 kb against 8 kb windows: 4 kb was twice as fast with the same validation AUROC. Vortex Triton kernels were rejected because they changed outputs.
   - **ESM-1v** (MIT; 5-model ensemble) for missense variants:
     - masked-marginal score from each model;
     - model-1 embeddings at the mutated site and averaged over the protein (mutant minus wild type, plus wild type).

     Proteins longer than 1,022 amino acids were windowed. Wild-type residues matched Ensembl peptides for 100% of the 173,342 variants scored.
   - VEP consequence terms, impact, variant type and indel length.
   - **Deliberately excluded:** gnomAD allele frequency, because ClinVar benign labels are often frequency-based and the model would learn that labelling rule; and gene identity.
6. **Models.** All selection was done on validation only.
   - Logistic regression and LightGBM on the scalar features plus PCA of the embedding blocks.
   - An MLP on the full embeddings (3 seeds).
   - **Final model: the mean of the LightGBM and MLP probabilities.**

The trained ensemble is frozen for this analysis. No retraining was performed on the patient or on the competition answer.

### Candidate preparation and automatic scoring

The prepared input contains 40,629 variants: PASS calls inside transcripts with gnomAD maximum AF below 1% or absent, together with additional PASS calls in BUB1B, CEP57 and TRIP13. This disease-specific preparation is an explicit limitation of the input set; the downstream ranking has no special boost for those genes.

Evo 2 and ESM-1v features were computed in a separate feature-extraction stage and stored locally. The automatic scoring stage applies the trained LightGBM/MLP ensemble to those features. Repeated scoring of all 40,629 variants produced identical values (maximum absolute difference 0.0). Validation-derived tiers are: high at approximately 0.946 (99% specificity on benign validation variants) and moderate at approximately 0.600 (95% specificity).

Before gene ranking, variants with depth below 10, genotype quality below 20, or known maximum AF of at least 1% are excluded. Missing depth or genotype quality is flagged rather than silently treated as measured quality. Allele frequency is used for eligibility and partner ordering, never as a pathogenicity-classifier feature.

### Phenotype–gene ranking

Seven observed proband features are supplied as HPO terms. One further supplied term is retained as family history and is not scored as a proband phenotype. The individual terms are not reproduced here; they remain in the access-controlled case inputs. The ranking inputs contain observed clinical features and candidate variants, without an inferred MVA diagnosis or a designated causal gene or pair.

Local HPO gene annotations are expanded through ontology ancestors. For each observed term and candidate gene, the algorithm finds the shared ancestor with the highest information content, normalises by the query term's information content and averages over the observed terms. Information content is computed from local gene annotation frequencies. Inheritance and clinical-modifier terms are excluded from the phenotype similarity calculation.

The gene ranking is:

`gene score = 0.55 × highest eligible variant score + 0.45 × phenotype similarity`

This is a transparent development heuristic, not a calibrated causal probability. Inheritance annotations and genotype completeness are reported separately. The ranking weights have not been validated on an independent patient cohort.

### Partner hypotheses and Track 2 handoff

Within each gene, the highest-scoring eligible variant is the lead. Remaining eligible variants are ordered by:

1. Missense, stop-gained, frameshift or canonical splice donor/acceptor consequence before other consequences.
2. Absence of an annotated population frequency before a recorded frequency.
3. Original pathogenicity score, with a deterministic variant-key tie-break.

This prioritises plausible second-allele hypotheses without increasing their pathogenicity scores. It is not a trained compound-heterozygosity model. Nearby damaging calls within 20 bp are flagged as potentially representing one event, rather than counted as confirmed independent alleles.

The machine-generated handoff contains the lead gene, variants, original scores, clinical features, alternative genes, inheritance annotations and requested follow-up. Phase and residual function are explicitly carried as unresolved. Track 2 can therefore investigate a conditional mechanism without converting an uncertain partner into established biallelic disease or assuming a rescuable protein product.

### Automatic submission rows and EPCR assignment

The workflow generates a CSV using the competition's exact column schema, with at most 10 rows; the submission file is named `working_variant_automatic_track1.csv`. The exporter derives all rows and EPCR values deterministically from the computed gene ranking and partner hypotheses.

For each gene, an autosomal-recessive HPO annotation and a heterozygous lead permit candidate pairs with eligible heterozygous partners in the partner order defined above. Partners must be on the same chromosome, more than 20 bp apart, and the gene must not have a clustered-call flag. These checks do not establish trans phase or resolve mixed inheritance annotations. The within-gene hypothesis order is the first eligible pair, the single lead, then remaining eligible pairs; without an eligible pair, only the single lead is proposed.

`EPCR = gene ranking score / within-gene hypothesis rank (1-based)`

All generated hypotheses are sorted globally by descending EPCR, identical singles or unordered pairs are deduplicated, and the first 10 are exported. Ties are broken deterministically. All rows are primary candidate explanations; secondary or incidental disease-risk classification is not inferred. The automatic shortlist can retain hypotheses in complex regions, including HLA genes, and other candidates needing artefact review; CSV export does not replace the unimplemented raw-read review.

EPCR values are **uncalibrated ranking estimates**, consistent with the competition's relative-ranking use of this field. The reciprocal-rank discount is an explicit engineering rule, not a learned causal probability or a measurement of partner function. The gene score and local hypothesis rank determine every exported value. Original pathogenicity scores remain unchanged and appear in the row notes. The scores need not sum to one because hypotheses can overlap.

The leading BUB1B pair receives EPCR **0.848188180**, calculated as 0.84818818 / 1. This is a ranking estimate, not an independently validated 84.8% probability of causality. `02_submission.json` records the calculation for every selected row, the selection policy and the number of omitted hypotheses; a submission-side copy is provided as `working_variant_automatic_track1_epcr_audit.json`. The workflow also writes checksums for the audit and CSV.

## 2. Results supporting the submission

### Held-out variant-classifier performance

The following results measure variant classification on held-out genes. They do not measure patient-level diagnostic accuracy.

**Model performance on held-out test genes** (30,088 variants; AUROC):

| Model | Overall | Missense | Intronic | Splice | Synonymous |
|---|---|---|---|---|---|
| Evo 2 alone, no training (Δ log-likelihood) | 0.913 | 0.886 | 0.893 | 0.832 | 0.858 |
| ESM-1v alone, no training (missense only) | – | 0.881 | – | – | – |
| Trained, scalar features only | 0.964 | 0.944 | 0.907 | 0.877 | 0.886 |
| **Trained, final ensemble** | **0.978** | **0.969** | **0.954** | **0.918** | 0.902 |

Test results matched validation (0.978 vs 0.978 overall). The embeddings add about +0.025 AUROC on missense, +0.05 on intronic and +0.04 on splice variants over the scalar-only model.


### Automatic case result

BUB1B ranked first with a gene score of **0.84818818**, comprising a lead variant score of **0.95075093** and phenotype similarity of **0.72283371**. ATP6V1E1 ranked second (0.81707125).

| Role | Variant, GRCh38 | Coding/protein change | Ensemble score | Genotype |
|---|---|---|---|---|
| Lead | chr15:40209701 T>G | BUB1B c.2210T>G, p.(Leu737Ter) | 0.950751 | Heterozygous |
| First partner hypothesis | chr15:40220612 T>G | BUB1B c.3006T>G, p.(Asn1002Lys) | 0.074798 | Heterozygous |

Transcript: ENST00000287598.11. The stop-gained variant was ranked 25th by variant pathogenicity alone. Adding phenotype compatibility moves its gene to first place. The missense partner remains low-scoring; its prioritisation arises from the explicit partner-ordering rule, not a hidden score adjustment. Its frequency is absent from the supplied annotation, which does not itself prove pathogenicity.

The pipeline exports this BUB1B pair as the **first row of `working_variant_automatic_track1.csv`**, with EPCR **0.848188180**. The remaining rows are alternative primary hypotheses selected by the same automatic rule. No competition score is claimed for this submission.

The automatic result retains `one_damaging_allele_second_unresolved`, `phase: unknown` and `residual_function: unestablished`. The ranking does not establish trans phase or the biological effect of the missense variant.

### Reproducibility and sensitivity

- Repeated ensemble scoring of the same feature inputs produced identical scores for all 40,629 variants.
- Permuting input variant order did not change the ranking.
- Removing either the single most influential observed phenotype term or the lead variant's strong pathogenicity score changed the top gene to ATP6V1E1.
- Output checks confirmed that the ranked BUB1B lead and first partner were carried consistently into the handoff and submission, while phase remained explicitly unknown.
- Repeated execution reproduced the complete gene ranking and EPCR audit. The export contains 10 rows matching the official CSV schema; CSV and checkpoint hashes passed verification.
- The workflow including submission export passed 28 regression tests. Checkpoint and source checksums passed; resuming the completed run reused results and reproduced its report without further model inference.

These are reproducibility and sensitivity checks, not independent diagnostic validation. Workflow development used knowledge of this case's outcome. Keeping expected answers outside inference prevents their direct use as inputs, but does not make the evaluation blinded or eliminate case-specific design bias.

## 3. Automated output or manual curation

Variant scoring, HPO matching, gene ranking, partner ordering, submission-row selection, EPCR assignment, CSV export and structured handoff are automatic. The case's clinical features were manually transcribed into HPO inputs; model training, annotation and feature generation were completed beforehand. The automatic scoring stage consumes those feature caches.

The submission contains 10 primary hypotheses selected from the automatic results. The exporter considers alternative genes, calculates each EPCR using the stated formula, and writes the final CSV without manual row selection or score assignment. No calibrated pair probability has been estimated.

## 4. Downstream manual review and scope

The automatic ranking requires no manual gene promotion or partner substitution. Review is still needed to interpret unresolved phase, possible hypomorphic function and alternative diagnoses. The workflow requests depth, structural-variant and mobile-element follow-up; it does not execute these analyses automatically.

The submitted ranking does not use raw-read realignment, exon-depth analysis, structural-variant calling or mobile-element detection to exclude alternative second alleles. Such analyses would be appropriate follow-up where indicated. Parental testing could resolve phase; RNA studies or targeted long-read analysis could investigate functional effects and variants missed by the prepared small-variant input.

## 5. Public vs proprietary data

Training and reference resources are public. The proband genome and clinical description are restricted hackathon-provided data, used under the applicable hackathon terms. No proprietary external training dataset was added.

## 6. Public data and tools used

- ClinVar GRCh38, 4 October 2026: supervised training labels; gene-grouped held-out evaluation.
- Ensembl release 115: reference assembly, peptides and offline VEP annotations, including bundled population-frequency annotations used outside the classifier.
- Human Phenotype Ontology: local ontology and gene-to-phenotype annotations for automatic clinical matching.
- Evo 2 7B and ESM-1v: local DNA and protein features.
- LightGBM, PyTorch and scikit-learn: trained classifiers and feature transformations.
- Python and Slurm: deterministic ranking, workflow orchestration, checkpoints and provenance.

The classifier excludes gene identity and population frequency as learned features. BUB1B, CEP57 and TRIP13 were assigned to held-out test groups during training. The HPO ranking intentionally uses public gene–phenotype knowledge; that is a distinct information source from variant training labels.

## 7. Proprietary data

No proprietary training data. Patient-derived inputs and outputs remain subject to the hackathon's data-use and deletion requirements.

## 8. Compound heterozygous output

The first pair is a **candidate compound-heterozygous hypothesis** comprising p.(Leu737Ter) and p.(Asn1002Lys). The workflow does not infer trans phase from the presence of two heterozygous variants. The first partner's low pathogenicity score is preserved, and functional rescue remains conditional on further evidence.

## 9. Secondary and incidental findings

Alternative genes compete as primary candidate explanations according to the same automatic rule. The exporter does not perform a dedicated incidental-finding classification, so it emits no secondary rows. An exported alternative-gene hypothesis, including LZTR1, represents a candidate explanation under the ranking rule; it is not a separate assessment of incidental disease risk.

## 10. Execution, provenance and connection to Track 2

A single command submits ensemble scoring and dependent downstream analysis through Slurm. The dependent workflow generates and seals the submission CSV before optional Track 2 reasoning; `--track1-only` also produces the CSV. Models and case data remain local. Checkpoints, input hashes, code snapshots and structured outputs support inspection and controlled resumption. The executed path begins with prepared annotated variants and precomputed features, not raw FASTQs; complete feature regeneration through the adapter has not been exercised for this case.

The connected reasoning job completed in 9 minutes 48 seconds with 60 logged tool calls. This is the downstream job's wall time, not total raw-data processing time or a Track 1-only benchmark. It excludes input annotation, feature extraction and model training. No monetary cost was estimated.

Track 2 receives the case directly from Track 1, performs unseeded candidate discovery and checks compound identity, source-linked mechanism, counter-evidence and safety context with local tools. A separate model context reviews the claims, followed by at most one automatic correction and re-review. Nicotinamide mononucleotide was proposed during candidate discovery, but did not pass the final claim/evidence audit; the final outcome remained **review required**, not a clinical nomination. This demonstrates an executable handoff while preserving the distinction between prioritising a variant hypothesis and validating a treatment hypothesis.

The audit trail can be accessed through the github repository linked to this submission, at https://github.com/mic-pik/rare_disease_real_kid_hackathon_2026.