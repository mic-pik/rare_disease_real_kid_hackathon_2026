# Working variant Track 2 report on NMN for preclinical investigation in MVA

**Rare Disease, Real Kid — MVA Hackathon 2026**  
**Team:** Working variant  
**Candidate:** Nicotinamide mononucleotide (NMN; local identity record PubChem CID 14180)  
**Analysis dates:** 12–16 October 2026  
**Report date:** 19 October 2026

We nominate **NMN for preclinical investigation of residual BubR1 rescue**. Our proposal asks whether increasing the abundance of remaining functional BubR1 can reduce chromosome-segregation errors while preserving cell viability and growth. The final automatic run discovered NMN without a supplied candidate list, assembled its supporting evidence and identified scientific uncertainties and claim-level audit failures. We select NMN from that output for hackathon submission, retaining the pipeline's `review_required` status and its caveats. Clinical suitability remains unestablished.

This proposal addresses the constitutional chromosome-instability disorder. It does not propose tumour cytotoxicity as a substitute for restoring normal tissue function. No patient-cell efficacy experiment was performed, and no treatment dose or clinical recommendation is made. NMN's qualification under the hackathon's market-approved-medication criterion has not been established by this work.

## Generative AI use

AI coding assistance was used to develop the workflow and prepare this methods description.

**Models used:** Claude Opus 5.5 (Anthropic API), GPT-6 Astra (OpenAI API), commercial terms, no training on customer content.

Case reasoning used **Muse-Glimmer-30B locally under vLLM on four GPUs**, with local evidence tools. The reviewer used a separate context of the same model. Evo 2 and ESM-1v supplied locally computed sequence representations upstream. Gene ranking, partner selection and EPCR assignment were deterministic calculations rather than conversational-model decisions.

## Method abstract

We developed a local workflow connecting variant prioritisation to mechanistic candidate discovery and evidence auditing for mosaic variegated aneuploidy. Track 1 combined Evo 2 DNA features, ESM-1v protein features and VEP annotations in a LightGBM/MLP ensemble trained on public ClinVar labels with gene-grouped splits. BUB1B, CEP57 and TRIP13 were held out of classifier training. Automatic phenotype-guided ranking placed BUB1B first and selected a truncating allele with a low-scoring missense partner hypothesis, preserving unknown phase and residual function.

Track 2 used locally hosted Muse-Glimmer-30B with local gene, pathway, literature, compound-identity and safety resources. No compound list was supplied to automatic discovery. The driver retained discovery sources, required compound-specific evidence and counter-evidence searches, and checked claims through source-membership rules and a separate context of the same model. One bounded correction and fresh review were permitted.

The completed automatic run proposed NMN and NPD11033. Neither passed the final audit: NMN retained evidence-attribution and review failures, while NPD11033 lacked resolved identity. We select the automatically discovered NMN hypothesis for preclinical submission, carrying forward the evidence, uncertainties and unresolved audit findings produced by the run. Published NMN experiments provide a reason to test BubR1 abundance rescue; they do not demonstrate restored chromosome segregation in patient cells or efficacy in human MVA.

Our proposed experiment would establish genotype and residual function, then compare NMN, nicotinamide and vehicle across controlled exposures. Functional chromosome-segregation endpoints, viability and proliferation would accompany protein and pathway measurements. SIRT2 perturbation would test dependency with appropriate controls.

The workflow's strengths are its explicit genetic handoff, exact compound identities, preserved sources and inspectable abstention. Its limitations include known-case development, disease-focused evidence, heuristic genotype interpretation and fallible model review. Adult exposure evidence does not establish chronic paediatric safety. The deliverable is a testable, qualified candidate hypothesis with a reproducible audit trail, not a validated therapy.

## Approach from variant to candidate

### Genetic starting point

The upstream classifier used a 271,170-variant ClinVar subset spanning training, validation and test. Evo 2 supplied DNA likelihood and embedding features; ESM-1v supplied protein scores and embeddings for eligible missense changes. VEP supplied consequence annotations. The validation-selected LightGBM/MLP ensemble excluded gene identity and population frequency as learned features. Its saved held-out AUROC was 0.9778 overall and 0.9694 for missense variants. Training, feature generation and QC are documented in `rare_real/NOTES.md`; the numerical results are in `rare_real/work/model/test_results.tsv`.

The case workflow freshly scored 40,629 prepared candidates using cached sequence features. HPO matching and deterministic partner ordering produced the BUB1B handoff below. Preparation retained additional calls in known MVA genes; development used knowledge of this case, so the complete workflow is not a blinded diagnostic evaluation.

| Handoff item | Recorded result |
|---|---|
| Leading gene | BUB1B, encoding BubR1 |
| Lead variant | c.2210T>G p.Leu737Ter; score 0.9507509286 |
| First partner hypothesis | c.3006T>G p.Asn1002Lys; score 0.0747982326 |
| Genotype interpretation | One damaging allele; second allele unresolved |
| Phase | Unknown |
| Residual function | Unestablished |

The low missense score was retained, not increased to fit a recessive diagnosis. The team reported 100/100 rank points and F-max 1.000 for the automatic Track 1 submission. That competition result is separate from confirmation of allelic phase, residual protein function or therapeutic response.

### Local discovery and review

The automatic Track 2 driver first retrieved gene, disease and protein-regulation context. Muse-Glimmer-30B then proposed candidates without being supplied molecule names. For each candidate, the workflow resolved identity, retrieved supporting and opposing literature, searched adult exposure and paediatric safety, and used local chemical liability models where a structure was available.

Evidence records distinguished the intervention actually tested, organism or tissue, exposure, endpoint, source passage and interpretation. Discovery citations were carried into later assessment so a new search could not silently replace the evidence that generated a candidate. Abstracts and full text were retained separately. Deterministic checks tested source membership and compound attribution; a separate context of the same model reviewed the scientific claims. This was model review, not an independent expert assessment. One bounded correction followed by fresh review was allowed.

The final automatic result was:

| Candidate | Result | Remaining issue |
|---|---|---|
| NMN | `review_required` | Compound-attribution mismatch and unsupported/issues-marked review; nomination withdrawn during correction |
| NPD11033 | `review_required` | Unresolved compound identity |

No candidate passed as a research lead; clinical suitability remained `unestablished`. We preserve that result alongside the team's subsequent decision to submit NMN. The audit outcome is a property of the generated claim set, not an experimental demonstration that NMN cannot work.

## Variant mechanism and therapeutic hypothesis

BubR1 participates in the spindle assembly checkpoint that restrains anaphase until chromosome attachment is appropriate. The truncating BUB1B allele supports a loss-of-function hypothesis. The candidate second allele, phase and remaining protein function are unresolved; the workflow therefore does not establish a complete biallelic molecular mechanism.

Our proposed intervention is conditional: **if sufficient functional BubR1 remains, increasing its stability or abundance might improve chromosome segregation**. The desired outcome is fewer abnormal divisions with continued cell growth. Increasing an inactive or harmful protein product would not meet that goal. This is a hypothesis about preserving function in constitutional disease, not a claim that checkpoint inhibition or killing deficient cells would benefit the syndrome.

NMN was considered through the NAD+/SIRT2/BubR1 regulatory connection. The retrieved primary study supplies compound-specific evidence for an abundance effect, while the transfer to this genotype and to functional rescue remains untested. [North et al., 2014](https://pubmed.ncbi.nlm.nih.gov/24825348/)

## Evidence supporting NMN and its limits

| Evidence | Finding relevant to the proposal | Limit on interpretation |
|---|---|---|
| North et al., 2014, PMID 24825348 | NMN increased BubR1 in mouse testes after 500 mg/kg/day intraperitoneally for seven days. In wild-type fibroblasts, 1 mM NMN for 24 hours increased BubR1; induction was almost completely absent in Sirt2-null cells. | These are experimental exposures, not proposed patient doses. SIRT2 dependence was demonstrated in cultured cells; NMN-mediated rescue of this patient's genotype was not shown. |
| Same primary study | Nicotinamide at 20 mM for 16 hours reduced BubR1 in the studied cell setting. | This was a sirtuin-inhibitory exposure, not a matched therapeutic comparison. It does not establish that all nicotinamide exposures are harmful. |
| Suematsu et al., 2014, PMID 25285631 | SIRT2 knockdown increased BubR1 K250 acetylation without a substantial increase in total BubR1 or delayed prometaphase-to-anaphase progression. | This constrains a simple “deacetylation always stabilises BubR1” explanation. K250 and K668 findings must not be conflated. |
| PMID 35479740 | Placebo-controlled oral NMN study in 30 healthy adults over 12 weeks reported blood NAD changes and short-term tolerability. | No evidence of MVA efficacy or lifelong paediatric safety. |
| PMID 36482258 | Placebo-controlled study in 80 healthy middle-aged adults over 60 days reported NAD changes and tolerability. | Adult exposure does not establish patient-specific suitability or a paediatric regimen. |

Sources: [North et al. full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4194088/), [Suematsu et al.](https://doi.org/10.1016/j.bbrc.2014.09.128), [12-week adult study](https://pubmed.ncbi.nlm.nih.gov/35479740/), [60-day adult study](https://pubmed.ncbi.nlm.nih.gov/36482258/).

NMN and nicotinamide are different molecules. The compound-specific evidence motivates choosing NMN over transferring its findings to nicotinamide. The lifespan benefit in the North study belongs to SIRT2 transgenesis, not NMN treatment. Protein abundance and survival effects from different interventions must remain separate.

Chemical identity resolution and liability screening are safeguards within the automatic workflow. Model predictions cannot establish clinical safety or rank patient tolerability. Missing NMN metadata in the local drug-safety table also cannot establish absence of human exposure or a regulatory status.

The team's nomination therefore rests on a compound-specific preclinical rationale and the ability to test it. It does not rest on demonstrated MVA efficacy, approval for MVA, or a completed automatic evidence-gate pass.

## Proposed validation and decision criteria

The first step would be to establish phase, second-allele interpretation and residual BubR1 expression/function in an appropriate patient-derived cellular model, with an isogenic corrected comparator where feasible. Candidate testing would compare NMN, nicotinamide and matched vehicle across a prespecified exposure range, with replication, exposure verification and blinded endpoint scoring.

Measure NAD-related pathway engagement, BubR1 abundance and half-life, and direct chromosome-segregation errors per division. Monitor viability, proliferation and aneuploidy over repeated divisions. Protein abundance should be measured alongside function so that apparent biochemical rescue cannot conceal impaired growth or persistently abnormal segregation.

SIRT2 inhibition and genetic perturbation would test mechanism dependency, with controls for their own effects. Loss of an NMN benefit after a valid SIRT2 perturbation would support dependence on that pathway; it would not itself falsify the proposal. Inhibitor specificity and cell viability must remain interpretable.

| Decision | Criterion |
|---|---|
| Support further investigation | Reproducible improvement in segregation with preserved viability and proliferation, accompanied by interpretable target-engagement and protein measurements |
| Abundance hypothesis weakened | No BubR1 increase despite adequate exposure and a sensitive assay |
| Functional rescue not established | More BubR1 without improved chromosome segregation |
| Stop advancement | Increased segregation errors, impaired viability/growth, or another adverse effect that defeats the constitutional-disease objective |

These experiments have not been performed. Translation would additionally require disease-relevant safety assessment, including effects on malignant cells and interactions with treatment. Evidence from unrelated tumour systems cannot establish this patient's risk. Adult studies and chemical target predictions do not resolve those questions.

## Automatic execution and team review

**Automated work:** variant scoring, phenotype-guided gene ranking, partner ordering and handoff; unseeded Track 2 discovery; tool dispatch; identity and source checks; model-based correction/review; report generation and conservative final gates. Version 0.5.0 subsequently added automatic competition-CSV export, verified separately against the same ranking and equivalent handoff.

**Team preparation and review:** the team prepared inputs, transcribed observed HPO terms, developed the pipeline and reviewed its outputs for submission. The final automatic run selected the gene and partner hypothesis, generated the NMN candidate and evaluated its evidence without manual candidate selection during execution. The team's submission decision is to advance that automatically generated hypothesis for preclinical investigation despite the unresolved audit findings; the saved gate result remains unchanged.

Report preparation included source-level checking of the generated claims, including separation of compound effects from genetic interventions and distinction between protein abundance and functional rescue. The report retains the automatic audit failures rather than presenting every generated assertion as established evidence. This review is not clinical or experimental validation.

The final run executed candidate discovery and review automatically from the genetic handoff. Workflow development and evidence collection used this known case, so automatic execution should not be interpreted as validation on an unseen patient. The completed scientific run is workflow 0.4.1, with its full discovery, correction and audit records preserved.

## Public and proprietary data sources

We used public biological databases and literature, public model resources, locally trained models, and the access-controlled hackathon case. No additional proprietary biological database or private compound-screening dataset was used in the documented workflow.

| Resource | Role |
|---|---|
| ClinVar, Ensembl GRCh38/VEP 115, gnomAD annotations | Variant labels, sequence/transcript context and population-frequency filtering upstream of the handoff |
| HPO ontology and gene annotations | Observed-phenotype matching and gene prioritisation |
| Open Targets and UniProt | Gene/disease relationships and protein-function context |
| Orphanet and the local disease-resource collection | Related-disease hypotheses and mechanistic context |
| PubMed/Europe PMC abstracts and cached primary full text | Supporting and opposing mechanism evidence, intervention attribution and exposure studies |
| PubChem identity records and ChEMBL 37 | Exact compound identification, structures and drug/target metadata |
| Locally trained chemical target models | Liability screening using 1,501 Laplacian naive-Bayes models over unfolded ECFP4; training resources included ChEMBL and ExCAPE-DB |

The broader tool collection also contains Reactome and DepMap resources. Their availability does not imply that every resource was decision-bearing in the final NMN assessment; the saved tool log identifies actual calls. The literature collection and some retrieval axes were MVA-focused, limiting generalisation and search completeness.

Patient inputs, case inference and case-run tool retrieval remained local. Public reference snapshots were prepared separately. Offline library settings were used but are not an operating-system network-isolation guarantee. The access-controlled case is not described as unrestricted public data, and its raw files are not required in a public methods repository.

## Time effort and reproducibility

The recorded development and analysis spanned 12–16 October 2026. The final automatic reasoning job made **60 logged tool calls** and completed in **9 minutes 48 seconds of Slurm wall time**. This excludes reference preparation, model training, feature extraction and earlier iterations. Total human effort and monetary cost were not recorded.

For upstream context, the training-subset ESM-1v extraction took approximately 26 minutes on four H100 GPUs; Evo 2 extraction took approximately 5 hours 40 minutes on four GPUs. The demonstrated workflow reused existing features and reran the trained scorer. These figures must not be summed as a measured complete-run runtime.

Verification recorded exact reproduction of all 40,629 scores, stable ranking under input permutation, validated checkpoint hashes and completed-run resume with an identical report. Version 0.4.1 had 21 recorded passing regression tests; the subsequent export checks recorded 28. These support implementation reproducibility, not therapeutic validity. Feature regeneration through the new adapter and independent disease-case generalisation remain unvalidated.

Reproducibility artifacts are retained under `rare_real/work/patient/unified_20261016/`: `e2e_enriched_v4/` contains the exact scientific-run source snapshot, manifest, structured results, prompts, tool logs and source ledger; `automatic_submission/` contains the later export and verification records. The source entry point is `mva_llm/unified/`, with usage in its README. Current code is 0.5.0; the historical scientific result came from the preserved 0.4.1 snapshot.

## Potential impact and limitations

The potential biological contribution is a concrete experiment connecting a ranked genetic lesion to an intervention on residual protein abundance. If functional rescue were demonstrated, it would motivate further disease-specific development. No effect size, patient benefit or clinical feasibility can yet be estimated.

The methodological contribution is automatic progression from genomic scoring through mechanism retrieval to compound discovery, evidence assessment and explicit uncertainty reporting. The system generated both a concrete candidate and reasons to scrutinise it. The team submits NMN as that qualified research hypothesis, preserving the distinction between candidate selection and a passed evidence audit.

The framework can be adapted to other cases, but broader use requires disease-independent input preparation, stronger inheritance reasoning, appropriate evidence coverage and evaluation on unseen cases. Neither a perfect score on one competition case nor held-out variant AUROC validates drug discovery. For NMN, the decisive gaps are residual function, functional rescue, tissue-specific response, chronic paediatric safety and patient-specific clinical context. The submission presents these gaps as testable constraints on the hypothesis.

## Methods template and acknowledgement

This report covers all Track 2 fields in the [official methods description workbook](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/resolve/main/static/templates/methods_description_form.xlsx), checked on 19 October 2026: team, detailed approach, AI disclosure, automation, manual review, public/proprietary sources, variant mechanism, effort and an abstract under 500 words. The [Track 2 submission instructions](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track2.py) accept Markdown or PDF and separately request a repository and three-minute pitch video.

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.
