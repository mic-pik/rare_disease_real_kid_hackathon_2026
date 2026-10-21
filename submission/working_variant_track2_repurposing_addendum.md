# Track 2 addendum: expression-reversal screen and nomination review

**Team:** Working variant  
**Review date:** 19 October 2026  
**Decision:** No replacement nomination justified by the four reviewed LINCS leaders. Add the public-data screen and its negative translational review to Track 2. The historical NMN report remains a qualified preclinical hypothesis with unresolved market-approved-medication eligibility; it is not upgraded by this work.

## Generative AI and provenance

This follow-up was implemented and reviewed with the coding assistant using public GEO/LINCS data, the frozen FDA/ChEMBL ingredient catalogue, regulatory labels and primary literature. It is a separate computational and source-review extension, not output from the original autonomous Muse-Glimmer discovery run. No new patient-cell experiment or clinical recommendation is claimed. The original nomination report and automatic audit records remain unchanged. This addendum has not been uploaded to the competition.

## New evidence

GSE22206 provides expression arrays from one affected donor carrying a homozygous BUB1B splice variant, controls and relatives. This is not the hackathon patient. Tissue-matched descriptive contrasts yielded 685 direction-consistent genes, of which 74 are directly measured LINCS landmarks. The FDA-intersected LINCS comparison scored 544 eligible ingredients across 2,630 quality-filtered profiles. Twenty-nine passed the predefined computational screen and twenty also passed a post-hoc matched-drug-background check. Methods and full results: [GSE22206 analysis](../docs/GSE22206_ANALYSIS.md), [LINCS comparison](../docs/GSE22206_LINCS_REVERSAL.md).

Those tests identify transcriptional reversal across separate experiments. They do not show drug rescue of the donor's cells. The background check does not adjust for measured cytostasis or cell death. This review examines the four numerical leaders only; the other sixteen background-supported ingredients have not received equivalent exposure/mechanism review.

## Candidate-specific review

### Nisoldipine — do not replace NMN with the current signal

The LINCS signal is at 10 µM in MCF7 for 6 hours and PC3 for 24 hours. Using the catalogue parent structure, molecular weight 388.420 g/mol, this nominal culture concentration equals 3,884 ng/mL.

A primary pharmacokinetic study in 20 healthy subjects receiving a single 5-mg sustained-release tablet reported mean total Cmax 1.4 ± 1.2 ng/mL (Table 2), or approximately 0.0036 µM. The nominal LINCS concentration is about 2,774 times that particular mean total peak. This is a formulation- and dose-specific reference, not the maximum possible exposure or an estimate of intracellular concentration. [Kang, 2006](https://www.jstage.jst.go.jp/article/analsci/22/12/22_12_1597/_pdf).

The label separately reports very high binding, less than 1% unbound in its tested concentration range, and substantial first-pass metabolism. We do not extrapolate that binding experiment into a precise low-concentration unbound estimate. [Nisoldipine label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=a52b773a-44f9-442f-a4be-39f79a558bd2).

No compound-specific BUB1B rescue experiment was identified in the bounded search. Calcium-channel pharmacology alone does not establish correction of chromosome segregation. The large exposure disconnect and missing mechanistic bridge prevent a replacement nomination from these data. A dose-specific reversal at much lower exposure would be new evidence; it is not present in the qualifying profiles used here.

### Nitazoxanide — defer; tested-species mismatch and growth-arrest alternative

The two qualifying LINCS profiles use nominal nitazoxanide at 10 µM for 6 hours in A549 and PC3. The label states that parent nitazoxanide is not detected in plasma after oral administration; circulating exposure is to tizoxanide and its glucuronide. It reports strong protein binding of tizoxanide (>99% in distribution and >99.9% in the interaction section). Thus, a metabolite plasma concentration cannot be presented as exposure to the LINCS parent compound. Conversion within LINCS cultures and the responsible intracellular species were not measured here. [Nitazoxanide label](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=c0dc0656-3c46-43f9-adea-5285c522cef8&version=8).

Primary work reports NTZ/TIZ-mediated mitochondrial uncoupling, reduced ATP production, E2F1 suppression and G1 arrest in prostate-cancer models. This supplies a concrete alternative explanation for reversing a proliferative expression signature, not evidence of restored segregation. [Prostate-cancer study, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11866447/). Separate colorectal experiments also found G1 arrest and proliferation inhibition; importantly, those experiments reported different effects in tumoroids and primary intestinal organoids, so they do not establish universal toxicity in normal cells. [Colorectal study](https://pubmed.ncbi.nlm.nih.gov/31844249/).

No direct BUB1B rescue was identified. A relevant metabolite-specific signal with preserved division and survival would be needed to distinguish a useful mechanism from growth suppression. We do not nominate it as a better-supported replacement.

### Methimazole — exposure overlap exists, but benefit rationale is missing

LINCS reversal appears at 10 µM for 6 hours in MCF7 and PC3. With molecular weight 114.173 g/mol, this is 1,142 ng/mL. An adult study reported mean peaks of 213 and 299 ng/mL after a 10-mg dose, approximately 1.87 and 2.62 µM. [Adult PK study](https://pubmed.ncbi.nlm.nih.gov/3830069/).

It would be incorrect to conclude that 10 µM can never be reached: a study in nine thyrotoxic children/adolescents reported peaks of 4.4–12.6 µM, median 9.2 µM, following its specified body-surface-area-based dose. These observations are from treatment of thyroid disease, not MVA, and a transient peak does not establish a sustained six-hour cellular exposure. No dose is proposed here. [Pediatric PK study](https://pubmed.ncbi.nlm.nih.gov/3591220/).

The label identifies hypothyroidism and potentially life-threatening agranulocytosis among material risks. [Methimazole label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=26ce3ca1-f2f2-7102-e063-6394a90ae819). These risks are not an automatic ban on research; the problem is that the review found no BUB1B-specific benefit mechanism or functional rescue to balance them. The expression signal alone does not justify replacing NMN with an antithyroid intervention for constitutional chromosome instability.

### Chlorhexidine — retain as a screen result, not a systemic MVA nomination

The numerical leader reverses the signature at 10 µM for 24 hours in MCF7 and PC3. The reviewed approved oral-rinse product acts locally, is poorly absorbed from the gastrointestinal tract, and is explicitly not intended for ingestion. [Chlorhexidine label](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=ec1c2a9f-f0e3-c560-360a-2d2cbdd0cb7c&version=9).

The label's historical ingestion study is not an approved systemic regimen. We do not convert its salt-related exposure into a parent-molecule therapeutic margin. There is no demonstrated route/exposure bridge from this local product to systemic BUB1B rescue. Its high numerical score should remain visible as a limitation of ranking by expression reversal alone.

## Submission consequence

The review does not establish any of these four as a better-supported Track 2 nomination. The defensible submission update is the completed FDA-focused screen plus this explicit candidate triage, with no forced new winner. NMN retains its narrower compound-specific BubR1-abundance rationale, but its functional-rescue and competition-eligibility gaps remain unresolved. Failure of these alternatives does not resolve those gaps.

If a market-approved medication is mandatory, the current work still does not supply a sufficiently supported replacement. Approval of an ingredient for another indication resolves only part of eligibility; it does not make its proposed MVA rationale credible. The competition's current acceptance rules and submission window were not independently verified in this review.

Suggested report wording:

> A subsequent public-data repurposing analysis identified FDA-approved ingredients whose LINCS profiles opposed an expression signature from a separate BUB1B-affected donor. Of 544 scored ingredients, 29 passed an exploratory screen and 20 retained support in a post-hoc matched-drug comparison. Focused review of four leading hits found unresolved exposure, molecular-species, mechanism or growth-confounding problems. We therefore report them as computational hypotheses without replacing the original qualified NMN proposal or claiming functional rescue. NMN's eligibility under the market-approved-medication criterion remains unestablished.

## Scope and verification

This is a bounded assistant-led desk review, not a systematic review or independent clinical assessment. Exact BUB1B/BubR1 rescue searches and compound-specific PK, label and cell-cycle searches were inspected; not finding a study does not prove it does not exist. Primary evidence supports the exposure and counter-evidence claims; the nomination decisions are our interpretation.

Molecular weights were calculated with RDKit from the frozen FDA catalogue parent structures. Conversions use `µM × molecular weight = ng/mL`; exposure ratios compare nominal culture concentration with study-specific total plasma peaks and are not safety margins. No culture free fractions, intracellular levels or measured-growth adjustments are available. The full nisoldipine FDA review PDF was archived but could not be text-inspected in the available environment; its search-index numerical excerpts were not used for the final exposure calculation. The smaller primary PK paper was text-inspected, including Table 2. Other sources were browser-inspected; they are linked rather than represented as archived full texts.
