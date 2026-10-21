# Repository verification

Verified on 19 October 2026 against the saved final-run artifacts.

## Completed checks

- All 28 original workflow unit tests and five portability tests passed.
- Python compilation, shell syntax checks and the complete staged dry-run passed.
- The portable workflow replayed Track 1 using the original prepared inputs and resources. Ranking, handoff and submission-audit JSON matched the saved automatic-export run exactly; the submission CSV matched byte for byte.
- All 11 local-tool schemas, including their order, matched the original final implementation. Checked public-resource responses for gene lookup, gene–disease associations, target drugs and drug safety also matched.
- The full-text parser reproduced all 75 paragraphs of the saved public source, including scientific notation.
- All 12 configured public source URLs responded successfully to availability checks. Availability alone does not establish that a fresh download matches the recorded checksum.
- Repository scanning found no likely credential values, machine-specific home paths, checkpoints or raw genomic inputs. The pre-existing reachable Git blobs were also checked for credential patterns.

The patient replay artifacts are retained outside this repository. Raw patient inputs remain excluded. The subsequently added `submission/` folder contains the explicitly requested original submission documents and selected historical scientific-run outputs; its artifact manifest records original and packaged hashes.

## Not rerun during packaging

The complete foundation-model feature extraction, supervised and chemical-model training, and final GPU reasoning execution were not repeated from scratch. The published numerical results describe the historical execution. Unit tests and the exact Track 1 replay do not establish identical output from a newly trained ensemble or a new LLM execution.

The repository pins model revisions, records dependency constraints and resource signatures, and fails on missing or changed evidence resources by default. A fresh rebuild still depends on upstream resource availability and compatible hardware/runtime libraries. Generated run logs and resource-check reports should accompany any future full reproduction claim.
