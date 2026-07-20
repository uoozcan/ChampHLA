# ChampHLA Manuscript — Submission Checklist

Last updated: 2026-07-16 (supersedes the 2026-06-04 MVHLA-era checklist)

Target: methods/software journal (*Bioinformatics* Application Note, *GigaScience*, or *BMC Bioinformatics*).

## Current manuscript state

- **Main manuscript:** `docs/CHAMPHLA_MANUSCRIPT_V3.md` / `.docx` — 7 main figures, 10 tables.
- **Supplementary:** `docs/CHAMPHLA_MANUSCRIPT_V3_SUPPLEMENTARY.md` / `.docx` — Supplementary Figures S1–S14, Supplementary Tables S1–S3, Supplementary Notes S1–S8 (each series numbered independently).
- **Abstract:** full abstract + short abstract (197 words, within a 200-word limit).
- **Headline result:** held-out 10-fold nested cross-validation; WGS Champion-Challenger 0.5012 vs MajorityVote 0.3844 (exact McNemar p<0.0001); WES/RNA parity; ablation shows the gain is routing, not weighting.
- **Reproducibility:** `bin/reproduce_manuscript_numbers.sh`, `slurm_nested_cv.sh`; shipped `conf/tool_weights_{wgs,wes,rna}.json` regenerated from the deposited `tool_confidence_weights.tsv` so production defaults match Table 5.
- **Repo:** https://github.com/uoozcan/ChampHLA, branch `roihu-migration-locityper-docs`.

## HUMAN TASKS (cannot be automated)

### Front matter
- [ ] **ORCID identifiers** for all authors (currently absent).
- [ ] **Author Contributions / CRediT** statement (currently absent).
- [ ] **Funding statement** — grant numbers + agencies (Acknowledgements currently: "Funding information will be provided prior to submission").
- [ ] **Corresponding-author contact** (email + postal address).
- [ ] Competing-interests statement — present ("The authors declare no competing interests"); confirm still accurate.

### Data & code deposition
- [ ] **Register Zenodo DOI** for the deposited analysis tables + figures; replace the single placeholder `[Zenodo DOI, to be inserted]` in Data Availability. (Data Availability now commits to deposit **at the time of submission** throughout.)
- [ ] **Confirm GitHub repository is public** and add a `LICENSE` (MIT/Apache-2.0 recommended for Nextflow workflows).
- [ ] **Tag a release** on the pipeline commit used for submission.

### Journal-specific
- [ ] **Confirm target journal** and reformat the abstract to its required structure (e.g. *Bioinformatics* wants Motivation / Results / Availability).
- [ ] **Cover letter.**
- [ ] Prepare the submission figure bundle (`bin/build_submission_figures.py`) at the journal's DPI/format.

## Notes for reviewers/co-authors
- Class II (DRB1/DQB1) is implemented but not truth-backed; stated as future work.
- External cohorts (GIAB trio, HG002, IHWG) are small-n confirmatory; NCI-60 is the powered external arm.
- The single significant accuracy gain (WGS) is per-gene champion routing; the override gate and calibration guardrail are auditability/soundness contributions, stated explicitly in Limitations.
