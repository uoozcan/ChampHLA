# MVHLA Manuscript — Submission Checklist

Last updated: 2026-06-04

---

## HUMAN TASKS (cannot be automated)

### Data & Code Deposition
- [ ] **Register Zenodo DOI** for benchmark data tables (`method_comparison.tsv`, `tool_confidence_weights.tsv`, `bimodal_accuracy_comparison.tsv`, `trimodal_accuracy_comparison.tsv`, etc.). Update placeholder in manuscript Data Availability section.
- [ ] **Tag and push pipeline commit** `4120c6c249f8dda1aedabc0db87b8abbce02ad75` to GitHub with a release tag (e.g., `v1.0.0`).
- [ ] **Confirm GitHub repository license** — MIT or Apache 2.0 recommended for Nextflow workflows. Add `LICENSE` file to repo. Update placeholder in manuscript.
- [ ] **Confirm IMGT/HLA 3.59.0 reference files** are deposited or linked (Zenodo or EBI). Add URL to Data Availability.

### Manuscript Front Matter
- [ ] **Author affiliations** — complete all institutional addresses and ORCID IDs.
- [ ] **Acknowledgements** — draft acknowledgements (funding agency, HPC centre CSC Puhti allocation, 1000 Genomes consortium).
- [ ] **Corresponding author contact** — email address and postal address.
- [ ] **Competing interests statement** — confirm no competing interests or disclose.
- [ ] **Funding statement** — grant numbers and funding agencies.

### Figures
- [ ] **Figure 1** — workflow architecture schematic is currently a manual diagram placeholder. Must be created as a publication-quality vector figure (SVG/PDF, Arial font, colourblind-safe palette). See `figure_hub/plot_style.py` for standards.
- [ ] **VENEX figure removal** — remove `mvhla_figures_v6/venex_fig1-5` from the submission figure set. Confirmed these are internal QC only (see Supplementary Note 1 in manuscript).
- [ ] **Final figure audit** — verify all 9 main + 3 supplementary figures exist in `analysis/figures_final_candidate/` (pdf + svg + png at 300 DPI). Run `ls analysis/figures_final_candidate/` to confirm.

### Journal Selection
- [ ] **Confirm target journal** — current best fit: *Bioinformatics* (Oxford). Alternatives: *GigaScience*, *Briefings in Bioinformatics*, *PLOS Computational Biology*.
- [ ] **Abstract word count** — Bioinformatics limit is 200 words. Current abstract is ~300 words — needs trimming before submission.
- [ ] **Check structured abstract requirement** — Bioinformatics requires Motivation / Results / Availability and Implementation sections (not Background/Results/Conclusions). Reformat if submitting there.

---

## MANUSCRIPT ITEMS REMAINING

### Content gaps
- [ ] **Table 1 — confidence type verification**: POLYSOLVER and Seq2HLA are listed as "no confidence" (no native confidence output). Verify against tool documentation that no score is available; add a footnote if a score exists but was not implemented.
- [ ] **DRB1/DQB1 scope decision**: Results §9 reports bimodal DRB1 (n=38) and DQB1 (n=35) as exploratory. Decide whether to (a) expand truth acquisition for these loci and report as headline results, or (b) formally move to a supplementary table and out-of-scope note.
- [ ] **Supplementary Methods**: expand weight formula derivation (currently in Methods §4 prose only); consider adding equations for `base_reliability`, `guardrail_factor`, and `final_weight`.
- [ ] **Supplementary data tables**: decide which benchmark TSVs to deposit as supplementary tables vs Zenodo-only.

### Cross-references to verify
- [ ] All figure references (Fig 2–10, S1–S3) match the actual figure filenames in `figures_final_candidate/`.
- [ ] Table 3 RNA-seq row for Seq2HLA: n=29, callable=1.0, accuracy=0.581. Confirm this reflects 29 samples with truth, not 107.
- [ ] Results §3 benchmark weight narrative: re-read against `tool_confidence_weights.tsv` (wave2) to confirm all numbers match.

---

## DATA/CODE READINESS CHECKS

Run these before submission:

```bash
# 1. All main + supp figures present
ls /scratch/project_2008084/pihla-publish/analysis/figures_final_candidate/*.png | wc -l
# Expected: ≥12

# 2. No stale numbers in manuscript
grep -n "n=129\|0\.9487\|0\.9533\|8\.0%.*WGS\|96\.16\|94\.79\|≤0\.1 pp" \
  /scratch/project_2008084/pihla-publish/docs/MANUSCRIPT_DRAFT_V1.md
# Expected: no output (all corrected)

# 3. Figure registry — all active figures have arial_font=yes
grep "partial\|	no	" \
  /scratch/project_2008084/pihla-publish/figure_hub/figure_registry.tsv | grep -v "removed\|legacy"
# Expected: no output

# 4. Confirm pipeline commit is accessible
cd /scratch/project_2008084/pihla-publish && git log --oneline | grep 4120c6c
```

---

## ESTIMATED REMAINING TIMELINE

| Task | Who | Effort |
|---|---|---|
| Figure 1 architecture diagram | Human | 2–4 h |
| Abstract reformat for Bioinformatics | Human + Claude | 30 min |
| Author affiliations + acknowledgements | Human | 1 h |
| Zenodo DOI registration | Human | 1 h |
| GitHub license + release tag | Human | 30 min |
| DRB1/DQB1 scope decision | Human | — |
| Final proofreading round | Human + Claude | 2–3 h |

---

## STATUS SUMMARY (as of 2026-06-04)

| Component | Status |
|---|---|
| All main figures (2–10) | ✅ Arial, hub-native, registered |
| All supplementary figures (S1–S3) | ✅ Arial, hub-native |
| FIMM concordance figures | ✅ Arial, de-identified |
| FIMM-LOH figure | ✅ Arial, de-identified |
| Figure 1 (architecture) | ❌ Manual — not yet created |
| Figure registry | ✅ 26 figures, all active=yes have arial_font=yes |
| FIMM pseudo-IDs | ✅ 1,227 IDs, 6 files de-identified |
| Manuscript data errors | ✅ 15 errors corrected |
| Tables 1–6 | ✅ All drafted (Tables 4–6 from authoritative TSVs) |
| Results §9/§10 per-gene numbers | ✅ Corrected from bimodal/trimodal_method_per_gene.tsv |
| Zenodo DOI | ❌ Not yet registered |
| GitHub license | ❌ Not yet confirmed |
| Abstract word count | ⚠️ ~300 words — needs trimming for Bioinformatics |
| Author affiliations | ❌ Not yet complete |
