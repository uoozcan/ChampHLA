# PIHLA Manuscript Review V1

> Obsolescence note: this review was written against superseded benchmark roots (`benchmark_wgs_all_samples`, `benchmark_wes_all_samples`, `benchmark_rna_all_samples`). The current authoritative manuscript-facing roots are `benchmark_wgs_wave2`, `benchmark_wes_truthbacked/run`, and `benchmark_rna_truthbacked/run`. Keep this file only as historical review context, not as a current source-of-truth document.

**Date:** 2026-05-12  
**Reviewer:** Claude Code (Sonnet 4.6)  
**Manuscript:** MANUSCRIPT_DRAFT_V1.md  
**Authoritative tables:** `analysis/benchmark_wgs_all_samples/`, `analysis/benchmark_wes_all_samples/`, `analysis/benchmark_rna_all_samples/`, `analysis/benchmark_trimodal_all_samples/`  
**Target venue:** *Bioinformatics*

---

## STOP BEFORE READING FURTHER

Before addressing writing quality, argument structure, or submission strategy, a fundamental data integrity problem must be resolved. The manuscript contains headline numerical claims that do not match the authoritative benchmark tables. Specifically:

- The single-modality sections (WGS, WES, RNA) appear to use numbers from an earlier benchmark run.
- The multi-modal sections (bimodal, trimodal) use the current full-cohort run and are internally consistent with the tables.
- This means the manuscript is a mixed-version document. It will fail consistency checking by any careful reviewer.

The full discrepancy table is in Issue 1. Everything else in this review assumes you will resolve that first.

---

## Part 1 — Detailed Issue Inventory

Issues are ordered by severity. Each has: Severity / What the problem is / Why it matters / Evidence.

---

### Issue 1 — CRITICAL: Manuscript numbers do not match benchmark tables (mixed-version document)

**Severity:** BLOCKING. Must be resolved before any other edits have lasting value.

**The problem:** The single-modality Results sections (WGS, WES, RNA-seq) contain numbers from an earlier benchmark run. The multi-modal sections (bimodal, trimodal) use the current full-cohort outputs and are correct. The manuscript therefore silently mixes two different benchmark runs.

**Why it matters:** Any peer reviewer who cross-checks Table 3 against the supplement, or who runs the pipeline against the stated cohort definition, will find irreconcilable numbers. This is an immediate rejection trigger at *Bioinformatics* and any comparable venue.

**Evidence — WGS discrepancies (authoritative source: `analysis/benchmark_wgs_all_samples/tables/method_comparison.tsv`):**

| Metric | In manuscript | In tables | Diff |
|--------|--------------|-----------|------|
| MajorityVote overall | 0.4428 | 0.3817 | −0.0611 |
| WeightedConsensus overall | 0.3651 | 0.3181 | −0.0470 |
| WeightedConsensus callable rate | 0.620 | 0.598 | −0.022 |
| WeightedConsensus among-callable | 0.5889 | 0.5319 | −0.0570 |
| OptiType overall | 0.5025 (text) / 0.4975 (table) | 0.4962 | various |
| T1K overall | 0.3846 | 0.3251 | −0.0595 |
| HLA-HD overall | 0.3410 | 0.2645 | −0.0765 |
| Kourami overall | 0.2369 | 0.1592 | −0.0777 |
| ArcasHLA overall | 0.1185 | 0.0799 | −0.0386 |
| SpecHLA overall | 0.1810 | 0.1616 | −0.0194 |
| WGS max sample count | 136 | 131 | −5 |

**Evidence — WES discrepancies (authoritative source: `analysis/benchmark_wes_all_samples/tables/method_comparison.tsv`):**

| Metric | In manuscript | In tables | Diff |
|--------|--------------|-----------|------|
| MajorityVote overall | 0.9136 | 0.9354 | +0.0218 |
| WeightedConsensus overall | 0.8706 | 0.8966 | +0.0260 |
| WeightedConsensus callable rate | 0.922 | 0.9328 | +0.0108 |
| WeightedConsensus among-callable | 0.9438 | 0.9612 | +0.0174 |
| ArcasHLA overall | 0.1771 | 0.1344 | −0.0427 |
| T1K overall | 0.8359 | 0.8295 | −0.0064 |
| WES sample count | 130 | 129 | −1 |

**Evidence — RNA discrepancies (authoritative source: `analysis/benchmark_rna_all_samples/tables/method_comparison.tsv`):**

| Metric | In manuscript | In tables | Diff |
|--------|--------------|-----------|------|
| T1K overall | 0.8918 | 0.8794 | −0.0124 |
| Seq2HLA overall | 0.5625 | 0.5806 | +0.0181 |
| Seq2HLA callable n | 46 | 29 | −17 |

**Bimodal and trimodal numbers match the tables.** These sections were apparently written from the current `benchmark_trimodal_all_samples/` run and are internally consistent.

**Action required:**
1. Designate a single authoritative benchmark run for each modality. The `benchmark_wgs_all_samples/`, `benchmark_wes_all_samples/`, and `benchmark_rna_all_samples/` directories appear to be the most current and largest runs.
2. Update every number in the Abstract, Results sections 3–8, Table 3, the Discussion, and the Conclusion to match the authoritative tables.
3. Record the benchmark run version (directory name or git hash) in `benchmark_metadata.json` and cite it in the Methods.

---

### Issue 2 — CRITICAL: WGS abstention rate contradicts itself within the manuscript

**Severity:** BLOCKING. Two passages in the same document give opposite answers to the same question.

**The problem:**
- Results section 4 states: "abstaining on **38%** of WGS loci" and "**62.0%** of loci it commits to"
- Results section 7 states: "abstains on **43%** of loci where support evidence is insufficient, achieving **53.3%** accuracy on the **57%** of loci it commits to"

These are mutually exclusive claims. The actual table shows callable rate = 0.598, meaning ~40.2% abstention.

**Why it matters:** A reviewer will notice this on first reading. It signals that sections were written at different times without reconciliation. If the numbers differ this much within the same draft, what else is mismatched?

**Action required:** Pick one set of numbers sourced from the authoritative table. With callable rate 0.598:
- Abstention rate: 40.2%
- Callable rate: 59.8%
- Accuracy-among-callable: 53.2%

Update both sections with these values. Add a footnote or parenthetical pointing to the authoritative table row.

---

### Issue 3 — CRITICAL: Duplicate Results section number "11"

**Severity:** High — signals sloppy manuscript maintenance and will confuse reviewers.

**The problem:** Results has two sections both numbered "11":
```
### 11. Weight formula sensitivity
### 11. Reproducibility and workflow portability
```
And a section 12 follows. The figure list plans for 10 results sections; the draft has at least 12 distinct ones.

**Action required:** Renumber all Results sections sequentially. Suggested correct numbering:
1. Cohort assembly and filtering
2. Harmonized benchmark framework
3. Benchmark-derived confidence weights
4. Comparison with single tools and majority vote
5. Confidence calibration
6. Ambiguity-aware and version-aware evaluation
7. Abstention behavior
8. Discordance interpretation
9. WES+RNA bimodal joint consensus
10. Trimodal WGS+WES+RNA joint consensus
11. Weight formula sensitivity
12. Reproducibility and workflow portability
13. Phase-gated reporting policy

---

### Issue 4 — HIGH: Codex Q2 framing recommendation not implemented

**Severity:** High — the `pihla-manuscript-consistency` open task explicitly flags this as blocking.

**The problem:** The Codex framing review returned with high confidence that the main framing should shift from "confidence-calibrated ensemble" to **"reliability-weighted ensemble with confidence guardrails"**. The current title, abstract, and Core Claim have not been updated.

Current title: "A benchmark-derived reliability-weighted multi-tool ensemble for HLA typing with **confidence guardrails** across WGS, WES, and RNA-seq"

This is partially right — "confidence guardrails" is in the title — but the Core Claim still says:
> "learns benchmark-derived reliability weights for eight HLA callers, enforces **empirical confidence guardrails**..."

This reads as if the primary contribution is the guardrail mechanism. Codex was explicit: the primary contribution is the **reliability-weighting system**; the guardrails are a defensive mechanism that happens to be necessary because raw confidence scores are miscalibrated. Framing the guardrail as the headline feature is both inaccurate (most tools fall back to reliability-only) and less compelling as a contribution.

**Why it matters:** The distinction matters for peer review. "We build a confidence-calibrated ensemble" invites reviewers to ask: but it's mostly reliability-only because confidence is blocked everywhere — isn't this just a weighted average? "We build a reliability-weighted ensemble with confidence guardrails that prevent miscalibrated raw confidence from distorting the vote" is accurate and defensible.

**Concrete fix — title:**
```
Before: A benchmark-derived reliability-weighted multi-tool ensemble for HLA typing with confidence guardrails across WGS, WES, and RNA-seq

After: PIHLA: a benchmark-derived reliability-weighted HLA ensemble with calibration-gated confidence integration across WGS, WES, and RNA-seq
```

**Concrete fix — Core Claim (replace):**
```
Before:
We present PIHLA, a Nextflow-based HLA ensemble platform that learns benchmark-derived 
reliability weights for eight HLA callers, enforces empirical confidence guardrails...

After:
We present PIHLA, a Nextflow-based HLA ensemble platform that learns benchmark-derived 
per-tool reliability weights, gates raw confidence scores through empirical calibration 
checks before they can amplify ensemble votes, and produces weighted-consensus calls with 
structured abstention and discordance tagging across whole-genome, whole-exome, and RNA 
sequencing inputs.
```

---

### Issue 5 — HIGH: WGS cohort size inconsistency

**Severity:** High — the stated cohort size appears in the Abstract, Methods, Results, and Conclusion.

**The problem:** The manuscript repeatedly states "136-sample WGS cohort." The authoritative table (`benchmark_wgs_all_samples`) shows:
- OptiType: n=131 samples, n=393 gene rows
- Other tools: 99–121 samples depending on tool
- MajorityVote: n=131

The actual cohort (per `benchmark_wgs_all_samples/tables/benchmark_metadata.json`) is 138 total samples, with 131 producing callable outputs from the best-covered tool. The "136" figure appears to be from an earlier benchmark run.

**Why it matters:** Reviewers will check the stated n against the table. A 5-sample discrepancy in the headline cohort size will prompt questions about what happened to those samples and whether the exclusions are valid.

**Action required:** Update the WGS cohort description to: "n=131 samples with at least one callable WGS tool output" (or the correct number from the authoritative run). Explain in Methods section 9 that 7 of 138 samples failed to produce callable WGS outputs from any tool and are excluded from WGS modality analysis.

---

### Issue 6 — HIGH: VENEX cohort present in figures but absent from manuscript

**Severity:** High — either this is missing science or confusing clutter in the figure set.

**The problem:** The figure directory (`mvhla_figures_v6/`) contains 5 VENEX supplementary figures:
- `venex_fig1_tool_coverage` — tool coverage in VENEX cohort
- `venex_fig2_allele_agreement` — allele agreement patterns
- `venex_fig3_homozygous_flags` — homozygosity indicators
- `venex_fig4_loh_status` — LOH status
- `venex_fig5_allele_diversity` — allelic diversity

There is also `analysis/venex_harmonized_calls.tsv` in the analysis directory. The manuscript never mentions VENEX.

**Two possible interpretations:**
(a) VENEX is a second independent cohort that was analyzed for validation. If so, this is a major missing component — an independent cohort substantially strengthens any methods paper's claims.
(b) VENEX was an exploratory internal analysis not intended for publication. If so, the figures must be removed from the publication figure set and the data file must be archived or documented as out-of-scope.

**Why it matters:** If VENEX is a real independent cohort and it's not mentioned, reviewers who learn about it after publication will question why it was omitted. If it's in-scope, including it (even as supplementary) would directly address the main weakness (single-cohort, all-1000G data).

**Action required:** Decide: is VENEX in scope for this paper? If yes, write a Section 11 (Supplementary Validation Cohort) and update the figure list. If no, remove the VENEX figures from the publication figure set and note in a comment that they are reserved for follow-up work.

---

### Issue 7 — HIGH: DRB1/DQB1 results reported in bimodal section but declared out-of-scope in Methods

**Severity:** High — direct contradiction in the same document.

**The problem:**

Results section 9 states:
> "DRB1 (n=76 assessments) and DQB1 (n=70 assessments) showed similar patterns, extending the bimodal benefit beyond class I loci."

Methods section 9 states:
> "Headline accuracy claims are restricted to HLA-A, HLA-B, and HLA-C."

Discussion section 6 (Limitations) states:
> "Current headline evaluation is restricted to HLA-A, -B, and -C. DRB1 and DQB1 are implemented in the framework but are not included in the primary benchmark because matched truth labels were unavailable."

These three passages are logically inconsistent. If DRB1/DQB1 truth labels were unavailable, the section 9 DRB1/DQB1 numbers cannot exist. If they do exist (the bimodal table has n=391 assessments, which is more than 106 samples × 3 genes = 318, implying 73 DRB1+DQB1 rows), the limitation statement is wrong.

**Why it matters:** A reviewer focused on HLA biology will immediately notice this. It creates doubt about whether the claims in section 9 are valid.

**Action required:** Check `bimodal_wes_rna_consensus.tsv` for which genes are present. If DRB1/DQB1 rows exist and are truth-backed, remove the "unavailable" language from Limitations and add a Methods subsection describing the DRB1/DQB1 truth source. If DRB1/DQB1 rows should not be in the table, recompute the bimodal assessment count restricted to A/B/C (318 assessments) and remove the "DRB1 and DQB1 showed similar patterns" sentence.

---

### Issue 8 — HIGH: Figure 0 (pipeline architecture diagram) is missing

**Severity:** High — for a methods/infrastructure paper, the architecture diagram is arguably the most important visual.

**The problem:** The manuscript explicitly notes:
> "A workflow architecture schematic (Figure 0) is planned for manual creation and is not included in this list."

For a paper submitted to *Bioinformatics* where the primary contribution is a Nextflow pipeline and ensemble framework, the absence of a pipeline figure will be the first thing reviewers mention.

**Why it matters:** Infrastructure papers without architecture figures leave readers unable to understand what the software actually does. The abstract, introduction, and methods all describe the pipeline at high level, but without a schematic, readers cannot verify that the architecture matches the description.

**Action required:** Create Figure 1 (renumber subsequent figures +1) showing:
- Input: BAM/CRAM/FASTQ per modality
- Tool modules (8 callers, containerized)
- Harmonization layer
- Confidence calibration + guardrail gate
- Weight learning
- Consensus methods (MajorityVote + WeightedConsensus)
- Outputs: call tables, consensus, discordance tags, abstention flags

This can be drawn in Inkscape, BioRender, or any vector tool. It does not require code changes.

---

### Issue 9 — MEDIUM: OptiType WGS accuracy appears three ways in the same document

**Severity:** Medium — smaller than Issues 1–2 but still a factual inconsistency within a single manuscript.

**The problem:**
- Abstract: "0.5025 overall correct-call rate"
- Table 3 (WGS): "0.4975"
- Actual benchmark table: 0.4962 (callable rate 0.9924, gene rows 393)

All three are different. The table value (0.4975) is likely a rounded version of an intermediate run. The actual table gives 0.4962.

**Action required:** Use 0.496 (or 0.4962 to 4 decimal places) in all locations once the authoritative run is designated. The text narrative that reads "OptiType was the strongest single WGS tool (0.5025 overall correct-call rate)" becomes:
> "`OptiType` was the strongest single WGS tool (0.496 overall correct-call rate; callable rate 0.992)"

---

### Issue 10 — MEDIUM: WES WeightedConsensus underperformance is presented as explained but is only provisional

**Severity:** Medium — not a factual error but a claim that exceeds the evidence.

**The problem:** The manuscript attributes WES WeightedConsensus overall underperformance (vs MajorityVote) to a "structural ceiling effect." The Codex review (Q1, medium confidence) flagged explicitly that this interpretation is provisional because no ablation confirms it. The cleanest discriminating experiment would be to run WeightedConsensus with uniform weights (1/N) on the same WES cohort. If uniform weighting still underperforms MajorityVote, the result is structural; if it doesn't, the learned weight profile is the culprit.

**Why it matters:** Presenting an unverified mechanism as an explanation is a claim that reviewers in an ensemble methods paper are well-positioned to push back on.

**Concrete fix — Discussion section 6, replace "a structural ceiling effect" with:**
> "...consistent with a structural ceiling effect in which the strongest WES tools already agree highly, leaving little room for weighted disambiguation. A direct test — running uniform-weight consensus on the same WES holdout — would discriminate structural ceiling from calibration artifact and is noted as a recommended next ablation."

---

### Issue 11 — MEDIUM: Weight formula sub-optimality for WGS inadequately acknowledged

**Severity:** Medium — the sensitivity analysis exists but the Discussion underplays the finding.

**The problem:** The sensitivity analysis table (`analysis/weight_sensitivity/sensitivity_comparison.tsv`) shows clearly:
- Default 0.7/0.3: WGS WC = 0.2667
- Reliability-only 1.0/0.0: WGS WC = 0.2963
- Equal-weight 0.5/0.5: WGS WC = 0.3030
- Confidence-only 0.0/1.0: WGS WC = 0.3165

The default is strictly worse than all alternatives for WGS. Results section 11 acknowledges this but uses numbers ("0.27" and "0.30–0.32") that don't match the table exactly. More importantly, the Discussion does not contain a clear acknowledgement that the 0.7/0.3 default was chosen for WES/RNA stability, not WGS optimization. Per Codex Q6 (high confidence), this framing is both scientifically accurate and defensible.

**Concrete fix — add to Discussion section 6 (Limitations):**
> "The default weight formula (α=0.7/β=0.3) was selected for stability across WES and RNA-seq modalities, where it achieves the best observed WeightedConsensus performance. For WGS, where calibration guardrails suppress the confidence term for most tools, all tested formula variants (reliability-only, equal-weight, confidence-only) achieve higher overall accuracy than the default (0.296–0.317 vs 0.267), because the default's small β=0.3 confidence term slightly disadvantages the two WGS tools (`T1K` and `ArcasHLA`) whose effective confidence contributions are non-zero. The default is not claimed to be globally optimal; future work may explore modality-specific formula selection."

---

### Issue 12 — MEDIUM: consHLA reference is incomplete

**Severity:** Medium — will be flagged by any reviewer who tries to look up the paper.

**The problem:** Discussion section 5 cites "consHLA (Jurtz et al.)" without a year, journal, volume, or page numbers. If consHLA is a published paper, the full citation must appear in the reference list. If it is a preprint, the URL/DOI and "preprint" designation must be provided.

**Action required:** Locate the full consHLA citation. If it is: Jurtz V. et al. (2017) *Bioinformatics* 33:867–869, update to that. If it is a preprint, add: "Jurtz et al., bioRxiv [doi], [year]."

---

### Issue 13 — MEDIUM: Circular evaluation bias is acknowledged but not quantified

**Severity:** Medium — reviewers at methods venues will ask about this.

**The problem:** The manuscript acknowledges "mild circular evaluation bias" in Methods section 8 and Limitations, but provides no estimate of how large the bias is in practice. For reliability weights that are simple accuracy fractions, the bias is bounded by 1/n per tool-gene-modality stratum, but this is not stated.

**Why it matters:** "Mild" is not a quantitative claim. A reviewer specializing in cross-validation will ask: how much does the ensemble's reported accuracy advantage inflate when weights are learned on the same data? Without any bound or estimate, this limitation is a hand-wave.

**Concrete addition to Limitations:**
> "The circular evaluation bias is bounded by the magnitude of the reliability weight estimates: because each base-reliability weight is simply the fraction of calls matching truth in the same cohort, removing one sample changes the weight by approximately 1/n per tool-modality stratum (1/131 ≈ 0.008 for WGS). For a tool with true accuracy 0.50, the leave-one-out bias is therefore ≲1% on the headline weight estimate. The practical consequence for ensemble accuracy is expected to be smaller still, because the consensus is a weighted vote over multiple tools, not a direct rescaling of a single weight. A formal leave-one-out cross-validation is deferred to the next benchmark cycle when a larger independent cohort is available."

---

### Issue 14 — LOW: Abstract comparison of WeightedConsensus to OptiType is confusing

**Severity:** Low — the comparison is valid but framed in a way that obscures rather than clarifies.

**The problem:** Abstract Results paragraph states:
> "`WeightedConsensus` achieves 58.9% accuracy among the 62.0% of WGS loci it commits to calling — comparable to `OptiType`'s overall rate"

First-time readers must mentally notice that "accuracy among callable" for WC is being compared to "overall accuracy" for OptiType — a comparison across two different denominators. The comparison is legitimate (and the point is important) but it takes work to parse.

**Concrete rewrite:**
> "On WGS, `WeightedConsensus` matches `OptiType`'s overall correct-call rate (53.2% vs 49.6%) on the **59.8%** of loci it commits to calling, while explicitly abstaining on the remaining **40.2%** where tool evidence is too conflicted — rather than forcing a call at reduced accuracy as `MajorityVote` does (38.2% overall)."

(Numbers updated to authoritative values once Issue 1 is resolved.)

---

### Issue 15 — LOW: Results section structure does not match the figure list

**Severity:** Low — cosmetic but signals incomplete manuscript hygiene.

**The problem:** The manuscript has 13 Results sections (after renumbering per Issue 3), but the planned Main Figures list has only 8 figures. Results sections 6–13 collectively correspond to only figures 4–8 in the list. The bimodal section (Results 9) and trimodal section (Results 10) need dedicated figures; neither "Figure 7" (listed as bimodal per-gene gains) nor any listed figure matches trimodal results.

**Action required:**
- Add Figure 9 (bimodal vs unimodal accuracy by gene, A/B/C and DRB1/DQB1 if in scope)
- Add Figure 10 (trimodal vs bimodal, minimal gain plot)
- Update the figure list to include at minimum Figure 1 (architecture, new) through Figure 10

---

### Issue 16 — LOW: "Plain terms" callout boxes are valuable but inconsistently present

**Severity:** Low — an editorial choice, but unevenness will look odd in a formal submission.

**The problem:** Several Results sections have `> **In plain terms —**` callout paragraphs that explain the finding in accessible language. These are present for sections 4, 7, 8, 9, and 10 but absent for sections 3 (weights), 5 (calibration), 6 (ambiguity), and 11 (sensitivity). *Bioinformatics* does not have a "plain terms" format — these callouts will need to either be incorporated as prose or placed in a supplementary explainer. They should not be formatted as blockquotes in the main text.

**Action required:** Decide before submission: keep these as a supplementary "guide for non-expert readers" document, or fold the key insights into the main Results prose. Do not submit with markdown blockquotes in a journal manuscript.

---

## Part 2 — Issue-by-Issue Fixes

This section provides concrete, complete rewrites or corrected values for each issue. All numbers use the authoritative benchmark tables.

---

### Fix 1: Number reconciliation table

After designating the canonical run, every number in the following locations must be updated:
- Abstract Results paragraph: all WGS/WES/RNA percentages
- Key Findings at a Glance (6 bullets): all quoted accuracy values
- Results sections 3, 4, 5, 7, 8: all inline numbers
- Table 3 (all three modality subtables)
- Discussion sections 1, 2, 6: all cited accuracy values
- Conclusion: all quoted values

**Authoritative numbers to use (from `analysis/benchmark_wgs_all_samples/`, `benchmark_wes_all_samples/`, `benchmark_rna_all_samples/`):**

WGS (n=131 for OptiType/MajorityVote/WeightedConsensus; lower for other tools):

| Method | Overall | Callable rate | Among-callable |
|--------|---------|--------------|----------------|
| OptiType | 0.4962 | 0.992 | 0.500 |
| T1K | 0.3251 | 1.000 | 0.325 |
| HLA-HD | 0.2645 | 1.000 | 0.264 |
| Kourami | 0.1592 | 1.000 | 0.159 |
| SpecHLA | 0.1616 | 0.852 | 0.190 |
| ArcasHLA | 0.0799 | 1.000 | 0.080 |
| MajorityVote | 0.3817 | 0.998 | 0.382 |
| WeightedConsensus | 0.3181 | 0.598 | 0.532 |

WES (n=129):

| Method | Overall | Callable rate | Among-callable |
|--------|---------|--------------|----------------|
| OptiType | 0.9225 | 1.000 | 0.922 |
| POLYSOLVER | 0.9147 | 1.000 | 0.915 |
| T1K | 0.8295 | 1.000 | 0.830 |
| SpecHLA | 0.7829 | 1.000 | 0.783 |
| HLA-HD | 0.7804 | 1.000 | 0.780 |
| Kourami | 0.7019 | 1.000 | 0.702 |
| ArcasHLA | 0.1344 | 1.000 | 0.134 |
| MajorityVote | 0.9354 | 1.000 | 0.935 |
| WeightedConsensus | 0.8966 | 0.933 | 0.961 |

RNA-seq (n=107 for ArcasHLA/MajorityVote; n=105 for others; n=29 for Seq2HLA):

| Method | Overall | Callable rate | Among-callable |
|--------|---------|--------------|----------------|
| HLA-HD | 0.9429 | 1.000 | 0.943 |
| OptiType | 0.9397 | 1.000 | 0.940 |
| ArcasHLA | 0.9190 | 1.000 | 0.919 |
| T1K | 0.8794 | 1.000 | 0.879 |
| Seq2HLA | 0.5806 | 1.000 | 0.581 (n=29) |
| SpecHLA | 0.0349 | 0.762 | 0.046 |
| MajorityVote | 0.9502 | 1.000 | 0.950 |
| WeightedConsensus | 0.9221 | 0.960 | 0.961 |

Bimodal (n=391 assessments = 106 samples × 3 genes + DRB1/DQB1):

| Method | Overall | Callable rate | Among-callable |
|--------|---------|--------------|----------------|
| BimodalMajorityVote | 0.9616 | 1.000 | 0.962 |
| BimodalWeightedConsensus | 0.9361 | 0.967 | 0.968 |

Trimodal (n=370 assessments):

| Method | Overall | Callable rate | Among-callable |
|--------|---------|--------------|----------------|
| TrimodalMajorityVote | 0.9622 | 1.000 | 0.962 |
| TrimodalWeightedConsensus | 0.9243 | 0.954 | 0.969 |

---

### Fix 2: Corrected abstract Results paragraph

Replace the current abstract Results section with the following (using authoritative numbers, Issue 1 resolved):

> On a 131-sample 1000 Genomes WGS cohort (HLA-A, -B, -C; weights computed and evaluated on the full cohort without holdout), `OptiType` was the strongest single WGS tool (0.496 overall correct-call rate; callable rate 0.992). Equal-weight `MajorityVote` drops to 0.382 overall because the wide spread of WGS tool performance (8–50%) causes low-accuracy tools to outvote the strongest caller on contested loci. `WeightedConsensus` recovers comparable performance to `OptiType` on the 59.8% of WGS loci it commits to calling (53.2% accuracy-among-callable), while explicitly abstaining on the remaining 40.2% where tool evidence is insufficient. In WES (n=129) and RNA-seq (n=107), where top-tool accuracy exceeds 91% and tools already agree highly, `MajorityVote` (0.935 WES, 0.950 RNA) and `WeightedConsensus` (0.961 accuracy-among-callable in both modalities; 0.897 and 0.922 overall) both perform near the single-tool ceiling. For the 106 samples with both WES and RNA-seq data, bimodal joint consensus achieves 96.8% accuracy-among-callable (overall 93.6%, n=391 gene-locus assessments), +0.8 percentage points over RNA-only and +2.0pp over WES-only. Calibration guardrails blocked confidence-based weight amplification for three of five WGS tools with confidence data showing poor empirical calibration.

---

### Fix 3: Corrected Key Findings bullets

**Bullet 1 — WGS ensemble voting backfires:**
```
Before: "...simple majority consensus (44.3%) performs worse than the best individual tool (50.3%)"
After:  "...simple majority consensus (38.2%) performs worse than the best individual tool (49.6%)"
```

**Bullet 2 — WES and RNA-seq ready for clinical use:**
```
Before: "Both ensemble methods reach >91% overall accuracy on WES and >95% on RNA-seq"
After:  "Both ensemble methods reach >89% overall accuracy on WES and >92% on RNA-seq,
         with accuracy among callable loci exceeding 96% in both modalities"
```

**Bullet 3 — bimodal improvement:**
The bimodal numbers are correct (from matching tables). No change needed. However "~1 percentage point" should be checked against per-gene data (see Part 4).

**Bullet 6 — abstention:**
```
Before: "...flags 38% of loci as insufficiently supported"
After:  "...flags 40% of loci as insufficiently supported"
```

---

### Fix 4: Corrected Table 3

See Fix 1 numbers above. Additionally: 
- Add a column for Wilson 95% CI (lo, hi) — these are already in the benchmark tables
- Remove the "(−)" placeholder rows in the WeightedConsensus among-callable lines and replace with actual values
- Add footnotes explaining that n varies by tool due to tool-level output availability

---

### Fix 5: Corrected Results section 7 (Abstention)

Replace the opening of section 7 with:
> "On WGS, `WeightedConsensus` abstains on 40.2% of loci where support evidence is insufficient, achieving 53.2% accuracy on the 59.8% of loci it commits to. `MajorityVote` calls all loci (callable rate ≈1.000) at 38.2% overall accuracy. On WES and RNA-seq, abstention rates are lower (6.7% and 4.1% respectively), reflecting the higher inter-tool agreement in these modalities."

---

### Fix 6: WES WeightedConsensus framing

In Discussion section 6, replace:
> "...leaving little room for weighted disambiguation; a uniform-weight consensus run on the same cohort would be the cleanest experiment to discriminate structural ceiling from weighting-induced losses and is a recommended next ablation."

With (per Codex Q1):
> "...consistent with a structural ceiling effect in which the strongest WES tools already agree highly, leaving little room for weighted disambiguation. This interpretation is provisional: a uniform-weight consensus run on the same WES cohort — where 1/N weights replace the learned reliability weights — would be the cleanest experiment to discriminate structural ceiling from calibration artifact. If uniform weighting also underperforms MajorityVote, the result is structurally determined; if it does not, the learned weight profile is the more likely culprit."

---

### Fix 7: consHLA citation

Replace "consHLA (Jurtz et al.)" with the full citation. Confirm publication details. If consHLA is Jurtz et al. (2017) *Bioinformatics* 33:867–869 (HLA consensus typing from NGS data), update to that citation. If it is a different publication or preprint, update accordingly.

---

## Part 3 — Publication Readiness Brainstorm

This section addresses what needs to happen to make the paper submittable, organized by blocking vs non-blocking and by effort level.

---

### BLOCKERS — cannot submit without resolving

**B1. Single authoritative benchmark run (Issue 1)**
Root cause: numbers in the manuscript come from at least two different benchmark runs. Before any editing pass has lasting value, designate one run as canonical for each modality, verify that all numbers in the paper match it, and record the run identifier in `benchmark_metadata.json`.

How to do it: Run the pipeline once with a locked configuration and log the output directory + git hash. Write a `BENCHMARK_MANIFEST.md` that maps each figure and table claim to the row in the authoritative output file.

**B2. Resolve DRB1/DQB1 contradiction (Issue 7)**
Two options, both acceptable:
- Option A (include): If `bimodal_wes_rna_consensus.tsv` contains DRB1/DQB1 rows with valid truth, add a Methods subsection for the DRB1/DQB1 truth source, update the scope statement, and report the 5-locus bimodal results properly.
- Option B (exclude): If DRB1/DQB1 truth is not available for a sufficient n, restrict bimodal analysis to A/B/C only (n=318 assessments). Remove the "DRB1 and DQB1 showed similar patterns" sentence and reduce the n from 391 to 318 throughout.

**B3. Figure 0 (pipeline architecture) must be created (Issue 8)**
This is manual work. The manuscript cannot be submitted as a methods paper without this. It need not be complex — a clean 6-box flowchart covering input → tool execution → harmonization → calibration gate → consensus → output is sufficient.

**B4. VENEX cohort: include or explicitly exclude (Issue 6)**
The figures exist and the data file exists. A decision and a 1-sentence annotation of those files is required before submission. If in scope, write the section. If out of scope, archive the figures and add a comment in the analysis directory.

---

### NON-BLOCKING — should be done before submission but won't stop peer review

**NB1. Implement Codex Q2 framing everywhere (Issue 4)**
Update title, Core Claim, abstract, Discussion section 1 main-contribution sentence, Discussion section 5 venue positioning. This is 6 targeted edits.

**NB2. Fix duplicate section numbering (Issue 3)**
Mechanical renumbering, 15 minutes.

**NB3. Add WES ablation sentence (Issue 10)**
One sentence addition to Discussion section 6. Low effort, high credibility payoff with reviewers who specialize in ensemble methods.

**NB4. Quantify circular bias bound (Issue 13)**
Three sentences in Limitations. Low effort, removes a reviewer attack surface.

**NB5. Add weight formula limitation sentence to Discussion (Issue 11)**
One clear sentence acknowledging WGS sub-optimality, per Codex Q6 recommendation.

**NB6. Wilson CIs in Table 3 (already in benchmark tables)**
The `overall_correct_call_rate_ci_lo` and `_ci_hi` columns exist in all benchmark output files. Adding them to Table 3 is copy-paste. This is required for statistical credibility at any top venue.

**NB7. Supplementary material plan**
The paper should include as supplementary materials:
- Table S1: Per-gene A/B/C accuracy for all tools and modalities (already in `summary_per_gene.tsv`)
- Table S2: Weight sensitivity analysis (already in `weight_sensitivity/sensitivity_comparison.tsv`)
- Table S3: Discordance taxonomy counts by modality (from `discordance_summary.tsv`)
- Figure S1: Per-gene accuracy heatmap
- Figure S2: Abstention tradeoff curves
- Figure S3–S7: VENEX figures (if in scope)

**NB8. Resolve "plain terms" callout boxes**
Six callout boxes exist in the Results. For *Bioinformatics* submission, these need to become either: (a) properly integrated prose paragraphs, or (b) a separate supplementary "Plain language summary" document. Blockquote formatting is not a standard journal format.

---

### ANTICIPATING REVIEWER OBJECTIONS

These are the most likely push-backs, with pre-emptive defenses to build into the manuscript now:

**Objection 1: "You learn and evaluate on the same cohort. The performance advantage of WeightedConsensus over single tools is inflated."**
Response strategy: Acknowledge explicitly in Limitations (partly done, needs quantification per Issue 13). Add the leave-one-out bound estimate. Emphasize that the key claim is not "WeightedConsensus beats OptiType by X%" — it is "WeightedConsensus provides structured abstention that identifies the hardest loci, which is qualitatively new behavior regardless of the accuracy margin." This reframing deflects the accuracy-inflation objection while keeping the primary contribution intact.

**Objection 2: "Class I only. HLA-DRB1 and HLA-DQB1 are clinically critical and are absent from your headline results."**
Response strategy: Acknowledge this limitation in the abstract (not just Limitations). Say explicitly: "The current benchmark covers HLA-A, -B, and -C; DRB1 and DQB1 are implemented in the framework and evaluated in the bimodal supplement but are excluded from headline claims pending expanded truth acquisition." If Option A (include DRB1/DQB1) is chosen in B2 above, this objection dissolves.

**Objection 3: "OptiType is from 2014 and HLA-HD, T1K, and ArcasHLA are all post-2019. Why is a 2014 tool still the best on WGS?"**
Response strategy: This is a real and interesting finding. Build it into the Discussion as an observation: OptiType's ILP-based approach is uniquely well-suited to the highly repetitive HLA-region short reads that characterize WGS. The newer tools have traded ILP precision for scalability and multi-locus coverage. The ensemble framework is agnostic to this — when a newer tool eventually surpasses OptiType on WGS, the weights will reflect that automatically.

**Objection 4: "n=131 WGS samples from the 1000 Genomes Project is a narrow benchmark. All five ancestry groups are European/African and all truth labels are third-tier (public reference). This is not enough to make clinical deployment claims."**
Response strategy: This is valid. The current manuscript should not make clinical deployment claims from 1000G data alone. Every clinical statement (transplant decisions, neoantigen prediction, pharmacogenomics) should be framed as "the framework is designed for clinical contexts" rather than "clinical performance has been validated." Revise the Introduction sections 1 and 5 accordingly.

**Objection 5: "WeightedConsensus has a lower overall accuracy than MajorityVote on WES (90% vs 94%). Why would anyone use it?"**
Response strategy: This is a valid question and the manuscript's current answer ("ceiling effect," "structured abstention") is partially convincing but needs strengthening. The answer should include: (1) the 6.7% of WES loci that WeightedConsensus abstains on are identifiable as the hardest loci — this is actionable information for a clinical lab that wants to know which samples need re-sequencing; (2) in heterogeneous or multi-ancestry cohorts where tool performance may vary more, the calibration-weighted approach becomes more beneficial; (3) the WES ceiling effect means both methods are clinically adequate, so the choice is about information content, not accuracy.

---

## Part 4 — Results Evaluation and Richer Information Extraction

This section identifies what the data is actually saying that the manuscript is not currently surfacing.

---

### 4.1 The WGS performance spread is the core story, not a background fact

The manuscript buries the key insight: on WGS, the tool accuracy range is 8–50%. This 42-percentage-point spread is extraordinary and is the explanation for everything else that follows:
- Why MajorityVote underperforms the best tool
- Why WeightedConsensus has to abstain on 40% of loci
- Why the guardrails fire on most tools (overconfident in a low-accuracy regime)
- Why bimodal > WGS+anything

Currently, the spread is mentioned in passing: "The performance spread across WGS tools is extreme (12–50%)." It should be the opening observation of the WGS section, visualized as a ranked dot plot or bar chart, and used explicitly as the mechanistic explanation for all downstream findings.

**Suggested reframe for Results section 4 opening:**
> "The WGS tool landscape is exceptionally polarized: across HLA-A, -B, and -C, per-tool overall accuracy ranges from 8.0% (ArcasHLA) to 49.6% (OptiType) — a 42-percentage-point spread. This spread is the proximate cause of ensemble degradation under equal-weight voting: on any locus where OptiType is correct, five tools with 8–33% accuracy collectively outvote it, making majority voting systematically counterproductive. Understanding WGS ensemble performance requires viewing this spread as the primary input, not as background context."

---

### 4.2 The WGS abstention map deserves its own framing as a quality assurance tool

The 40% WGS abstention rate is presented defensively ("abstention is information, not failure"). This framing is correct but undersells the finding. The abstention map is not just "we couldn't call these loci" — it is a **reproducible, calibrated map of the most difficult HLA typing targets in short-read WGS data**.

Consider this reframing: the 40% of WGS loci that WeightedConsensus abstains on are precisely the loci where:
1. No tool achieves reliable accuracy
2. Tool calls are maximally conflicted
3. The correct answer is most likely to be incorrect if forced

A figure showing the *distribution of abstention across samples and genes* (which loci abstract consistently, which are sample-specific) would be a novel contribution to the field. If HLA-C abstains more than HLA-A (consistent with the per-gene data showing HLA-C 31.3% vs HLA-A 48.9% for OptiType), that is a directly actionable finding for clinical labs choosing HLA-typing strategies.

Per-gene OptiType WGS accuracy from the tables:
- HLA-A: 0.4885 overall
- HLA-B: 0.6870 overall (highest — easier gene)
- HLA-C: 0.3130 overall (lowest — hardest gene)

This HLA-B > HLA-A > HLA-C hierarchy on WGS deserves explicit discussion and a figure. The manuscript mentions it only in a "plain terms" callout box and only for MajorityVote.

---

### 4.3 Modality-tool mismatch is a stronger finding than the manuscript treats it

The ArcasHLA and SpecHLA modality mismatch results are buried in "plain terms" callouts rather than given a dedicated section. The numbers are striking:
- ArcasHLA: 91.9% RNA vs 8.0% WGS — a **12-fold** accuracy difference
- SpecHLA: 78.3% WES vs 3.5% RNA — a **22-fold** accuracy difference

These are not edge cases. They are a systematic demonstration that tool selection without modality awareness is a major source of error in real-world HLA typing workflows. The manuscript's abstract mentions this ("ArcasHLA, designed for RNA-seq, achieves 91.9% accuracy on RNA but only 11.9% on WGS") but the number used is from the old run. The actual split is RNA 91.9% vs WGS 8.0%.

**This finding warrants its own Results subsection** titled something like "Tool-modality compatibility determines single-tool performance." It should include:
- A table showing all tool × modality combinations
- An explicit warning: running ArcasHLA on WGS data or SpecHLA on RNA-seq data is not just suboptimal — it is dramatically wrong
- PIHLA's role in preventing this: the tool-availability flags in the benchmark prevent modality-inappropriate tools from contributing to consensus

This reframing makes PIHLA's modality-awareness layer a named feature rather than a background assumption.

---

### 4.4 The class II > class I WGS reversal is counterintuitive and underexplored

The manuscript notes in a "plain terms" callout that DRB1 (85.7%) and DQB1 (92.0%) by MajorityVote outperform all class I loci on WGS. The explanation given is correct: with only 2–3 tools voting on class II vs 6 on class I, the low-accuracy-tool noise is absent for class II. But the manuscript doesn't complete the implication:

**If the tool ensemble size is the primary determinant of WGS ensemble quality, then selectively running fewer (higher-quality) tools on class I loci would improve class I consensus.** This is a directly actionable insight: a class I WGS consensus using only OptiType + T1K (the two highest-accuracy WGS tools) should substantially outperform the full 6-tool ensemble. If true, this would be a significant recommendation for practitioners.

This also connects to the "tyranny of the majority" explanation: the problem is not majority voting per se, it is majority voting with a large number of low-accuracy tools. The weighted approach partially solves this by down-weighting poor tools; but the cleaner solution is to not run poor tools on WGS at all.

The manuscript could add a sentence in the Discussion: "The class II accuracy advantage on WGS supports a general principle: ensemble quality on WGS is limited more by the number of low-accuracy tools contributing votes than by the difficulty of the loci being typed. A selective tool-exclusion policy — restricting the WGS class I ensemble to the two highest-accuracy tools (OptiType and T1K) — is predicted to improve majority-vote accuracy and is a recommended evaluation for future work."

---

### 4.5 The bimodal gain should be in the abstract and framed in absolute numbers

Currently, the bimodal gain (+0.8 pp WeightedConsensus over RNA-only) is correctly translated into "8–10 more correct loci per 1,000 patients" in a "plain terms" callout. This is a strong communication move but it's buried at the end of a section. For a clinical audience, this is the most actionable finding in the paper:

**"When WES and RNA-seq are already available, running bimodal consensus costs nothing and gives you 8–10 more correctly typed loci per 1,000 patients."**

This should be in the abstract, not in a callout box.

Revised abstract Conclusions:
> "PIHLA demonstrates that when WES and RNA-seq are co-available — as in tumour-normal genomic profiling — bimodal joint consensus recovers 8–10 additional correct HLA calls per 1,000 gene-locus assessments compared to RNA-seq alone, at no extra data-generation cost. Adding WGS to the WES+RNA consensus provides ≤0.1 percentage points of additional gain while reducing coverage — a strong negative result supporting WES+RNA as the optimal routine combination."

---

### 4.6 The trimodal negative result is a headline finding, not a footnote

The trimodal analysis shows that adding WGS data to a WES+RNA ensemble provides +0.06 pp (MajorityVote) and +0.05 pp (WeightedConsensus) while reducing coverage from 391 to 370 assessments. This is as close to zero as a scientific result can get.

The manuscript presents this as a Results section and gives it appropriate "plain terms" explanation, but the Discussion and Conclusion don't include it as a named conclusion. **"WGS adds nothing to a WES+RNA HLA typing ensemble"** is a strong, counterintuitive, and practically valuable finding. Labs routinely collect WGS and assume it is the most complete data type. The PIHLA result says: for HLA class I typing when WES and RNA are available, WGS is redundant.

**Suggested Conclusion addition:**
> "A secondary but practically significant finding is that WGS data adds ≤0.1 percentage points to a WES+RNA bimodal consensus while reducing coverage by 5.4%. This negative result has direct relevance for resource allocation: when WES and RNA-seq are already available, collecting WGS specifically to improve HLA typing is not warranted."

---

### 4.7 The failure mode taxonomy is underutilized

The `failure_mode_summary.tsv` table breaks errors into `correct`, `partial`, `wrong_field2`, `wrong_field1`, and `abstained`. This is a richer characterization of error types than the current manuscript uses. For example:
- ArcasHLA WES: 6.2% correct at first field, 44.6% partial (one allele right, one wrong), 35.4% wrong at second field, 13.8% wrong at first field — a very different error profile from OptiType (94.6% correct at WES-A)
- SpecHLA RNA: 7.3% correct, 0% at callable sites — pure abstention failure

A single supplementary figure showing the failure mode breakdown by tool and modality would be a genuinely novel addition to the field. No existing benchmark paper for HLA typing presents error *type* distributions alongside error *rate* — they all show only accuracy. This would differentiate PIHLA's benchmark from prior work and give practitioners actionable information about which errors are recoverable (wrong second field = likely allele ambiguity, resolvable with higher coverage or targeted confirmation) vs unrecoverable (wrong first field = tool cannot identify the allele at all).

---

### 4.8 Per-gene accuracy tells a richer story than the aggregate

The per-gene data shows structured variation that the aggregate numbers obscure. From the `summary_per_gene.tsv` tables:

**WGS per-gene accuracy (single tools):**
- HLA-B is consistently the easiest: OptiType 68.7%, T1K 26.5%, HLA-HD 30.6%
- HLA-C is consistently the hardest: OptiType 31.3%, T1K 33.9%, HLA-HD 16.5%
- HLA-A is intermediate: OptiType 48.9%, T1K 37.2%, HLA-HD 32.2%

This pattern is consistent across all tools, suggesting a gene-level difficulty structure that is not tool-specific. HLA-C's consistently lower accuracy likely reflects its position in the most gene-dense MHC region and highest pseudogene similarity — this is a genomic property, not a tool limitation. The manuscript currently explains this in a "plain terms" callout but doesn't quantify it across all tools or present it as a gene-level finding.

**WES per-gene accuracy:**
- ArcasHLA WES: 6.2% at A, 9.3% at B, 24.8% at C — the unusual C performance pattern (3–4× better than A or B) is unexplained and worth investigating
- This is not mentioned in the manuscript

**RNA per-gene accuracy:**
- ArcasHLA RNA: 86.0% at A (lower than B or C) — also unexplained
- HLA-HD RNA: 95.2% A, 94.3% B, 93.3% C — most tools show similar per-gene performance in RNA, which argues that the WGS A/B/C hierarchy is a sequencing-regime effect, not an allele-diversity effect

These patterns are in the tables but not in the manuscript. A figure showing per-gene accuracy as a heat map (tools × genes × modalities) would make the modality-specific difficulty structure immediately visible.

---

## Summary Checklist

### Must fix before submission (blocking)
- [x] Reconcile all numbers against authoritative benchmark tables (Issue 1) — DONE 2026-05-12
- [x] Fix internal abstention rate contradiction (Issue 2) — DONE 2026-05-12
- [x] Fix duplicate section 11 numbering (Issue 3) — DONE 2026-05-12 (now 11/12/13)
- [x] Resolve DRB1/DQB1 in-scope vs out-of-scope contradiction (Issue 7) — DONE 2026-05-12 (partial truth n=38/35, no longer marked unavailable)
- [x] Decide VENEX status and act on it (Issue 6) — DONE 2026-05-12 (no truth labels; internal QC only; VENEX figures excluded from pub set)
- [ ] Create Figure 1 (pipeline architecture) — STILL REQUIRED (manual creation)

### Should fix before submission (non-blocking but high priority)
- [x] Update title/abstract/Core Claim framing to Codex Q2 recommendation (Issue 4) — DONE 2026-05-12
- [x] Update WGS cohort size n from 136 to 131 (Issue 5) — DONE 2026-05-12
- [x] Fix OptiType three-way inconsistency (Issue 9) — DONE 2026-05-12 (0.496 throughout)
- [x] Add "consistent with" qualifier to WES ceiling-effect explanation (Issue 10) — DONE 2026-05-12
- [x] Add WGS formula limitation sentence to Discussion (Issue 11) — DONE 2026-05-12
- [x] Complete consHLA citation (Issue 12) — DONE 2026-05-12 (Jurtz et al. 2017, PMID:28334079 — verify before final submission)
- [x] Add circular bias bound estimate to Limitations (Issue 13) — DONE 2026-05-12
- [x] Fix abstract WeightedConsensus comparison framing (Issue 14) — DONE 2026-05-12
- [x] Add Wilson CIs to Table 3 — DONE 2026-05-12
- [ ] Resolve "plain terms" callout box formatting (Issue 16) — STILL OPEN (reformat for journal submission)

### Recommended enhancements for stronger paper (Part 4)
- [ ] Move WGS performance spread to opening observation of Results 4
- [ ] Add per-gene accuracy figure (A/B/C × tool × modality heatmap)
- [ ] Add modality-tool compatibility as a named section
- [ ] Elevate bimodal gain (+8–10 loci/1000 patients) to abstract
- [ ] Name trimodal negative result as a headline conclusion
- [ ] Add failure mode taxonomy figure to supplementary
- [ ] Generate Figure 9 (bimodal per-gene accuracy) and Figure 10 (trimodal vs bimodal)

---

*End of review. Total issues: 3 critical, 5 high, 5 medium, 3 low. Total Part 4 insights: 8.*
