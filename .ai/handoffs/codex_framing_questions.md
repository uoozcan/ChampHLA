---
handoff_id: pihla-framing-2026-04-18
from: claude
to: codex
task_type: decision-support
priority: high
created: 2026-04-18
updated: 2026-04-19
status: open
---

# Codex: PIHLA Manuscript Framing Questions

This handoff requests Codex evaluation on five high-stakes scientific and editorial
decisions before manuscript revision proceeds. These questions involve tradeoffs where
a second-opinion evaluation is more valuable than direct implementation.

Claude owns the implementation follow-up once Codex recommendations are received.

---

## Context

PIHLA is a calibrated, uncertainty-aware HLA ensemble framework benchmarked on:
- WGS wave1: 42 1000G samples, 25/8/9 train/val/holdout, loci A/B/C, IMGT 3.59.0
- Full tri-modal: WES holdout n=12, RNA holdout n=7 (separate benchmark)

Confidence integration is active. Tools with confidence parsers (OptiType, ArcasHLA,
HLA-HD, T1K, Kourami) contribute real confidence data. Runtime guardrails
(poor_calibration: ECE > 0.35 or Brier > 0.35) block confidence from amplifying
weights for miscalibrated tools. In the WGS wave, HLA-HD, Kourami, and OptiType all
triggered poor_calibration and fell back to base_reliability weights.

ArcasHLA WGS: 0% accuracy → 98.7% shrinkage → final weight 0.0101 (guardrail example).
T1K retained a partial confidence boost (guardrail_factor 0.70, guardrail_status applied).

WGS holdout results (n=9):
- OptiType: 0.5185 overall correct-call rate
- WeightedConsensus: 0.5185 overall, 0.9259 callable rate
- MajorityVote: 0.4815 overall, 0.7778 callable rate

WES holdout results (n=12):
- MajorityVote: 0.9444 overall
- WeightedConsensus: 0.8611 overall (underperforms)

---

## Question 1: WES WeightedConsensus Underperformance Mechanism

WeightedConsensus (0.8611) underperforms MajorityVote (0.9444) in WES holdout (n=12).
The manuscript attributes this to ceiling effects (all single tools already high accuracy),
but no ablation supports this interpretation.

**Q**: Is this underperformance best explained as:
(a) a structural ceiling effect — when tools already agree at high accuracy, any weighting
    that reduces one tool's vote can only hurt?
(b) a weight calibration artifact — the training set weights are not well-matched to the
    holdout distribution?
(c) something else?

What experiment would cleanly distinguish (a) from (b)? Specifically: would running
WeightedConsensus with uniform weights (1/N per tool) on the same WES holdout distinguish
structural from calibration effects? Or does it need confidence ablation (confidence
removed entirely, weights = base_reliability only)?

**Stakes**: If (a), the result is a feature not a bug and should be framed as "calibrated
weighting is most beneficial in low-signal settings (WGS)." If (b), the training
procedure needs revision before submission.

---

## Question 2: Framing "Confidence-Calibrated" Given Guardrail Dominance

The manuscript's working title includes "calibrated uncertainty-aware ensemble." The core
confidence claim is that benchmark-derived guardrails prevent miscalibrated tools from
distorting consensus.

However, in the WGS wave, 3 of 5 tools with confidence data triggered poor_calibration
and fell back to base_reliability. The weighted consensus is therefore primarily
reliability-weighted, not confidence-weighted, in the current benchmark.

**Q**: Given that guardrails fire on the majority of tools in the current WGS benchmark,
is "confidence-calibrated ensemble" still the accurate core claim, or should the framing
shift? Two candidate alternatives:

(a) Keep "confidence-calibrated" — the guardrail behavior IS the calibration story.
    The paper shows that naive confidence would harm the ensemble; the guardrail system
    prevents this. The negative result (confidence blocked) is a positive finding about
    the system's integrity.

(b) Shift to "reliability-weighted ensemble with confidence guardrails" — a more accurate
    description that avoids implying confidence is actively boosting performance when it is
    mostly being blocked.

**Stakes**: This affects the abstract, title, and methods framing. Option (a) is stronger
scientifically but requires the paper to be clear that "confidence-calibrated" means
"calibration-gated" not "confidence-boosted." Option (b) is more conservative but
potentially less compelling as a methodological contribution.

---

## Question 3: Statistical Defensibility of Small Holdout Cohorts

WGS holdout: n=9 gene-pairs per locus (27 total gene rows across A, B, C).
WES holdout: n=12.
RNA holdout: n=7.

The plan is to add Wilson score 95% CIs to all key performance numbers.

**Q**: For a methodology paper at a top bioinformatics venue (Bioinformatics, PLOS Comp Bio,
Genome Biology), are these sample sizes defensible? Specifically:

(a) Will reviewers accept Wilson score intervals as sufficient statistical reporting
    for a framework paper (i.e., the paper is selling the method, not claiming large
    effect sizes)?

(b) Are there precedents in HLA typing benchmark papers for similarly small holdout
    cohorts? (e.g., OptiType paper, T1K paper, ArcasHLA paper — what were their
    holdout sizes?)

(c) Is there a minimum holdout size threshold below which reviewers are likely to
    reject on statistical grounds alone?

**Stakes**: If holdout sizes are insufficient, the 50-sample WGS trimodal run (gate
just cleared) needs to complete and generate new holdout numbers before submission.

---

## Question 4: Publication Scope With A/B/C Loci Only

Current benchmark evaluates only HLA-A, HLA-B, HLA-C. The Gourraud et al. 2014 truth
does not include DRB1 or DQB1. These are implemented in the framework but truth-blocked.

**Q**: At which venue tier does A/B/C-only evaluation become a rejection risk?

(a) Is this acceptable at Bioinformatics or PLOS Comp Bio as a framework/methods paper
    where the contribution is the ensemble infrastructure, not the locus coverage?

(b) Would Genome Biology or Nature Methods reviewers expect DRB1/DQB1 to be included
    before accepting?

(c) Is there a complementary truth resource (e.g., IPD-IMGT/HLA reference panel,
    another 1000G effort, DKMS cohort) that could provide DRB1/DQB1 truth for even
    a subset of samples to demonstrate the framework's 5-locus capability?

**Stakes**: This determines target venue and whether a DRB1/DQB1 truth acquisition
step is needed before submission.

---

## Question 5: Weight Formula (0.7/0.3) Sensitivity Requirement

The current weight formula is `final_weight = 0.7 × base_reliability + 0.3 × effective_confidence`.
This split was chosen empirically. No ablation or sensitivity analysis has been run.

**Q**:
(a) Will reviewers at top venues require a sensitivity analysis (e.g., showing results
    at 1.0/0.0, 0.5/0.5, 0.0/1.0) to establish that the 0.7/0.3 choice is not
    cherry-picked?

(b) Is there a published precedent for this type of reliability-confidence blending
    in ensemble genomics tools that could provide a prior justification for the split?

(c) If a sensitivity analysis is required, is the correct framing:
    — Show that 0.7/0.3 outperforms alternatives (empirical justification)
    — Or show that results are robust across the range (robustness justification)?

**Stakes**: If a sensitivity analysis is required, it needs to be added to `run_1000g_benchmark.py`
as a weight-variant run mode before submission.

---

---

## Question 6: Weight Formula Sub-Optimality for WGS — Correct Framing (ADDED 2026-04-19)

Weight sensitivity analysis (α=1.0/β=0.0, α=0.5/β=0.5, α=0.0/β=1.0) has now been
completed on the 100-sample tri-modal cohort. Results from
`analysis/weight_sensitivity/sensitivity_comparison.tsv`:

| Weight variant | WES WC | WGS WC | RNA WC |
|---|---|---|---|
| 0.7/0.3 (default) | 0.9167 | 0.2667 | 1.0000 |
| 1.0/0.0 (reliability-only) | 0.8627 | 0.2963 | 0.8867 |
| 0.5/0.5 (equal weight) | 0.8693 | 0.3030 | 0.8867 |
| 0.0/1.0 (confidence-only) | 0.8693 | 0.3165 | 0.8867 |

Key finding: 0.7/0.3 is optimal for WES (beats all variants) but is beaten by ALL three
sensitivity variants for WGS (0.296–0.317 vs 0.267). Root cause: 4 of 5 WGS tools
trigger `poor_calibration` guardrail, falling back to base_reliability regardless of β.
ArcasHLA WGS has near-zero effective confidence (0.0068); higher β reduces its
contribution slightly, which is why confidence-only (β=1.0) gives the best WGS result.
RNA drops from 1.0 (default) to 0.887 (all variants) — confidence scores from some
tools slightly hurt the RNA ensemble.

**Q**:
(a) How should this result be framed in the manuscript? The 0.7/0.3 split was not chosen
    to optimise WGS performance — it was chosen for WES/RNA stability. The WGS result
    is dominated by guardrail behaviour. Is "the formula was selected to maximise WES
    stability; WGS performance is guardrail-dominated" a defensible framing?

(b) Does the WGS sub-optimality of 0.7/0.3 require a change to the default formula
    before submission, or is documenting the sensitivity analysis in supplementary
    sufficient?

(c) Given that RNA drops from 1.0 to 0.887 under any sensitivity variant, does the
    default formula's RNA advantage reflect a real weighting benefit, or an artifact
    of RNA tool agreement being so high that any deviation from reliability-only weights
    can reduce it?

**Stakes**: If (b) requires a formula change, the trimodal benchmark and wave1 must be
re-run before submission. If the sensitivity analysis is supplementary-only, the
manuscript framing needs to clearly acknowledge the WGS sub-optimality without
undermining the overall claim.

---

## Requested Output

For each question:
1. Primary recommendation (one of the candidate options, or a third alternative)
2. Reasoning (scientific, editorial, or precedent-based)
3. Confidence level (high / medium / low)
4. Any suggested follow-up action for Claude to implement

Return your response as a new file in `.ai/handoffs/codex_response_framing_2026-04-18.md`.

---
