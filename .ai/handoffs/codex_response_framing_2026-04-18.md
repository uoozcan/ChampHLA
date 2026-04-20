# Codex Response: PIHLA Framing Questions

Source handoff: `pihla-publish/.ai/handoffs/codex_framing_questions.md`
Responder: `Codex`
Role: `Scientific/editorial second-opinion only`
Status: `complete`

## Question 1
Primary recommendation:
The current WES underperformance is best framed as a structural ceiling effect, but only provisionally because the existing evidence is ablation-limited.
Reasoning:
The local benchmark tables support the direction of the current interpretation, but not its strongest version. In the WES holdout of the full cohort benchmark, `MajorityVote` reaches `0.9444` overall while `WeightedConsensus` reaches `0.8611`, despite the best single tools already being very strong (`OptiType` and `POLYSOLVER` at `0.9167`). In the 100-sample trimodal sensitivity run, WES `WeightedConsensus` stays high but below `MajorityVote` across all tested formulas (`0.8627` to `0.9167` depending on the setting), which argues against a narrow `0.7/0.3` calibration artefact and toward a regime where weighting cannot add much once the strongest tools already agree. That said, the manuscript should not claim this mechanism is proven, because there is still no direct ablation showing whether the loss comes from weighting in general versus from the learned reliability profile. The cleanest discriminating experiment is to run uniform-weight consensus on the same WES holdout. If uniform weighting still underperforms majority vote, the result is more plausibly structural; if it does not, the learned weighting scheme is the more likely culprit. A second-stage base-reliability-only ablation is then useful to separate learned-weight effects from confidence-gating effects, but it is not the first discriminating check.
Confidence:
medium
Suggested follow-up for Claude:
Describe the current WES result as "consistent with a structural ceiling effect" rather than as a settled mechanism, and add a sentence that uniform-weight consensus on the same WES holdout would be the cleanest next ablation if the team wants to convert this from a plausible interpretation into a demonstrated claim.

## Question 2
Primary recommendation:
Shift the main framing to `reliability-weighted ensemble with confidence guardrails`, and treat `confidence-calibrated` as a stronger but riskier secondary description only if it is explicitly defined as calibration-gated rather than confidence-boosted.
Reasoning:
The current repo evidence is clear that guardrail dominance, not broad confidence amplification, is the primary story in the WGS wave. In the authoritative runtime table for `benchmark_wgs_wave1`, `WeightedConsensus` reaches `0.5185` overall correct-call rate and `0.9259` callable rate, while the benchmark metadata shows that `OptiType`, `HLA-HD`, and `Kourami` all fail the calibration guardrail and fall back to reliability-only behavior. The manuscript and discussion are already moving in the right direction: the system’s contribution is that it prevents badly calibrated confidence from distorting runtime voting. That is a valuable methodological result, but it is more accurately described as reliability-weighting constrained by calibration checks than as confidence actively improving consensus in the current benchmark. Keeping `confidence-calibrated` as the lead label would be defensible only if the paper repeatedly clarifies that the calibration story is mostly about blocking unsafe confidence rather than boosting well-calibrated tools. Editorially, the safer and more precise headline is therefore the more conservative option.
Confidence:
high
Suggested follow-up for Claude:
Keep the current conservative placeholder framing in the title, abstract, discussion, and conclusion. If Claude wants to preserve some form of `confidence-calibrated` language, demote it to supporting text and define it explicitly as "calibration-gated confidence integration" rather than implying broad confidence-driven performance gains.

## Question 3
Primary recommendation:
The current holdout sizes are defensible for a conservative methods/framework paper if Wilson intervals and limitations are explicit, but they are not strong enough for aggressive effect-size framing or a high-ambition general-performance narrative.
Reasoning:
The main distinction here is not whether the sample sizes are ideal, but whether the manuscript overclaims what they can support. The current local manuscript posture is already conservative: WGS wave1 has a 9-sample holdout, WES has 12, and RNA has 7, with Wilson intervals reported or planned. That is small, but not automatically disqualifying for a methods paper if the paper’s value is framed around infrastructure, transparency, guardrails, and reproducibility rather than around large definitive effect sizes. This is especially true because benchmark precedent in this domain is heterogeneous rather than uniformly large. OptiType reported a much larger pooled benchmark (`361` benchmark runs) and arcasHLA evaluated `358` RNA-seq 1000 Genomes samples, so those papers set a stronger empirical scale than PIHLA currently has. T1K also reports a substantially larger 1KGP-based comparison (`463` samples with HLA genotype, RNA-seq, and WES). So the current PIHLA holdouts are publishable only under a more cautious methods-paper argument, not because they are competitive in scale with the strongest tool papers. I do not think there is a universal minimum threshold below which reviewers reject on sample size alone, but there is a real risk threshold below which reviewers reject if the manuscript uses small cohorts to make sweeping comparative claims. Wilson intervals are necessary, but by themselves they are not enough unless the claims are correspondingly modest.
Confidence:
medium
Suggested follow-up for Claude:
Keep Wilson intervals in all headline tables, keep the WGS wave as the primary quantitative claim, and explicitly position WES/RNA as secondary evidence. Do not imply that the present cohort sizes support strong cross-modality superiority claims.

## Question 4
Primary recommendation:
A/B/C-only scope is likely acceptable for a conservative Bioinformatics or PLOS Computational Biology methods framing, but it becomes a meaningful rejection risk at higher-tier venues that are likely to expect broader locus coverage or stronger validation of the framework’s full intended scope.
Reasoning:
The local manuscript is already being steered toward the right tier. Bioinformatics scope guidelines emphasize realistic biological scenarios, public methods, and comparison on standard datasets, which aligns well with the current PIHLA emphasis on reproducibility, truth-backed benchmarking, and explicit scope control. PLOS Computational Biology’s benchmarking and methods guidance also emphasizes transparency, reproducibility, and clear validation rather than demanding any one fixed locus set. Under that posture, A/B/C-only reporting can be defended if the paper clearly says that DRB1/DQB1 are implemented in the framework but not yet in the headline truth-backed evaluation. By contrast, a higher-tier venue that expects the framework contribution to be demonstrated at the breadth of its advertised capability is more likely to view A/B/C-only reporting as a limitation that needs to be closed before acceptance. There is also a useful opening here: external benchmark precedent indicates that broader 1000 Genomes-derived truth resources do exist. The arcasHLA paper describes 1000 Genomes individuals with `HLA-A`, `-B`, `-C`, `-DRB1`, and `-DQB1` typing, suggesting that complementary truth acquisition may be possible even if it is not available in the current PIHLA matched benchmark slice. That means DRB1/DQB1 truth expansion is value-adding and strategically important, but not necessarily mandatory for the conservative venue posture.
Confidence:
medium
Suggested follow-up for Claude:
Position the manuscript explicitly as a conservative methods/framework paper for Bioinformatics/PLOS-style venues, and add a forward-looking note that DRB1/DQB1 truth expansion is a high-value extension rather than an already-satisfied validation target.

## Question 5
Primary recommendation:
Yes, reviewers are likely to expect a sensitivity analysis for the `0.7/0.3` split, but the right justification is robustness and interpretability, not a claim that the chosen split is globally optimal; in practice, that expectation is now already satisfied by the sensitivity results in Question 6.
Reasoning:
An empirical blend between reliability and confidence is exactly the sort of design choice that reviewers may suspect is cherry-picked unless its behavior across nearby settings is shown. The current repo state is already ahead of the question as originally posed: the sensitivity analysis now exists, and it shows a nuanced result rather than a uniformly dominant optimum. That is actually better for credibility. It supports the argument that the default was chosen as a conservative engineering tradeoff rather than as an overfit global optimum. I did not find strong HLA-specific precedent for an identical reliability-confidence blend that would let PIHLA lean on prior literature as a direct parameter justification. That means the manuscript should not oversell literature support for the exact `0.7/0.3` ratio. Instead, it should say that sensitivity was examined and that the default is justified by the stability profile the authors wanted in the current benchmark regime.
Confidence:
high
Suggested follow-up for Claude:
Remove any suggestion that `0.7/0.3` was self-evidently correct a priori. Present the sensitivity analysis as a robustness/interpretability check that supports the conservative default rather than as proof that the ratio is universally optimal.

## Question 6
Primary recommendation:
Frame `0.7/0.3` as a conservative default chosen for WES/RNA stability in the current benchmark regime, acknowledge that WGS behavior is guardrail-dominated and sub-optimal under that default in the trimodal sensitivity run, and treat this as a limitation that can be documented without forcing a pre-submission formula change.
Reasoning:
The runtime evidence supports a conservative rather than crisis-level interpretation. In the `weight_sensitivity` table, the default is best for WES (`0.9167`) and RNA (`1.0`) but worst for WGS (`0.2667`) relative to the three alternatives (`0.2963`, `0.3030`, `0.3165`). The manuscript’s own sensitivity section already points to the key mechanistic interpretation: most WGS tools in the trimodal benchmark fail calibration guardrails, so the WGS result is largely governed by fallback reliability behavior plus the residual effect of a few tools rather than by a general confidence-boosting principle. That means the manuscript should not claim the default is best for WGS. But it also does not mean the formula must be changed before submission, because the present paper is not strongest when cast as "we found the single globally optimal blend." It is strongest when cast as "we built a guarded ensemble framework and examined how the default behaves across modalities." The RNA pattern should also be described cautiously: the default’s RNA advantage is likely real in the current benchmark regime, but because all non-default variants collapse to the same lower value, the safer interpretation is that RNA is a high-agreement setting in which perturbing the default weighting can reduce consensus quality. That is a useful empirical finding, not a universal law.
Confidence:
high
Suggested follow-up for Claude:
Keep the default formula for this manuscript cycle, add a concise main-text acknowledgement that WGS trimodal performance is sub-optimal under the default and likely guardrail-dominated, and push the detailed sensitivity table into supplementary or results text rather than rewriting the whole method around WGS-specific optimization.

## Editorial Bottom Line
- Safest framing now:
  A benchmark-derived, reliability-weighted HLA ensemble with empirical confidence guardrails that protects runtime voting from miscalibrated tool confidence while improving decision coverage in the strongest current WGS benchmark.
- Stronger but riskier framing:
  A confidence-calibrated HLA ensemble, but only if the manuscript repeatedly clarifies that the key positive result is calibration-gated confidence blocking rather than broad confidence-driven performance gains.
- Mandatory caveats:
  Small WES and RNA holdouts; A/B/C-only headline scope; incomplete RNA tool coverage for Seq2HLA and SpecHLA; WES underperformance is only provisionally structural pending direct ablation; WGS trimodal sensitivity shows the default formula is not globally optimal; `ANALYSIS_2026-04-01.md` appears to contain a WGS metric typo because the authoritative runtime table gives `WeightedConsensus` `0.5185` overall and `0.9259` callable in `benchmark_wgs_wave1`.

## Blocking vs Non-Blocking
- Blocking:
  Final high-level terminology for the main contribution (`confidence-calibrated` versus `reliability-weighted with confidence guardrails`); final venue-facing positioning if the authors still want to target above a conservative methods-paper tier.
- Non-blocking:
  The exact WES mechanism can remain caveated rather than fully resolved; the `0.7/0.3` sensitivity analysis can be handled as a documented limitation/robustness result; DRB1/DQB1 truth expansion can remain future work for the current conservative submission posture; the WGS discrepancy in `ANALYSIS_2026-04-01.md` can be corrected editorially without changing the scientific conclusions.
