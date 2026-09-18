# WGS Interim Notes

These notes summarize provisional observations from already-available WGS benchmark outputs while refreshed WES and RNA jobs are still queued. They are not final tri-modal scientific claims.

## Data Basis
- Existing WGS tool outputs under `/scratch/project_2008084/hla_calibration/wgs/results`
- Current benchmark tables under `analysis/latest_benchmark/tables/`
- Current headline focus limited to `HLA-A`, `HLA-B`, and `HLA-C`

## Provisional Observations
### 1. Current WGS method comparison structure is usable
The existing benchmark tables already contain WGS rows for single-tool, majority-vote, and weighted-consensus comparison. This means the reporting structure for the final manuscript is in place even though the final tri-modal refresh is still pending.

### 2. Weighted consensus is distinguishable from majority vote on WGS
In the currently generated WGS benchmark rows, `WeightedConsensus` reaches `overall_correct_call_rate = 1.0` with `callable_rate = 1.0`, while `MajorityVote` shows `overall_correct_call_rate = 0.8333` with `callable_rate = 0.8333`. This is a useful signal that the calibrated ensemble can behave differently from equal-weight voting on real WGS outputs.

### 3. The clearest provisional gain is at HLA-C
The current `per_gene_gain.tsv` indicates no gain over majority vote for WGS `A` or `B`, but a positive gain at `C` (`gain_vs_majority = 0.5`). In the same provisional table, the best single tool for WGS `C` is `SpecHLA` at `1.0`, while weighted consensus also reaches `1.0`. This suggests that disagreement resolution at `HLA-C` may be where calibration is most informative in the present subset.

### 4. SpecHLA is the strongest current WGS single-tool signal in the available subset
The existing WGS method rows show `SpecHLA` at `overall_correct_call_rate = 1.0`, while `OptiType` is at `0.8333`. The current weight table also preserves `SpecHLA` at `final_weight = 1.0`, with no confidence coverage but full fallback to reliability. This is consistent with the current weighting behavior and supports describing SpecHLA as a strong WGS contributor in the provisional subset.

### 5. Confidence weighting is operational on WGS
The current WGS confidence-weight table distinguishes tools numerically: `OptiType` has `final_weight = 0.8788`, while `SpecHLA` remains `1.0`. This means the calibrated runtime pathway is not just scaffolding; it already produces differentiated weights on real WGS outputs.

## Interpretation Boundaries
- These notes are provisional and should not be used as final tri-modal claims.
- The current benchmark metadata in `analysis/latest_benchmark/tables/benchmark_metadata.json` is incomplete for final reporting and should be replaced by the queued phase-gated refresh.
- WES and RNA may change the final comparative framing, especially around modality-specific calibration behavior and abstention patterns.

## Figure/Table Readiness
The current files already support the structure of:
- `Figure 2` single-tool vs majority-vote vs weighted-consensus comparison
- `Figure 3` per-gene gains
- `Figure 4` confidence calibration
- `Figure 5` abstention tradeoff
- `Figure 6` discordance taxonomy
- `Figure 7` learned confidence weights

The remaining task is not figure design but final real-data refresh and number replacement.
