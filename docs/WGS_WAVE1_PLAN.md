# WGS Wave 1 Benchmark Plan

## Scope
- Real-data 1000 Genomes WGS-only expansion using the on-disk result pool.
- Truth source: Gourraud et al. 2014 1000 Genomes HLA calls from `/scratch/project_2008084/ozcanumu/hla_calibration/conf/ground_truth_data.csv`.
- Headline loci remain `HLA-A`, `HLA-B`, and `HLA-C`.
- Required tool set for Wave 1: `OptiType`, `T1K`, `ArcasHLA`, `SpecHLA`, `HLA-HD`, `Kourami`.

## Cohort Size
- Truth-backed WGS samples with the full six-tool result set on disk: `42`.
- This is the first larger cohort that avoids the three-sample pilot split instability.

## Population Composition
- `CEU`: `7` total, `4` training, `1` validation, `2` holdout.
- `CHB`: `8` total, `5` training, `2` validation, `1` holdout.
- `GBR`: `7` total, `4` training, `1` validation, `2` holdout.
- `TSI`: `9` total, `5` training, `2` validation, `2` holdout.
- `YRI`: `11` total, `7` training, `2` validation, `2` holdout.

## Confidence Coverage Status
- Native confidence is now available in the pilot benchmark for `OptiType`, `T1K`, and `ArcasHLA`.
- `OptiType` uses native `Objective` values from `*_result.tsv`.
- `T1K` uses genotype-table read support aggregated from allele support columns.
- `ArcasHLA` uses `*.genes.json` support plus raw score values where available.
- Wave 1 should therefore report both accuracy and confidence-calibration outputs for at least these three WGS tools.

## Recommended Execution Order
1. Generate a WGS-only benchmark config from the existing phase-gated config, limiting runs to the six WGS tools and the Wave 1 manifest.
2. Benchmark all 42 samples with sample-level grouped splits from `conf/1000g_wgs_wave1_manifest.tsv`.
3. Report single-tool, majority-vote, and weighted-consensus performance at 2-field resolution.
4. Report per-gene gains for `A`, `B`, and `C`, plus confidence calibration for `OptiType`, `T1K`, and `ArcasHLA`.
5. Use the larger WGS cohort to stabilize runtime weights before adding more WES/RNA samples.

## Immediate Deliverables
- `conf/1000g_wgs_wave1_manifest.tsv`: frozen deterministic split manifest for the 42-sample WGS wave.
- A future WGS-only benchmark config should consume this manifest directly to avoid ad hoc sample selection.
