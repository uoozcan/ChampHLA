# ChampHLA publication recovery status

## Completed locally

- Immutable source DOCX copied and checksummed.
- Nine manuscript problems registered; four are critical.
- Source DOCX audit fails as intended: eight verification markers and invalid
  historical WGS performance remain.
- Benchmark manuscript copy passes the automated claim audit.
- Conditional method manuscript is blocked by a pending-confirmation gate.
- Canonical result registry created; invalid WGS rows cannot be abstract-eligible.
- Cross-project exposure ledger covers 138 unique previously exposed 1000G
  subjects across ten training, fold, pilot, prediction, and evaluation sources.
- Official WGS/WES/RNA assay manifest created without reading HLA genotypes.
- The 2014 laboratory truth parser never selects an arbitrary representative
  from a cross-two-field ambiguity.
- Frozen initial WGS (120) and WES (89) rosters are feasible.
- Strict RNA feasibility is 99 unexposed donors, below the locked minimum of
  130; RNA three-modality confirmation is therefore currently infeasible.
- Core release-candidate tests pass. RefFormer is isolated as a negative
  ablation with a separate pinned CPU test environment.

## Roihu gates

1. Sync this Git-controlled release candidate.
2. Run the environment/reference/container preflight.
3. Complete manual review of the two-subject full-CRAM WGS pilot.
4. Require five nonempty caller-native outputs and 100% parser agreement.
5. Run the frozen WGS/WES rosters in rolling batches without truth.
6. Freeze predictions, then join sealed truth exactly once.

No new WGS accuracy or three-modality method claim is currently authorized.

