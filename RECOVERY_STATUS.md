# ChampHLA publication recovery status

## Completed locally and on Roihu

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
- The canonical D-drive and Roihu mirrors are clean at Git commit `b9770d6`
  before this status update.
- The Roihu preflight passes with Python 3.11.15, samtools 1.21, Nextflow
  25.10.2, Apptainer 1.4.5, the required reference/index, caller containers,
  configuration checksums, and adequate scratch availability.
- A repository-local Roihu Python 3.11 environment reproduced all 43 tests at
  the pre-update commit; the current release adds two regression tests for
  environment initialization and artifact integrity and contains 45 tests.
- Two discarded Python 3.9 environment attempts are retained under
  ignored audit names; they contain no analysis output.
- Both bounded full-CRAM WGS pilots (HG00096 and HG00097) completed extraction,
  all five caller runs, compact collection, and a valid truth-blind freeze.
- Automated native-source, checksum, normalization, and parser checks pass for
  all 30 pilot caller/locus records. Independent source-line review also found
  exact agreement for all 30; named human sign-off and the locked 50-record
  production review gate remain incomplete.
- The only discovered general WES/RNA Nextflow runner is in the dirty,
  read-only `pihla-publish` tree and uses a permissive task-failure policy.
  It has not been copied or run for confirmation; a fail-closed runner must be
  isolated and stub/integration-tested before the frozen WES roster is started.

## Roihu gates

1. Obtain named human sign-off for the 30 reviewed pilot records.
2. Accumulate and manually review 20 additional valid records so the frozen
   ten-record-per-caller (50 total) production gate can be met.
3. Require five nonempty caller-native outputs and 100% parser agreement for
   every new batch.
4. Run the frozen WGS/WES rosters in rolling batches without truth.
5. Freeze predictions, then join sealed truth exactly once.

No new WGS accuracy or three-modality method claim is currently authorized.
The two-subject capacity gate also fails by design (only two discordant
subjects), and the strict unseen RNA roster remains below its locked minimum.
