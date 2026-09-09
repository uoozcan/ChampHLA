# ChampHLA robustness implementation status

Status date: 2026-09-09.

## Implemented in the canonical repository

- The active route is the plurality-centered benchmark; the Guarded CC method
  draft is a superseded provenance archive.
- `SimplePluralityLex` is the primary machine identifier and `MajorityVote` is
  a documented legacy alias.
- Native `-` second alleles remain partial, retain allele 1, never imply
  homozygosity, do not vote, and are incorrect for genotype accuracy.
- Plurality output includes complete/partial/missing caller accounting,
  support, ties, supporting callers, version, and source hashes.
- Raw CC and all other integration outputs are optional; plurality can be
  generated from caller calls alone.
- Evaluation is reference-method configurable and reports the three explicit
  plurality gates, modality-local Holm families, simultaneous donor-clustered
  bounds, allele and called-only accuracy, caller call states, and unpooled
  independence strata. The v1 Guarded CC evaluator is read-only legacy code.
- The result registry resolves each entry to one canonical long-format row and
  validates values and both canonical/source checksums. Development,
  same-resource, independent, exploratory, and invalid roles are explicit.
- The imported discovery program preserves source checksums for 11 dataset
  candidates, 10 truth sources, 391 crosswalk records, and 111 truth-free pilot
  records. Per-dataset lane states are generated fail closed.
- Main and supplementary manuscripts now center the workflow, assay-dependent
  callers, selection-free plurality, negative integration ablations, and the
  WGS boundary. Invalid WGS diagnostics are confined to the supplement.
- Freeze paths are POSIX-normalized across operating systems. Ubuntu and
  Windows CI cover Python 3.10 and 3.11 and both pytest entry forms.
- Roihu interfaces now validate and freeze truth-free manifests, enforce exact
  production matrices, inventory repositories/tools/references/containers,
  maintain caller-level state ledgers, gate scaling on pilot disk use, submit
  dependent Slurm waves, reparse native outputs, quarantine failures, plan
  cleanup dry-run first, and export only compact allow-listed artifacts.
- The 11-subject NCI-60 RNA lane has an official ENA checksum-manifest builder.
  The HPRC lane has truth-blind 120-subject selection, a checksummed short-read
  manifest gate, and conservative dual-method assembly-truth reconciliation.
- Benchmark-release readiness is separate from independent three-modality
  readiness; missing independent WES/RNA evidence cannot invalidate an honest
  corrected same-resource benchmark release.

## Deliberately blocked

- The primary amendment remains `UNSIGNED_AWAITING_AUTHOR`.
- Several comparator versions and executable hashes remain unresolved, so the
  comparator manifest is not frozen.
- The WGS, WES, and RNA-seq end-to-end corrected reruns have not been executed.
- The production WGS audit lacks a complete batch, environment manifest, and
  named ten-record-per-caller review.
- No corrected prediction freeze or one-time truth join exists.
- Current WES and RNA-seq estimates remain retrospective development evidence;
  historical WGS results remain invalid.
- No independent locus-level evaluator has reproduced corrected headline
  counts, and remote/CI attestations for this revision do not yet exist.
- The benchmark claim audit and release audit remain fail closed until the
  prospective outcome selects the predeclared wording.
- The `roihu-cpu` SSH certificate expired on 2026-09-02, so no current Roihu
  inventory, transfer, or job submission can be performed yet.
- This checkout has no GitHub remote. Cross-platform CI cannot be triggered
  until the author supplies the existing repository URL.
- HLA-ASM and Immuannot remain intentionally `PREPARED_NOT_FROZEN`; their Roihu
  versions and artifact hashes must be verified before HPRC truth jobs.
- Independent public WES (minimum 89 donors) and RNA (minimum 130 donors)
  cohorts with direct complete two-field HLA-A/B/C truth remain evidence gaps.

## Next authorized execution sequence

1. Renew the CSC certificate, supply the existing GitHub remote, and confirm
   local/GitHub/Roihu commit and tracked-tree hashes.
2. Obtain the author signature and freeze all comparator artifact/version
   metadata from the live Roihu inventory.
3. Execute the corrected truth-blind runs using
   `docs/PLURALITY_RERUN_RUNBOOK.md`.
4. Pass the complete WGS technical and named human-review audit.
5. Freeze predictions, validate the freeze, and join truth once.
6. Run the same-resource evaluator and independent count reproduction; then
   execute NCI-60 and HPRC as separately frozen evidence strata.
7. Regenerate the registry, main table, supplement, claim audits, test
   attestations, clean Git release manifest, and DOI archive bundle.
