# ChampHLA plurality-centered publication candidate

This repository is the canonical source for the ChampHLA research workflow and
cross-assay benchmark. Its recommended default is pair-level plurality
consensus, serialized as `SimplePluralityLex`. The repository does not claim a
novel voting algorithm. It supplies reproducible caller execution,
harmonisation, provenance, a truth firewall, multi-assay evaluation, and
explicit evidence boundaries.

The plurality-centered benchmark is active. The Guarded Champion–Challenger
manuscript and all historical integration analyses are retained as superseded
or supplemental evidence. Historical WGS performance is invalid because the
input slice omitted relevant alternate-contig mappings and mates. A final
three-modality performance claim is blocked until the partial-call correction
is rerun end to end on WGS, WES, and bulk RNA-seq.

## Primary contract

- One equal vote is cast for each complete canonical unordered allele pair from
  an intended-use caller.
- Partial and missing calls do not vote. A missing second allele is never
  inferred to be homozygous.
- Maximum pair support wins; lexicographic order resolves an exact top tie; no
  complete caller pair yields `no_evidence`.
- Exact unordered two-field HLA-A/B/C accuracy at the fixed eligible-truth
  denominator is primary. Callability, called-only and allele accuracy, ties,
  partial calls, missingness, stratified results, and resources are secondary.
- Predictions are truth blind, checksum-frozen, and joined to truth once.

The detailed contract is in
[`shared_methods.md`](manuscripts/shared/shared_methods.md), the comparator
freeze is in [`comparator_manifest.json`](configs/comparator_manifest.json),
and remote execution is described in
[`PLURALITY_RERUN_RUNBOOK.md`](docs/PLURALITY_RERUN_RUNBOOK.md).

## Current evidence state

Development WES and RNA-seq plurality results are registered but remain
provisional because the primary framing followed inspection of development
outcomes. Invalid WGS rows remain visible for diagnosis and cannot enter gates,
pooled estimates, the abstract, or primary figures. Dataset discovery has 11
candidate assay rows and 10 truth-source rows, each tracked through per-lane
states rather than a single readiness label.

The author must sign
[`20260908_consensus_primary_amendment.json`](decisions/20260908_consensus_primary_amendment.json)
before an external prediction freeze. The comparator manifest must also contain
verified executable hashes and status `FROZEN`.

## Verification

Install with a supported Python version and run the complete local check:

```bash
python -m pip install -e '.[test]'
bash scripts/run_local_verification.sh
```

The script runs both pytest entry forms, reconciles registered development
counts, regenerates tables and discovery states, audits the main text and
supplement together, and verifies that unfinished prospective gates fail
closed. Ubuntu and Windows CI run Python 3.10 and 3.11.

The release readiness audit can be run directly:

```bash
python -m champhla_recovery.cli audit-release-readiness \
  --project-root . \
  --config configs/release_requirements.json \
  --output artifacts/release_readiness.json
```

An exit status of 2 is expected until the signed amendment, corrected remote
runs, WGS audit and named review, prediction freeze, one-time truth join,
prospective evaluation, cross-platform attestations, claim audits, and clean
Git freeze all pass.

## Roihu-first execution

Roihu under CSC project `project_2008084` is the authoritative compute and
large-data environment. The laptop checkout must never contain sequencing
reads, alignments, containers, reference bundles, AGC assemblies, or Nextflow
work directories; `.gitignore` enforces those classes. Only compact text,
figures, logs, checksums, tests, and release archives return to the laptop.

The live preflight is currently blocked because the CSC SSH certificate expired
on 2026-09-02. After the author renews it, install this exact reviewed commit at
`/scratch/project_2008084/champhla_plurality`, create the external roots defined
in `configs/roihu_paths.env`, and run:

```bash
bash scripts/roihu_preflight.sh
bash scripts/setup_roihu_test_env.sh
```

Production arrays cannot start until the amendment is manually signed, every
comparator is checksum-pinned and the comparator manifest is `FROZEN`, the
environment inventory passes, and pilot disk measurements clear the storage
gate. See [`PLURALITY_RERUN_RUNBOOK.md`](docs/PLURALITY_RERUN_RUNBOOK.md) for
the exact staged commands and separate benchmark/independent readiness gates.
