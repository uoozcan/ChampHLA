# Roihu-first plurality closure runbook

This is the operational boundary for CSC project `project_2008084`. Roihu is
authoritative for raw data, extracted HLA inputs, references, containers,
Nextflow work, native outputs, truth, and evaluation. Historical
`champhla_cc_confirmation` paths are read-only; the new checkout is
`/scratch/project_2008084/champhla_plurality`.

## Manual blockers

No production prediction is authorized until all are true:

1. The author renews the CSC certificate and manually signs
   `decisions/20260908_consensus_primary_amendment.json` with name and UTC time.
2. The existing GitHub remote is supplied. Local, GitHub, and Roihu point to
   the same reviewed commit and tracked-tree SHA-256.
3. `configs/comparator_manifest.json` contains observed versions and immutable
   hashes for every intended-use caller and frozen integration artifact, and
   its status is `FROZEN`. A missing required caller stops production.
4. `bash scripts/roihu_preflight.sh` passes and writes a clean repository,
   environment, storage, reference, software, and container inventory.

Never place truth in a run manifest or prediction-visible directory. Never
download CRAM/BAM/FASTQ, SIF, references, Nextflow work, or assemblies to the
laptop.

## Truth-free manifests and pilot waves

Canonical columns are:

```text
cohort sample_id donor_id modality input_type input_uri index_uri index_checksum source_checksum reference_build read_layout independence_stratum evidence_role
```

Use tab separation. Checksums may be explicit `md5:` or `sha256:` values.
Validation rejects truth columns, duplicate samples, missing pairs/checksums,
unknown roles/modalities, and reference incompatibility.

```bash
freeze_run_manifest --manifest RUN.tsv --expected-counts configs/roihu_storage_targets.json --output RUN.freeze.json
bash scripts/roihu_submit_wave.sh RUN.tsv wgs pilot
bash scripts/roihu_submit_wave.sh RUN.tsv wes pilot
bash scripts/roihu_submit_wave.sh RUN.tsv rnaseq pilot
```

After the two-sample pilots finish, populate their ledger disk fields and run
`assess_roihu_storage`. Scaling requires the modality-specific 95th percentile
and at least `1.25 × projected peak + 100 GiB`. Then run `capacity` (10 samples)
and only afterward `production`. Production enforces 137 WGS, 130 WES, and 107
RNA samples for a same-resource manifest. Expected caller-locus records are
2,055, 1,950, and 1,284; expected plurality rows total 1,122.

The staging job checks source/index hashes. WGS extraction is mate-aware and
uses chr6:28–34 Mb plus all HLA-A/B/C alternate contigs from the pinned
GRCh38DH index. Caller jobs have `afterok` dependencies, reject failed/empty or
placeholder native results, reparse all A/B/C calls, and create completion
markers only after hashes validate. Caller-level ledgers record job, commit,
manifest hash, state, exit, runtime, peak memory/disk, output hash, and retry
lineage.

## Same-resource closure

Collect native outputs and require exactly 5,289 caller-locus records. Generate
`SimplePluralityLex` and all frozen comparators without truth. Audit every WGS
row and create a stratified 50-record packet—ten records per caller spanning
partial, high-field, homozygous, heterozygous, and missing cases where present.
A named reviewer must approve all 50.

Freeze code, protocols, panels, callers, references, sources, audit, and
predictions. Validate from a fresh checkout. Only then run the non-overwriting
`join_truth_once`. Evaluate with 100,000 donor-clustered iterations,
`configs/consensus_evaluation_same_resource.json`, and modality-local Holm
families. Recount with `independent_recount`, which does not import the primary
evaluator. This lane may set `same_resource_benchmark_ready`; it can never set
the donor-independent three-modality gate.

Cleanup remains a dry-run until the prediction freeze validates. Execution may
delete only listed `work/` and `tmp/` paths. Failed results go to quarantine.
Use `scripts/roihu_export_compact.sh` to create a checksummed allow-listed
archive; genomic inputs, containers, references, assemblies, and work data are
excluded.

## Independent validation lanes

For NCI-60, run `scripts/roihu_prepare_nci60_manifest.sh`. It retrieves compact
official ENA file metadata, selects the predeclared 11 strict RNA subjects, and
freezes both FASTQ MD5 values per subject. Submit through the frozen RNA panel,
freeze/join separately, label `exploratory`, and never pool with 1000G.

For HPRC Release 2, build candidates from official metadata, select 120
non-overlapping subjects before truth, verify each CRAM header and checksum,
then use `build_hprc_run_manifest`. Assembly truth is a segregated workflow.
`configs/hprc_truth_protocol.json` must be changed to `FROZEN` only after the
AGC archive, HLA-ASM, Immuannot, and wrappers are pinned. Extract one subject's
two haplotypes at a time. `build_hprc_assembly_truth` accepts a locus only when
both haplotypes are complete, both methods agree at two-field resolution,
coding exons 2/3 are gap-free, and no equal conflicting allele remains.

Independent WES ≥89 and RNA ≥130 with direct complete A/B/C truth remain
blocked if the bounded public search finds no qualifying cohorts. Do not lower
truth or sample-size criteria. `independent_three_modality_claim_ready` remains
false until prospective donor-independent WGS/WES/RNA all pass the 2-point
noninferiority gate.

## Final verification

Both `pytest -q` and `python -m pytest -q` must pass on Roihu Python 3.11.15 and
GitHub Ubuntu/Windows for the identical commit. Regenerate registry rows,
manuscript tables, independent recount, main/supplement claim audits, and the
release-readiness report. A staged benchmark release may proceed when its own
gate passes even while independent three-modality confirmation remains visibly
blocked.
