# ChampHLA completion-goal progress

## 2026-09-15T10:16:07Z — goal initialized from attached contract

Checkpoint: the 808-line completion contract was read completely and activated
as the durable project objective. The canonical checkout is
`/mnt/d/Users/uonur/Downloads/champhla_publication_candidate_stage` on branch
`roihu-first-closure-20260909` at commit
`9b95de2ab8f63362e1062bcc62502a83765dd835` (Git tree
`22ea2f9861eea1b1dae40942177b8bf57cfdf77c`). The only pre-existing local
worktree change is the user's untracked
`docs/CHAMPHLA_COMPLETION_GOAL_PROMPT.md`; it is textually identical to the
reattached source when end-of-line differences are ignored.

Verified evidence:

- Roihu checkout: clean and at the same commit.
- Fundamental supported-environment test: Roihu Python 3.11.15,
  `tests/test_consensus.py`, 15 passed in 0.13 seconds.
- Frozen workflow and comparator manifests were read; the comparator manifest
  is `FROZEN` with no blocking deployable methods.
- WES `pilot-wes-stagev3` is validated.
- RNA `pilot-rna-stagev3` scheduler jobs and finalizer completed successfully;
  its final ledger still requires a post-terminal validation check.
- WGS `pilot-wgs-stagev3` remains active under the frozen commit and must not be
  restarted or promoted early.

Material findings entered in `artifacts/goal/issue_ledger.tsv` include the
unfinished same-resource firewall/evaluation, unresolved independent WGS/WES/RNA
lanes, missing figure and submission systems, human signature/review/declaration
boundaries, incomplete historical-ablation provenance, and the contradictory
repository-wide IMGT/HLA 3.59.0 field in `confirmation_protocol.json` despite
caller-specific releases spanning 3.11.0–3.63.0.

Next action: poll the existing WGS stagev3 jobs to terminal state. Then validate
all six authoritative pilot samples and calculate the three-modality storage
gate exclusively from the stagev3 ledgers. Do not sign the amendment or submit
capacity/production work before that gate passes.

## 2026-09-15T10:20:40Z — authoritative RNA pilot validated

Roihu scheduler jobs `1332826` and `1332827` completed successfully. The
`pilot-rna-stagev3` staging ledger contains two `validated` rows with verified
paired-source checksums and extracted hashes. Its caller ledger contains eight
`validated` rows, and both samples have `CALLERS_COMPLETE` markers. This closes
the RNA portion of the technical-pilot matrix without using any accuracy or
truth information.

WES remains validated. WGS staging job `1332816_1` remains active; its caller,
finalizer, and second-sample chain are correctly dependency-held. Next action:
continue polling that existing WGS chain to terminal state, then validate the
combined expected pilot matrix before running the storage assessor.

Read-only baseline inspection also found that the local verification script
currently writes the combined main-plus-supplement audit to the main-only audit
path and never regenerates the distinct combined audit. Tracked platform test
attestations cite commit `1b3209f`, and the release freeze predates the current
workflow. These are now explicit issues `MS-010`, `REPRO-002`, and
`RELEASE-002`; fixes are deferred until the active WGS pilot has finished.

The local data-location audit found no BAM, CRAM, FASTQ, SIF, or AGC payloads.
It did identify two fail-closed release defects: the remote compact exporter is
extension-only rather than path-allow-listed, and the release freezer can
include ignored/untracked prohibited files outside its small directory deny
set. Issues `RELEASE-003` and `RELEASE-004` record the required shared,
path-aware export policy and regression coverage.

## 2026-09-15T10:33:48Z — supported Roihu baseline verified

Both complete test entrypoints passed on the clean Roihu checkout at
`9b95de2` under Python 3.11.15: each reported 229 passed and 3 skipped. The
first combined SSH invocation lost its closing transcript after its pytest
process had already exited; process inspection showed no live pytest, and a
separate bounded invocation of `python -m pytest` passed in 0.76 seconds. The
local WSL Python remains unsupported at 3.8.10, while the existing ignored
Windows virtual environment uses 3.13, so a supported local 3.10/3.11
environment remains an explicit reproducibility issue rather than being
papered over by syntax changes.

Both authoritative WES POLYSOLVER wrapper audits were inspected directly on
Roihu. Each verified source hash
`bcc2c4c5a1c7b2f3547111002aa2f4e7e244115aac3dd2c263f0b47d9e7b7fab`,
spec hash `eb904bedeacb7c223f836c112c7c9aa173113b04c799a4cf226b76ba5e691afa`,
`chr6`, three interval substitutions, one temporary-directory substitution,
and transformation version `polysolver-hg38-contig-compatibility-1`. The
remaining pending lifecycle labels are issue `PROV-005`; updating them must not
erase the separately disclosed region-plus-unmapped method change or its
historical-comparability caveat.

Rendering inventory found no available Pandoc, LibreOffice, SVG renderer,
Matplotlib, Pillow, or python-docx in the checked local/Roihu environments.
There is also no bibliography, submission directory, or figure subsystem.
Issues `MS-011` and `FIG-003` separate the scientific citation work from a
pinned minimal rendering-environment plan; no package is installed and no
download is authorized by this checkpoint.

The baseline claim/gate dependency map is now recorded in
`artifacts/goal/baseline_artifact_map.tsv`. It distinguishes current valid
frozen inputs from provisional development evidence, remotely dynamic pilot
records, stale generated audits, and missing blocked outputs. Narrative status
is not used as proof that any scientific gate passed.

Storage-code inspection found that preflight overwrites fixed audit/inventory
paths and that sequential transient headroom is computed by subtracting two
marginal p95 values instead of taking p95 of each sample's measured
`peak-retained` difference. Issues `PROV-006`, `WF-003`, and `WF-004` require
versioned post-pilot allocation evidence and the correct per-sample transient
quantile before the authoritative storage gate can run.

## 2026-09-15T16:21:31Z — local fail-closed design audit while WGS remained unobserved

The renewed-certificate Roihu status request was not executed because the
remote-access approval service rejected it after reaching its usage limit. No
claim is therefore made that the WGS chain completed, and no pilot job,
workflow file, tracked configuration, or manuscript was changed. The same
Slurm job chain remains the authoritative object to re-query; it must not be
restarted merely because observation was unavailable.

The local audit resolved the implementation shape of the storage correction.
`assess_storage` already collapses repeated caller rows into one internally
consistent measurement per sample and enforces exactly two validated samples
for each modality. Its sequential calculation is nevertheless unsafe when the
sample having the largest peak is not the sample having the largest retained
footprint. The correction must form `peak_disk_bytes - retained_disk_bytes`
for each accepted sample first and then apply nearest-rank p95 to those paired
differences. A regression fixture must use crossed peak/retained maxima so the
old marginal subtraction demonstrably underestimates the required headroom.

The allocation evidence also needs a time and identity boundary. The current
CLI accepts any passed `environment_inventory.json`, even if its `free_bytes`
predates the pilots. The final design will require a newly generated,
non-overwriting inventory identified by an evidence ID, carrying its UTC
creation time and repository/workflow identity. The storage gate will record
and validate that exact inventory hash. Capacity and production submission
will consume an explicitly promoted current pointer rather than whichever
fixed filename was most recently overwritten.

Two release paths were confirmed to need a shared allow-list policy rather
than more filename extensions. The compact Roihu exporter currently accepts
every small Markdown/TSV/JSON/image/log/text file below a run root, which can
include truth-bearing or secret-bearing files. The Python release freezer
recursively includes untracked and ignored files and then asserts
`truth_data_included: false` without proving it. The required correction is a
single path-aware policy used by both tools: explicit roots and artifact
classes, deny rules for truth/secrets/raw sequencing/containers/references/work
trees, rejection of symlinks and unexpected untracked files, and independent
archive extraction verification. This work remains deferred until the active
pilot chain is authoritatively terminal, honoring the instruction not to
change the repository during the six-sample run.

A temporary, non-writing execution of the current claim auditor established a
second concrete manuscript defect. The main-only audit remains correctly
blocked only by pending claim `BM-007`. The combined audit has an additional
failure because supplementary S9 cites the historical POLYSOLVER `357/390`
result without a registry reference, and no corresponding caller-level row
exists in `result_registry.tsv`. The number is scientifically necessary to
explain why the corrected region-plus-unmapped configuration is not directly
comparable, so it must not simply be deleted. Issue `MS-007` now requires an
immutable development-role registry row for that historical configuration and
an explicit source artifact. The main-text sentence saying the final compact
table comes from “prospective evaluation output” must also distinguish the
corrected same-resource table from later donor-independent prospective
inference. No manuscript content was changed during the active pilot.

The staged benchmark-readiness implementation was also compared with the full
completion contract. `configs/release_requirements.json` and
`audit_release_readiness` correctly separate an honest same-resource benchmark
release from the donor-independent three-modality claim, but they are not a
full goal-completion audit. They currently have no required gates for the
figure manifest and named visual review, bibliography lineage, Markdown/DOCX
structural parity, extracted-archive verification and prohibited-payload scan,
or completion of donor-independent WGS, WES, and RNA. Issue `RELEASE-005`
therefore requires a separate machine-readable completion schema and auditor;
the benchmark gate must retain its narrower documented meaning rather than be
silently repurposed.

## 2026-09-15T16:30:16Z — paused at approval-gated Roihu observation

The read-only Roihu query remains unavailable for a third consecutive goal
turn after the approval service's explicit rejection. Because elapsed time is
not terminal-state evidence, neither the WGS pilot nor the combined six-sample
matrix has been promoted. All safe pre-terminal local audit work is now
consolidated. `reports/HUMAN_ACTION_REQUIRED.md` requests one narrowly scoped,
read-only SSH observation and gives the exact validation and resume sequence.
The goal is paused, not complete; no job or tracked scientific artifact was
changed.

## 2026-09-15T22:38:30Z — local evidence and submission controls completed

The renewed CSC identity/certificate pair has matching public-key fingerprints
and is valid through 2026-09-17T00:38:03Z. A read-only Roihu check confirmed the
remote checkout remains clean at the unchanged base commit. No remote file,
Git reference, or Slurm job was changed.

The local candidate now includes a shared path-aware release/export policy,
Git-tracked-only release freezing, deterministic archive creation and safe
extraction verification, and a twelve-gate full-completion auditor that remains
strictly separate from the staged same-resource benchmark gate. Historical
comparator rows now have exact artifact/source/workflow hashes; caller-specific
database releases are authoritative; and the workflow lock includes all
imported confirmation modules. Focused release, manifest, workflow, figure, and
submission tests pass.

Two non-performance figures are generated reproducibly in SVG, PDF, and
300-dpi PNG: the workflow/truth firewall and the cohort/evidence map. The three
result-dependent main figures remain explicitly blocked and have no placeholder
outputs. The figure validator reports zero malformed generated rows while
correctly leaving the overall gate false.

A journal-neutral submission package now contains byte-deterministic DOCX and
Markdown versions of the main text and supplement, embedded figures, an
11-entry primary-source bibliography, and checksum manifest. Semantic
Markdown/DOCX parity and bibliography audits pass. The declaration audit stays
false because six author-only fields are intentionally marked for input; no
authorship, conflicts, funding, ethics, or acknowledgements were inferred.
Both main-only and combined claim audits have zero failures and only the
evidence-dependent BM-007 warning.

Next action: complete the final local diff/release-policy review, create one
reversible reviewed commit, and prepare a narrowly scoped authorization packet
for GitHub push/CI, Roihu synchronization and preflight, followed by exactly six
stagev4 technical-pilot samples. The amendment remains unsigned until the
post-pilot storage gate passes.

## 2026-09-15T23:18:00Z — public independent-cohort search repeated

The required bounded pre-submission search was repeated against primary papers
and official repository records without opening any controlled dataset. No
qualifying public donor-independent WES or RNA-seq cohort was found, and no
criterion was relaxed. The superficially promising 829-subject WES benchmark is
entirely 1000 Genomes and therefore same-resource. GeT-RM supplies public direct
HLA truth for 108 cell lines, but its 108-sample exome benchmark was explicitly
private; the currently linked public sequencing is WGS for 70 samples, not the
required WES lane. The Finnish 93-94-subject dataset uses dedicated MHC capture,
not WES.

For RNA-seq, public AFGR reads do not close the gate: independent MKK retains
126 samples after QC and the high-resolution truth is controlled, while the
remaining AFGR subjects overlap 1000 Genomes and cannot be pooled. FNLCR has 96
controlled subjects and DICE has 91 controlled donors. The exact dispositions
and sources are versioned in
`decisions/20260909_public_validation_gap_assessment.json`; DATA-006 and
DATA-007 remain genuine external evidence gaps rather than benchmark-release
failures.

## 2026-09-15T23:42:00Z — HPRC truth gate corrected before acquisition

A static audit found that the HPRC protocol validator treated
`source_commit` as a SHA-256 field, which would reject a normal 40-character
Git commit, while the assembly-truth builder did not itself require the
protocol to be frozen. Both defects are corrected before any HPRC download or
execution. Frozen protocols now accept exact 40- or 64-character Git commits,
require an immutable AGC file rather than a directory URL, and continue to
require SHA-256 identities for archives, tool artifacts, wrappers, and
references. The truth builder requires the frozen protocol as an input, records
its hash, rejects ambiguous Boolean audit fields, and cannot run on the current
draft protocol.

The regenerated 49-file workflow lock is
`e8d7b33f5ccb4f5b69fd03083ba5f7046c33342c03b317f74911019cce4a36f3`,
and all comparator workflow hashes agree with it. The selected local suite now
passes 127 tests. HPRC remains correctly blocked on separately authorized tool,
reference, and AGC acquisition; no remote data or software was accessed.

## 2026-09-15T17:21:49Z — access restored; WGS stagev3 failed closed

The user explicitly authorized read-only Roihu access. The existing CSC
identity/certificate pair validated successfully, so renewal is not needed.
Roihu remained clean at commit `9b95de2ab8f63362e1062bcc62502a83765dd835`
and tree `22ea2f9861eea1b1dae40942177b8bf57cfdf77c`.

The WGS chain is terminal. HG00096 staging, all five callers, its finalizer, and
`CALLERS_COMPLETE` marker validated. HG00097 staging job `1332819` failed with
exit 1 after 7,141 seconds; caller job `1332820` was dependency-cancelled and
finalizer `1332821` failed because staging/caller/marker evidence was
incomplete. The quarantine preserves a 116,748,342-byte partial BAM with hash
`277b1449f177200ed527f21200cf2dc4d6ae5fc6fbfc6b38819309e57d08d59e`.
Its only samtools diagnostic was `error reading file` followed by the quoted
HTTPS CRAM URI. The classifier labelled this `deterministic` as an unclassified
exit 1, so the permitted second in-job attempt never ran.

This is a narrow transport-classifier coverage defect: the allowlist handled
CRC, closing, timeout, reset, DNS, HTTP 429/5xx, seek, and EOF forms but not
samtools' generic quoted remote-URI read-error form. A local correction now
classifies only `error reading file "http(s)://..."` as transient; a local-path
read error remains deterministic. Targeted classifier tests pass.

The local correction set also computes the sequential-storage p95 from paired
per-sample transient measurements, introduces environment-inventory schema v2
with safe evidence IDs, UTC capture time, freshness and post-pilot constraints,
and creates non-overwriting `preflight/<evidence_id>` evidence with a separate
explicit `CURRENT` promotion script. Targeted tests currently pass (23 tests in
the selected control/script set). The workflow lock and full test suite still
need refreshing and verification.

Because this changes the candidate commit, the validated WES/RNA stagev3 rows
and successful HG00096 WGS row remain calibration evidence. The final storage
gate requires new stagev4 runs of all two samples in all three modalities on
one reviewed commit. No retry or new remote job has been submitted. The earlier
access action packet was removed after verified resolution; the next packet
will separately scope the push, Roihu synchronization, preflight, and six pilot
jobs.
