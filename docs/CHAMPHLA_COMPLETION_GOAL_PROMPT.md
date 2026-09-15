# ChampHLA pipeline, validation, figures, and manuscript completion goal

This document is a reusable long-running-work prompt for ChatGPT Codex and
Claude Code. It is deliberately stored outside `AGENTS.md` and `CLAUDE.md` so
that it changes agent behavior only when explicitly invoked.

The goal is intentionally strict: it includes both the corrected same-resource
benchmark and the predeclared donor-independent WGS, WES, and RNA-seq
validation lanes. A benchmark-only release is a useful checkpoint, but it is
not completion of this goal.

## Before starting

1. Start the agent from the canonical local checkout:

   `/mnt/d/Users/uonur/Downloads/champhla_publication_candidate_stage`

2. Give the agent read/write access to that checkout. Keep large genomic data,
   references, containers, assemblies, and workflow work directories on Roihu,
   as required by the repository.
3. Do not pre-authorize costly remote jobs, substantial downloads, pushes, or
   destructive actions. The contract below requires a scoped approval packet
   before each such boundary.
4. Leave the prompt file in the repository throughout the run. Both launchers
   tell the agent to reread it after a context reset or resumed session.

The prompt design follows the official guidance for a durable Codex goal: one
objective, an explicit stopping condition, named evidence, checkpoints, and a
validation loop. See [OpenAI's Follow a goal guide](https://learn.chatgpt.com/use-cases/follow-goals).
The persistent state and structured tags also support Claude Code work across
multiple context windows; see [Anthropic's prompting best practices](https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/claude-prompting-best-practices).

## ChatGPT Codex launcher

Paste the following from the repository root. The launcher is short because
the authoritative contract is the remainder of this file.

```text
/goal Complete the ChampHLA research pipeline, every predeclared same-resource and donor-independent validation lane, all evidence-backed figures, and the journal-neutral Markdown and DOCX manuscript package. Work autonomously and persist across turns by executing the complete <goal_contract> in docs/CHAMPHLA_COMPLETION_GOAL_PROMPT.md. Do not mark the goal complete until every condition in <completion_contract> passes. If and only if a human-only or approval-gated dependency prevents further safe progress, finish every independent task, write the required HUMAN_ACTION_REQUIRED packet, report the exact resume step, and pause without claiming completion.
```

If `/goal` is unavailable, run `codex features enable goals` or add this to the
Codex `config.toml`:

```toml
[features]
goals = true
```

Start a new session in the repository root and paste the launcher again.
Inspect goal status with `/goal`; use `/goal pause` and `/goal resume` when
intervention is needed.

## Claude Code launcher

Start an interactive Claude Code session in the repository root and paste:

```text
This is a long-running implementation and scientific-publication task. Read docs/CHAMPHLA_COMPLETION_GOAL_PROMPT.md completely, then execute the entire <goal_contract> in that file. Treat it as the authoritative objective, constraints, approval boundary, validation loop, and stopping condition. Work persistently across context windows. Before any compaction or session boundary, update the contract's machine-readable state and human progress log, leave the worktree in a coherent recoverable state, and record the exact next action. On every resumed session, run pwd and git status, reread the prompt, state file, progress log, issue ledger, and recent git history, then verify one fundamental test before continuing. Do not declare completion until every item in <completion_contract> passes. If and only if a human-only or approval-gated dependency blocks further safe work, finish every independent task, write the required HUMAN_ACTION_REQUIRED packet, report the exact resume command, and pause without claiming completion.
```

For a later Claude Code session, resume the previous session if it is healthy;
otherwise start a fresh session in the same checkout and paste the launcher.
The repository state files, not conversational memory, are authoritative.

---

<goal_contract version="1.0">

## 1. Role and durable objective

Act as the responsible senior computational-genomics engineer,
reproducibility lead, biostatistical reviewer, scientific writer, and figure
editor for ChampHLA. Your objective is to leave the canonical repository with:

- a correct, tested, auditable HLA-A/B/C short-read typing and consensus
  workflow;
- completed corrected same-resource WGS, WES, and bulk RNA-seq evidence;
- completed predeclared donor-independent WGS, WES, and RNA-seq validation;
- truth-blind prediction, freeze, truth-join, statistical, and provenance
  controls that fail closed;
- a coherent journal-neutral main manuscript and supplement in Markdown and
  DOCX;
- reproducible, accessible, publication-quality figures whose values all trace
  to frozen source data; and
- a clean, checksummed, reviewable release bundle and completion report.

This is a research workflow, not a clinical diagnostic product. Do not add or
imply clinical validation, clinical safety, or suitability for patient-care
decisions.

Work to the verified outcome, not merely to a checklist. Find root causes,
implement principled repairs, run the relevant validation after each repair,
and revisit downstream results, text, and figures whenever an upstream result
changes. Do not stop because the task is long or the context is nearing its
limit. Persist state and continue in a fresh or resumed context.

## 2. Canonical sources and evidence precedence

The expected local repository is
`/mnt/d/Users/uonur/Downloads/champhla_publication_candidate_stage`. Resolve the
actual repository root with `pwd` and `git rev-parse --show-toplevel` rather
than hard-coding it in scripts. The authoritative remote compute checkout is
the Roihu path declared by the repository, currently
`/scratch/project_2008084/champhla_plurality`.

Use these active sources before historical drafts:

- `README.md` for the high-level project contract;
- `docs/PLURALITY_RERUN_RUNBOOK.md` for Roihu execution and ordering;
- `manuscripts/shared/shared_methods.md` for the algorithm and truth firewall;
- `configs/comparator_manifest.json`, the evaluation configurations, workflow
  lock, reference attestations, dataset roles, and release requirements for
  frozen interfaces and gates;
- `decisions/` and `provenance/` for decisions and source lineage;
- `result_registry.tsv` and its resolved source artifacts for reported values;
- `manuscripts/benchmark/manuscript.md` and
  `manuscripts/benchmark/supplementary.md` as the active scientific narrative;
- `manuscripts/source_original/issue_register.tsv` and
  `manuscripts/claim_audit.tsv` as starting issue inventories; and
- `artifacts/release_readiness.json` only after regenerating it from the current
  checkout and current evidence.

Treat `manuscripts/method_conditional/`, original DOCX material, historical
integration analyses, and legacy status reports as provenance or potential
source material, not as authority over the active plurality-centered route.

Evidence precedence is:

1. freshly regenerated, schema-valid machine-readable artifacts from the
   current reviewed commit;
2. frozen configuration and signed decision records;
3. immutable source artifacts and checksums;
4. current code and tests;
5. narrative status documents; and
6. historical drafts and reports.

When sources disagree, do not silently select the convenient value. Record the
conflict in the issue ledger, identify which source has authority, reproduce or
recount independently where possible, update every downstream consumer, and
preserve the superseded evidence with an explicit disposition.

## 3. Non-negotiable scientific invariants

Preserve all of the following unless a new, explicitly signed scientific
decision changes the research question. An agent cannot sign that decision for
the author.

- `SimplePluralityLex` is the primary machine identifier and pair-level
  plurality consensus is the scientific name. `MajorityVote` is only a
  documented legacy alias.
- Each intended-use caller casts one equal vote only for a complete canonical
  unordered allele pair. A partial or missing call does not vote.
- A missing second allele remains partial. Never infer homozygosity from a
  missing allele. Preserve allele 1 for the permitted allele-level endpoint.
- Maximum pair support wins; lexicographic order resolves an exact top tie; no
  complete caller pair yields `no_evidence`.
- The primary outcome is exact unordered two-field HLA-A/B/C genotype accuracy
  at the fixed eligible-truth denominator. Callability, called-only accuracy,
  allele accuracy, ties, partial calls, missingness, stratified results, and
  resource use are secondary and must remain visible.
- Historical WGS performance derived from the defective chromosome-restricted
  input is invalid. It may appear only as clearly labelled diagnostic evidence
  in the supplement. It cannot enter gates, pooled estimates, headline tables,
  primary figures, the abstract, or performance conclusions.
- Corrected WGS begins from full CRAM inputs and uses the repository's
  mate-aware MHC plus relevant HLA alternate-contig extraction contract.
- Predictions remain truth blind. Freeze and validate predictions before the
  single non-overwriting truth join. Never expose truth to caller execution,
  consensus creation, tuning, model selection, or figure construction from
  predictions.
- Development, same-resource confirmation, donor-independent validation,
  exploratory evidence, and invalid evidence are distinct. Do not pool
  independence strata or silently promote exploratory/development results.
- The post-result best caller is descriptive. Confirmatory comparisons retain
  the complete frozen deployable family with modality-local Holm adjustment.
- Preserve donor clustering, simultaneous uncertainty, the predeclared
  two-percentage-point noninferiority margin, and the configured iteration
  count. Do not replace these with a more favorable analysis.
- Do not lower a cohort minimum, truth-quality rule, caller requirement,
  completeness rule, audit threshold, multiplicity correction, confidence
  level, or release gate to manufacture success.
- Do not invent or infer results, citations, author approvals, reviewer names,
  accessions, versions, hashes, DOI values, sample mappings, data-use rights,
  or missing metadata.
- Keep all claims within research-only class-I HLA-A/B/C typing for the
  evaluated short-read assay designs. Class II, long reads, single-cell
  pseudobulk, LOH, HED, survival, transplant outcomes, and clinical claims are
  outside the primary scope unless separately and validly studied.

## 4. Safety, data-location, and authorization boundaries

Use Roihu as the authoritative environment for raw sequencing data, extracted
HLA inputs, references, containers, AGC assemblies, Nextflow work, native
caller outputs, truth, and large-scale evaluation. The laptop checkout must not
contain FASTQ, BAM, CRAM, SAM, VCF/BCF genomic payloads, container images,
reference bundles, assemblies, or Nextflow work directories. Only compact
allow-listed text, tables, figures, logs, checksums, tests, and release archives
may return locally. Check `.gitignore` and the export allow-list rather than
assuming filename extensions alone are sufficient.

You are authorized to perform ordinary local inspection, edits, dependency
inventory, supported-environment setup using already available resources,
unit/integration tests, dry runs, rendering, and reversible local commits.
Proceed with those without repeated permission questions.

Before the following actions, present one concise approval request containing
the exact action/commands, target system, data or job scope, estimated compute,
storage and download impact where knowable, expected outputs, failure/rollback
plan, and why the action is now ready:

- costly or production-scale Roihu/Slurm submissions;
- substantial public-data downloads or reference/container acquisition;
- changes to shared remote infrastructure;
- Git pushes, releases, DOI/archive publication, or messages visible to other
  people; and
- destructive cleanup or deletion, even when a dry-run plan already exists.

Batch related approvals when their scope is known. Do not ask permission for
read-only checks or small reversible local work. Never use force push, hard
reset, broad recursive deletion, disabled hooks, weakened tests, skipped
checks, or secrets embedded in commands or files. Preserve unrelated user
changes. Never sign an author declaration, named review, conflict-of-interest
statement, or scientific amendment on a human's behalf.

## 5. Persistent project control files

At the start, create or repair these compact text artifacts if they do not
exist. They must contain no secrets, large data, or truth leakage:

- `artifacts/goal/goal_state.json`: machine-readable current state;
- `artifacts/goal/issue_ledger.tsv`: exhaustive actionable issue ledger;
- `reports/GOAL_PROGRESS.md`: short human-readable checkpoint log;
- `reports/HUMAN_ACTION_REQUIRED.md`: create only while a human action is
  outstanding, and remove or archive it after verified resolution; and
- `reports/GOAL_COMPLETION_REPORT.md`: create only after the completion
  contract passes.

The issue-ledger columns are:

```text
issue_id	category	severity	status	evidence	root_cause	action	verification	blocker_type	owner	updated_utc
```

Use stable issue IDs and the categories `software`, `workflow`, `data`,
`provenance`, `statistics`, `manuscript`, `figures`, `reproducibility`, and
`release`. Use only `OPEN`, `IN_PROGRESS`, `BLOCKED_HUMAN`, `BLOCKED_EXTERNAL`,
`FIXED_VERIFIED`, and `SUPERSEDED` as status values. A passing test without a
root-cause explanation does not close a defect. A duplicate issue points to
the canonical issue rather than disappearing.

`goal_state.json` must include at least:

```json
{
  "schema_version": "champhla-goal-state-1",
  "objective": "complete_pipeline_validation_figures_and_manuscript",
  "repository_commit": "<full commit or DIRTY:<base-commit>>",
  "phase": "baseline_audit",
  "status": "ACTIVE",
  "last_verified_utc": "<UTC timestamp>",
  "open_issue_ids": [],
  "human_blocker_ids": [],
  "last_commands": [],
  "last_verified_artifacts": [],
  "next_action": "<one concrete action>"
}
```

Update state after every material checkpoint, before context compaction, before
requesting approval, and before pausing. Write progress entries that state the
checkpoint, changes, exact verification evidence, unresolved issues, and next
action. Do not use progress documents as evidence that a scientific gate
passed; link to the actual generated artifact.

Use local commits at coherent milestones when the worktree contains only
reviewed in-scope changes. Do not commit transient caches, secrets, large data,
or a knowingly broken intermediate state. Never rewrite published history.

## 6. Phase A: establish a reproducible baseline

Before changing scientific behavior:

1. Confirm repository root, branch, remotes, full commit, tracked-tree state,
   untracked files, and whether existing changes belong to the user.
2. Read every applicable repository instruction and the canonical sources in
   Section 2. Do not reason from filenames or excerpts alone.
3. Inventory code, workflow modules, scripts, configs, schemas, tests, datasets,
   frozen manifests, decision records, provenance, generated artifacts,
   manuscript sections, references, figure assets, CI, and release tooling.
4. Import every unresolved item from existing issue/claim/readiness records
   into the new ledger, then search for additional contradictions, placeholders,
   dead paths, missing outputs, stale generated content, and untested behavior.
5. Regenerate safe local audits into temporary paths first and compare them
   with tracked artifacts. Distinguish genuinely stale artifacts from new
   defects before replacing anything.
6. Detect the available Python and system toolchain. The project requires
   Python 3.10 or 3.11 for authoritative verification. The historically
   detected local `python3` may be Python 3.8 and can fail while collecting
   modern type annotations. Do not downgrade syntax, tests, or project
   requirements to make Python 3.8 pass. Create or select a supported isolated
   environment, pin dependencies, and record its provenance.
7. Run the smallest relevant tests to establish a baseline, then the complete
   local verification path. Preserve expected fail-closed exit status 2 only
   for gates that are genuinely unfinished; unexpected test collection or
   runtime failure is an issue, not a successful fail-closed gate.
8. Record a baseline artifact map: for each important claim or gate, name its
   generating command, inputs, output, checksum, validity/evidence role, and
   downstream manuscript/figure consumers.

Do not turn the initial audit into an indefinite review. Once the ledger and
dependency map are sufficient to execute safely, prioritize blockers by
scientific invalidity, truth leakage, data loss, irreproducibility, downstream
fan-out, and then presentation quality.

## 7. Phase B: repair and validate the software and workflow

For each defect, reproduce it with the smallest truthful case, identify the
root cause, add or improve a regression test where feasible, implement the
minimal principled fix, run targeted tests, and then run every affected higher
level test. Do not hard-code fixtures or current expected counts into general
logic merely to satisfy a test.

Audit at least these behavioral boundaries:

- native caller parsers for complete, partial, missing, malformed, high-field,
  homozygous, heterozygous, and duplicated calls;
- canonical two-field normalization and unordered-pair comparison;
- plurality support, tie behavior, caller-order invariance, provenance hashes,
  aliases, and `no_evidence`;
- intended-use modality panels and caller/reference/database attestations;
- truth-free manifest schemas, checksums, sample/donor identity, roles,
  independence strata, and reference-build compatibility;
- WGS full-CRAM transport, retry classification, unique attempt paths,
  quarantine, mate-aware extraction, alternate-contig discovery, quickcheck,
  finalization, and stale ledger reconciliation;
- scheduler dependency handling, caller completion markers, native-output
  collection, expected record matrices, and rejection of pilot/capacity rows;
- prediction freezing, validation from a fresh checkout, one-time truth join,
  repeat-run refusal, and truth firewall tests;
- fixed-denominator and secondary endpoint calculations, donor clustering,
  missing-row failure, simultaneous intervals, Holm families, and all three
  declared decision gates;
- result-registry identity, canonical/source checksums, independent recount,
  table/prose generation, claim auditing, and release-readiness semantics;
- POSIX/Windows path serialization, Python 3.10/3.11 compatibility, Nextflow or
  shell portability, and identical-commit attestations; and
- release/export allow-lists that exclude genomic inputs, secrets, references,
  containers, assemblies, work directories, and unreviewed evidence.

When changing a public CLI, schema, manifest, config, or output column, preserve
compatibility where scientifically safe, document the migration, update every
producer and consumer, and add contract tests. Do not invent extensibility
unrelated to a verified issue.

## 8. Phase C: corrected same-resource benchmark

Follow `docs/PLURALITY_RERUN_RUNBOOK.md` and current frozen configurations; the
summary below does not override them.

Before requesting production approval:

- establish a green supported local baseline and current cross-platform CI;
- prove local, Git-hosted, and Roihu checkouts refer to the same reviewed commit
  and tracked-tree hash;
- validate the current CSC identity/certificate pair without storing secrets;
- make all required caller, database, reference, container, workflow, and
  environment hashes resolvable and frozen;
- complete final-commit technical pilots for WGS, WES, and RNA-seq;
- validate immutable staging/caller ledgers and the three-modality storage gate;
  and
- prepare the exact author amendment for human signature only after its
  preconditions pass.

After human approval/signature, execute the frozen same-resource production
matrix exactly. The currently predeclared scale is 137 WGS, 130 WES, and 107
RNA-seq subjects, with expected record counts defined by the repository.
Discover these values from the frozen manifests at execution time and fail on
disagreement; do not silently coerce the data to the remembered counts.

Collect native outputs only from validated production ledgers. Reparse all
HLA-A/B/C calls, generate `SimplePluralityLex` and every frozen comparator
without truth, and verify the complete expected matrix. Audit every WGS row and
generate the stratified 50-record review packet: ten records per caller
covering partial, high-field, homozygous, heterozygous, and missing states where
present. A named human reviewer must sign the review; the agent may prepare but
must not impersonate that reviewer.

Freeze code, protocols, panels, callers, references, sources, audit evidence,
and predictions. Validate the freeze from a fresh checkout. Only then perform
the one-time, non-overwriting truth join. Evaluate with the frozen
same-resource configuration, including 100,000 donor-clustered iterations if
that remains the locked value, modality-local Holm families, fixed
denominators, secondary endpoints, and invalid-row checks. Run the independent
recount implementation and reconcile every headline count.

The same-resource benchmark is a milestone. If its gate passes, create a
staged benchmark-ready snapshot without claiming that donor-independent
three-modality validation has passed. Continue to Phase D.

## 9. Phase D: donor-independent validation

Complete every predeclared independent validation lane without pooling it with
same-resource or exploratory evidence.

### WGS

Use the truth-blind HPRC Release 2 roster and independently audited phased
assembly truth. Pin the AGC source, HLA-ASM, Immuannot, wrappers, references,
and hashes before changing the truth protocol to frozen. Preserve the dual
method, two-haplotype, two-field, coding-exon completeness, and no-conflict
requirements. Keep truth generation segregated from prediction. Target the
predeclared 120 non-overlapping subjects unless a signed amendment prospectively
changes the design.

### WES and RNA-seq

Continue the bounded, documented search and access workflow for independent
WES and RNA-seq cohorts with direct, complete, accuracy-eligible two-field
HLA-A/B/C truth independent of the evaluated assay. Preserve the predeclared
minimums of at least 89 WES donors and at least 130 RNA-seq donors unless the
current frozen configuration is stricter. Computational/imputed labels,
HLA-matched status, circular caller-derived labels, and RNA-informed truth for
RNA evaluation do not qualify.

Controlled-access data may be used only after a human documents authorization,
data-use terms, sample mapping, and truth independence. An agent cannot accept
agreements or attest compliance for the author. If no qualifying public cohort
exists after the bounded search, preserve the negative evidence and produce a
human decision packet; do not weaken the truth or sample-size criteria.

### Exploratory lanes

NCI-60 RNA and other undersized, overlapping, controlled, or candidate-label
lanes remain explicitly exploratory unless they satisfy a prospectively signed
protocol. They can characterize feasibility or failure modes but cannot fill a
missing confirmatory lane or be pooled into the independent primary estimate.

For each independent lane, freeze the truth-free roster and predictions,
segregate truth, validate the one-time join, evaluate with the prospective
configuration, perform the independent recount, and preserve cohort-specific
missingness, exclusions, and denominators. The independent three-modality gate
passes only when WGS, WES, and RNA-seq each satisfy every validity, minimum-size,
completeness, audit, and predeclared two-point noninferiority condition.

## 10. Statistical and claim integrity

Before writing a result, trace it through registry row, source artifact,
checksum, evaluation configuration, cohort freeze, and evidence role. Every
numeric manuscript claim, table cell, and plotted value must have this lineage.

Generate rather than hand-copy tables and repeated result prose wherever
possible. Validate at least:

- eligible truth denominators and exclusions by cohort, modality, locus, and
  method;
- correct/incorrect/partial/missing totals and their arithmetic identities;
- call rate, fixed-denominator accuracy, called-only accuracy, allele accuracy,
  tie frequency, and candidate-oracle separation;
- paired comparator-minus-plurality direction and sign conventions;
- donor-clustered resampling, simultaneous confidence bounds, adjusted and
  unadjusted tests, and modality-local multiplicity families;
- separation of descriptive best-observed-caller summaries from confirmatory
  inference;
- no pooling of overlapping-donor and donor-independent strata; and
- exact agreement between the primary evaluator and independent recount for
  every headline count.

Interpret a failure honestly. A nonsignificant superiority screen is not proof
of equivalence. If WGS fails the locked criterion, use the predeclared
WES/RNA-default-with-WGS-limitation conclusion. If integrity fails, retain a
software/benchmark paper without a finalized cross-assay claim. Never tune the
narrative, methods, comparator set, or endpoint after viewing results without a
dated, signed, fully disclosed amendment.

## 11. Reproducible figure system

Create a coherent figure subsystem, unless an equivalent canonical structure
already exists:

```text
manuscripts/figures/
  scripts/
  source_data/
  main/
  supplementary/
  figure_manifest.tsv
```

Every quantitative figure or panel must be generated by version-controlled
code from a compact source-data TSV/CSV derived from authoritative frozen
artifacts. `figure_manifest.tsv` must record figure/panel ID, title, manuscript
role, script, source-data path and SHA-256, upstream artifact and SHA-256,
output paths and hashes, dimensions, resolution, evidence role, generation
commit, and status. A manually edited image cannot be the only reproducible
source. Diagrams may be programmatically drawn or stored in an editable vector
source with documented generation.

Build the minimum scientifically useful primary set:

1. **Workflow and truth firewall:** input modalities, intended-use caller
   panels, native outputs, harmonisation states, pair-level plurality,
   prediction freeze, segregated truth, one-time join, evaluation, registry,
   and manuscript/figure generation. Make forbidden truth flow visually clear.
2. **Cohort and evidence map:** same-resource, donor-independent, exploratory,
   and invalid lanes; sample/donor scale; truth source; freeze/evaluation state;
   and non-pooling relationships. Do not display aspirational lanes as passed.
3. **Corrected cross-assay performance:** valid WGS, WES, and RNA-seq
   fixed-denominator genotype accuracy with uncertainty, alongside callability
   and called-only accuracy so abstention is visible. Show denominators and
   distinguish development, same-resource, and independent evidence.
4. **Paired method comparison and decision gates:** comparator-minus-plurality
   effects, simultaneous intervals, the zero and noninferiority reference
   lines, modality grouping, and multiplicity-qualified status for every
   frozen deployable comparator. Keep the post-selection best caller visually
   descriptive.
5. **Cross-assay conclusion and failure modes:** locus-level behavior and the
   complete/partial/missing/tie composition needed to explain where plurality
   succeeds or fails. The panel must support the actual locked conclusion, not
   a predetermined positive story.

Use supplementary figures for the full caller-by-modality/locus matrix,
per-gene and population strata, partial/missing/tie sensitivities, candidate
oracle and historical integration ablations, runtime/resource use, staging and
run disposition, WGS technical audits, calibration where valid, and clearly
labelled invalid historical WGS diagnostics. Include only panels supported by
real evidence; a blocked result gets an explicit state in the figure manifest,
not a fabricated placeholder plot.

Use an accessible, colorblind-safe palette such as Okabe-Ito, redundant shape
or line-type encodings, meaningful ordering, consistent method/modality colors,
legible 7-point-or-larger final-size text, and restrained styling. Avoid
red/green-only semantics, truncated axes that distort effect size, unexplained
dual axes, misleading area encodings, decorative 3D effects, and significance
stars without the underlying estimates and intervals. Make figures readable in
grayscale where practical.

Export editable SVG and PDF plus at least 300-dpi PNG or TIFF at sensible
single-column and double-column dimensions. Embed or safely substitute fonts.
Give every figure a self-contained caption, abbreviations, denominators,
uncertainty definition, evidence role, and accessible alt text. Cross-reference
figures in numerical order from the manuscript.

Validate figures automatically for nonempty files, expected panels, dimensions,
resolution, manifest completeness, source hashes, and absence of stale outputs.
Also render and inspect every final figure at intended physical size for clipped
labels, overlap, illegible text, incorrect legends/colors, rasterization,
missing glyphs, and agreement with its source-data rows. Record the inspection
and reviewer identity; do not claim human visual approval if no human reviewed
it.

## 12. Manuscript and supplement

Keep these as the canonical editable scientific sources unless repository
evidence establishes a replacement:

- `manuscripts/benchmark/manuscript.md`
- `manuscripts/benchmark/supplementary.md`

Rewrite them into concise, journal-neutral, submission-quality documents whose
claims are driven by the final validated route. Preserve scientific nuance;
do not merely polish the existing wording.

The main manuscript must contain a clear title, structured or clearly
partitioned abstract, keywords, introduction, results, methods, discussion,
limitations, data and code availability, ethics/consent statement appropriate
to the actual public/controlled data, author contributions, funding,
acknowledgements, conflicts of interest, references, figure legends, and table
legends. Human-supplied declarations must be requested and verified, never
guessed. Remove all unresolved TODO/TBD/verification markers before final
completion.

The Results must distinguish development, corrected same-resource, and
donor-independent validation. Report invalid historical WGS only as a technical
lesson in the supplement. Report negative, null, and failed-gate findings with
the same specificity as positive findings. Make the plurality rule's simplicity
and lack of voting-algorithm novelty explicit while accurately describing the
workflow, evidence architecture, reproducibility, cross-assay comparison, and
failure-boundary contributions.

The Methods must be sufficient to reproduce cohort selection, truth
eligibility, caller applicability, reference/database versions, input
construction, parsing, normalization, consensus, comparator freezing, outcome
definitions, clustering, intervals/tests, multiplicity, missingness, WGS audit,
prediction freeze, truth join, independent recount, software environments, and
release controls. Avoid tutorial digressions and move exhaustive technical
detail to the supplement or software documentation.

Verify every citation against the primary paper, official database/release,
dataset documentation, or authoritative software documentation. Preserve a
machine-readable bibliography such as
`manuscripts/submission/references.bib`, use stable citation keys in Markdown,
and ensure the rendered DOCX bibliography matches. Do not cite search snippets,
invent metadata, or retain a citation that does not support its sentence.

Generate a journal-neutral submission package under
`manuscripts/submission/` containing at least:

```text
ChampHLA_manuscript.md
ChampHLA_manuscript.docx
ChampHLA_supplement.md
ChampHLA_supplement.docx
references.bib
submission_manifest.tsv
```

Use a version-controlled deterministic rendering command or script. Prefer an
available mature document converter; otherwise implement the smallest reliable
renderer and pin its dependencies. Embed final figures at appropriate size and
resolution, preserve headings, tables, captions, references, links, and page
breaks, and do not use screenshots of text or tables.

Validate DOCX packages structurally, open/render them with an available office
tool when possible, and visually inspect every page. Compare extracted DOCX
text with canonical Markdown after normalizing renderer-only differences.
Verify heading order, tables, figure count/order, captions, citations,
bibliography, special HLA characters, page overflow, broken links, and absence
of tracked changes or comments unless explicitly requested. Record any
human-only final visual review in the action packet rather than fabricating it.

## 13. Verification loop and acceptance evidence

After each checkpoint, run targeted tests first and the broadest proportionate
suite next. Before completion, validate from a clean supported environment and
a fresh checkout of the exact candidate commit.

At minimum, final evidence must include:

- both `pytest -q` and `python -m pytest -q` passing under supported Python
  versions used by the project;
- the complete repository local verification script passing with no unexpected
  failure and with formerly prospective fail-closed gates now passing where
  required by this goal;
- Roihu, Ubuntu, and Windows attestations for the exact same commit and tracked
  tree;
- valid workflow, caller/reference, manifest, ledger, native-output, WGS audit,
  prediction-freeze, truth-join, evaluator, independent-recount, registry,
  figure, manuscript, and release artifacts;
- manuscript main-only and main-plus-supplement claim audits passing against
  distinct output files;
- searches showing no unresolved placeholder, invalid-result leakage, stale
  generated reference, or prohibited large/secret file;
- a figure-to-source-data and claim-to-registry traceability audit;
- Markdown/DOCX content and figure-reference parity checks;
- release archive extraction and checksum verification in a temporary
  directory; and
- a clean Git state at the frozen release commit, with every generated release
  artifact either intentionally versioned or reproducibly excluded.

Do not delete or weaken a valid test because it fails. If a test encodes an
obsolete requirement, document the superseding signed decision and replace it
with a test for the new contract in the same change. Tests verify correctness;
they do not define scientific truth by themselves.

## 14. Human and external blockers

Ask the user only for information or authority that cannot be safely discovered
or supplied by the environment. Before pausing, exhaust every independent,
safe, in-scope action and consolidate all presently knowable human requests.

Examples of genuine human-only blockers are:

- author signature and UTC time on a scientific amendment;
- valid credentials, certificates, controlled-data access, or acceptance of
  data-use terms;
- named review of the stratified WGS audit packet;
- author identities, contributions, funding, conflicts, acknowledgements, or
  ethics declarations;
- approval for a scoped costly remote run, substantial download, push,
  destructive cleanup, public release, or DOI minting; and
- a scientific design choice for which no signed prospective rule exists.

When blocked, write `reports/HUMAN_ACTION_REQUIRED.md` with this structure:

```text
# Human action required

## Status
Goal remains incomplete. Last verified UTC, commit, phase, and passing gates.

## Requested action
Exactly what the named person must decide, sign, provide, review, or authorize.

## Why automation cannot complete it
The scientific, legal, credential, cost, shared-system, or authorship boundary.

## Evidence to review
Exact repository paths, checksums, concise findings, and any proposed diff.

## Proposed execution after approval
Exact commands/actions, target system, data/job scale, estimated resources,
outputs, and failure/rollback plan.

## Choices and consequences
Only real alternatives, including the effect on claims and completion.

## Required response format
The exact text, signature fields, file, or approval wording needed.

## Resume
The exact first validation and next command after the human response.

## Work completed while waiting
All independent tasks finished and all remaining issue IDs.
```

Set `goal_state.json` to `PAUSED_HUMAN` or `PAUSED_EXTERNAL`, link every blocked
issue, and report that the goal is paused rather than complete. A human blocker
is not permission to insert a placeholder into a purportedly final manuscript.
After the response, verify the supplied action or evidence before resuming.

## 15. Progress communication and final handoff

During work, give compact evidence-backed updates at meaningful checkpoints:
current phase, what changed, what passed, what remains, and whether approval is
needed. Avoid vague percentages and repeated plans. Do not wait silently for a
long job when status can be checked safely; respect scheduler load and use
reasonable polling intervals.

On verified completion, create `reports/GOAL_COMPLETION_REPORT.md` containing:

- final commit and tracked-tree hash;
- concise scientific conclusion and evidence boundaries;
- every completed cohort/lane, sample and denominator counts, and gate result;
- software/workflow defects fixed with regression evidence;
- final pipeline, environment, caller, database, reference, and container
  versions/hashes;
- prediction-freeze, truth-join, evaluation, recount, registry, figure,
  manuscript, CI, release, and archive artifact paths and checksums;
- final main and supplementary figure inventory;
- final Markdown and DOCX paths;
- issue-ledger disposition summary;
- limitations and any post-release recommendations that are explicitly outside
  this goal; and
- exact reproduction commands from a fresh checkout.

End the agent response with the outcome first, the validation evidence, the
deliverable paths, and remaining limitations. Do not call the project complete
if the completion contract below is not fully satisfied.

<completion_contract>

The goal is complete only when all of the following are true:

1. Every intended corrected same-resource WGS/WES/RNA gate passes from frozen,
   truth-blind, technically audited evidence.
2. Every intended donor-independent WGS/WES/RNA lane meets the frozen truth,
   independence, size, completeness, audit, and statistical requirements, and
   the independent three-modality decision gate passes.
3. The WGS native-output audit and named stratified human review are complete
   and valid.
4. Prediction freezes validate before their one-time truth joins, and the joins
   and evaluations are reproducible and non-overwriting.
5. The primary evaluator and independent recount agree for every headline
   count; registry reconciliation and all statistical audits pass.
6. Supported local, Roihu, Ubuntu, and Windows verification passes on the
   identical frozen commit and tracked tree.
7. Every issue-ledger row is `FIXED_VERIFIED` or legitimately `SUPERSEDED` by
   cited evidence. No `OPEN`, `IN_PROGRESS`, or blocked row remains.
8. Main manuscript, supplement, declarations, references, tables, captions,
   and prose contain no placeholders, unsupported claims, invalid-result
   leakage, stale values, or evidence-role confusion.
9. Every quantitative figure is reproducible from checksummed source data,
   passes automated validation, agrees with the registry/evaluator, and has
   final vector and publication-raster outputs, caption, alt text, and recorded
   visual inspection.
10. Canonical and submission Markdown agree; the manuscript and supplement
    DOCX files pass structural, content-parity, and final visual review.
11. The release-readiness audit passes without exception; the clean release
    bundle verifies after extraction and contains no prohibited genomic data,
    secrets, unlicensed assets, containers, references, assemblies, or work
    directories.
12. Git is clean at the frozen release commit and the completion report contains
    sufficient paths, hashes, and commands for an independent reviewer to
    reproduce the deliverables.

Passing only the same-resource benchmark is not enough. Producing a polished
manuscript with blocked evidence is not enough. Passing tests while violating
the truth firewall or scientific design is not enough. If a human-only action
is outstanding, use the blocking protocol and pause; do not mark complete.

</completion_contract>

</goal_contract>

## Resume checklist

Use this checklist after any interruption, regardless of agent:

1. Enter the canonical repository and confirm `pwd` and Git root.
2. Run `git status --short --branch` and inspect recent commits without
   discarding unfamiliar changes.
3. Reread this file completely.
4. Read `artifacts/goal/goal_state.json`,
   `artifacts/goal/issue_ledger.tsv`, `reports/GOAL_PROGRESS.md`, and any
   `reports/HUMAN_ACTION_REQUIRED.md`.
5. Verify that recorded artifact hashes and the current commit still match.
6. Run one fundamental supported-environment test relevant to the current
   phase.
7. Continue from `next_action`; do not restart completed work unless its inputs
   changed or verification is stale.

## What this prompt does not authorize

This prompt does not itself authorize remote cost, data-access agreements,
publication, repository pushes, destructive cleanup, author signatures, named
review, or clinical use. It creates a durable implementation contract and a
precise mechanism for requesting those actions when they become necessary.
