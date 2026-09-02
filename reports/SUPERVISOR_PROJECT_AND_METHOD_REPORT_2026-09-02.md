# ChampHLA project and method report for supervisors

**Date:** 2 September 2026  
**Purpose:** Summarise the consensus baselines and ChampHLA methods that were
tested, explain what worked and what did not, and define the scientifically
defensible publication route.

## One-page summary

ChampHLA began as a multi-caller HLA typing workflow with a
Champion–Challenger (CC) decision layer. During development, we tested several
simple consensus baselines and more complex extensions because raw CC did not
consistently outperform the strongest simple baseline, unweighted
plurality/MajorityVote, in WES and RNA-seq.

The main findings are:

1. **The software and harmonisation framework worked.** The project integrates
   multiple HLA callers, normalises their outputs, preserves audit information,
   and supports WGS, WES, and RNA-seq with modality-specific caller panels.
2. **MajorityVote is already very strong in WES and RNA-seq.** On the current
   1000 Genomes development data it produced 365/390 correct WES loci (93.59%)
   and 306/321 correct RNA-seq loci (95.33%). There is little remaining room for
   a caller-output-only method to improve these modalities.
3. **Raw CC did not demonstrate superiority over MajorityVote in WES or
   RNA-seq.** Depending on the frozen artifact, WES raw CC is reported as
   364/390 or 365/390, and RNA-seq as 304/321. This one-call WES discrepancy
   must be reconciled. In either case, raw CC does not clearly beat the
   always-call plurality baseline.
4. **CC significantly outperformed an abstaining two-thirds consensus in
   internal development analyses.** The simple baseline calls only when at
   least two-thirds of callable tools agree. The advantage comes from resolving
   loci where this baseline reports `no_consensus`; it is not improved accuracy
   on the high-agreement loci already called by the baseline.
5. **TwoThirdsGuardedCC is the best current CC variant in the internal data.**
   It preserves every two-thirds consensus call and uses CC only at
   low-agreement loci. It achieved 369/390 in WES (94.62%) and 306/321 in
   RNA-seq (95.33%). This is +10.51 and +3.12 percentage points versus the
   abstaining baseline, but only +1.03 and 0.00 points versus always-call
   plurality. These are development results, not external confirmation.
6. **The more complex experimental methods did not produce a publishable new
   accuracy method.** WeightedConsensus gave only a small, non-significant WES
   gain; MetaConsensus regressed in NCI-60; RefFormer failed its nested gate;
   and EvidenceGatedCC failed to rescue the MajorityVote errors using read-level
   evidence.
7. **Historical WGS performance is invalid.** The old extraction omitted HLA
   alternate-contig reads and discarded mates outside the chromosome-6 slice.
   The manuscript's current WGS headline must be removed. A repaired
   full-CRAM pilot now passes technical validation for two subjects, but it is
   not a performance study.
8. **The default publication route should be a benchmark/software paper.** A
   method-focused paper should be activated only if frozen unseen-cohort tests
   confirm Guarded CC without hiding MajorityVote and other strong comparators.

## What the project currently contains

The project has two related contributions:

- **Software contribution:** reproducible execution of several HLA callers,
  canonical two-field allele normalisation, modality-specific caller panels,
  consensus generation, failure tracking, and per-call audit traces.
- **Method contribution:** Champion–Challenger and later variants intended to
  choose a reliable final allele pair when callers disagree.

The intended-use panels currently frozen for class-I evaluation are:

- WGS: HLA-HD, Kourami, OptiType, SpecHLA, and T1K.
- WES: HLA-HD, OptiType, POLYSOLVER, SpecHLA, and T1K.
- RNA-seq: ArcasHLA, HLA-HD, OptiType, and T1K.

All primary evaluation is restricted to HLA-A, HLA-B, and HLA-C at exact,
unordered, two-field resolution. Missing calls count as incorrect in the
fixed-denominator endpoint.

## Baselines and comparators that were evaluated

| Method | What it does | What happened | Current role |
|---|---|---|---|
| Plurality/MajorityVote | Each callable tool casts one equal vote for its complete allele pair; the most frequent pair wins with deterministic tie-breaking. | Strongest simple always-call baseline in WES and RNA-seq: 365/390 and 306/321 correct. Raw CC did not improve it consistently. | Mandatory main comparator. |
| Per-allele independent vote | Votes separately for individual alleles and combines the two most supported alleles. | Slightly worse than pair-level plurality: 363/390 WES and 305/321 RNA-seq. | Exploratory baseline; not preferred. |
| Tool-priority cascade | Uses the first callable tool in a development-ranked priority order. | Worse than plurality in WES and RNA-seq: 360/390 and 303/321. It is also partially trained rather than a purely naive baseline. | Secondary comparator only. |
| Fixed single caller | Uses one caller, such as OptiType, HLA-HD, or T1K. | Demonstrated strong modality dependence and did not provide a universally best policy. | Required descriptive comparator. |
| Nested best-single-tool | Selects the strongest single tool within training folds. | Useful for testing whether an ensemble exceeds the best caller ceiling; did not provide a universal improvement. | Required secondary comparator. |
| SimpleStrictMajority | Calls only when one pair receives more than half of callable votes. | CC improved fixed-denominator performance by filling no-consensus loci, while preserving the baseline calls. | Exploratory abstaining baseline. |
| SimpleTwoThirdsConsensus | Calls only when one pair has at least two-thirds support; otherwise reports `no_consensus`. | Produced 328/390 WES and 296/321 RNA-seq correct calls. It has very high called-only accuracy but lower callability. | Frozen baseline for the narrow “resolution of consensus failures” claim. It cannot be the only comparator. |

The baseline ladder was useful because it showed that changing the form of a
simple vote does not create a better general method. Pair-level plurality
remained the strongest simple always-call baseline. The two-thirds rule is
useful for measuring whether CC can resolve disagreements, but it was selected
during development after several agreement thresholds were examined.
Therefore, its present significance is discovery evidence and requires frozen
confirmation.

## ChampHLA methods and experimental extensions

| Method | Intended improvement | Result | Decision |
|---|---|---|---|
| Raw Champion–Challenger | Use a benchmark-selected champion by default and allow a weighted challenger to replace it after support and margin gates. | Did not beat plurality in WES/RNA-seq. The old apparent WGS advantage was based on invalid inputs and cannot be used. | Retain as the core explanatory ablation, not a confirmed universal winner. |
| WeightedConsensus | Weight callers by training-derived reliability rather than equal votes. | In the subject-separated development subset it rescued two WES loci without harm, a +0.725-point gain, but McNemar P=0.50. RNA-seq gain was zero. | Exploratory; no headline claim. |
| MV-floor | Use raw CC only in training-supported vote strata and otherwise fall back to MajorityVote. | It equalled MajorityVote in WES/RNA-seq. Its reported WGS gain is invalid because of the historical input extraction. | Remove as the manuscript's operating headline. |
| TwoThirdsGuardedCC | Preserve calls accepted by the two-thirds baseline; use raw CC below two-thirds agreement and deterministic plurality only when raw CC is unavailable. | Best internal CC variant: 369/390 WES and 306/321 RNA-seq. Its benefit is resolution of abstentions. | Freeze for confirmation; do not claim external superiority yet. |
| MetaConsensus | L2-regularised candidate-level logistic scorer using vote, caller, modality, gene, and confidence features. | Passed an internal gate mainly because of apparent WGS improvement, but NCI-60 WES and RNA-seq each regressed by 2.44 points versus MajorityVote. | Negative/exploratory ablation, not the central method. |
| RefFormer | Transformer candidate reranker using pinned IMGT sequence representations and caller-context tokens. | Failed the prespecified nested gate. It was below the strongest baseline in WGS, tied in WES, and lower in RNA-seq; the IMGT-aware model did not outperform its no-IMGT neural ablation. | Negative ablation; no further GPU scaling. |
| EvidenceGatedCC | Use independent read-level IMGT 31-mer evidence to override MajorityVote only when evidence is strong. | Rescued 0/23 development WES MajorityVote errors and 0/9 RNA-seq errors. Candidate expansion increased the theoretical oracle, but the evidence scorer could not identify the correct candidate reliably. | Prespecified gate failed; method development stopped. |
| Algorithm-family and correlation-aware gates | Reduce the effect of callers that may share algorithms or reference blind spots. | Did not materially change held-out performance; reported differences were at most approximately 0.5 points. | Sensitivity analysis only. |

## What genuinely worked

### 1. Reproducible multi-caller integration

The strongest current contribution is the workflow rather than a universally
superior accuracy algorithm. The project provides harmonised caller outputs,
canonical pair handling, explicit missingness, intended-use panels, checksums,
and auditable consensus decisions. These are valuable because HLA callers have
incompatible formats, references, dependencies, and modality restrictions.

### 2. A clear understanding of consensus saturation

The work demonstrates that caller agreement is already high in WES and
RNA-seq. MajorityVote reaches 93.59% and 95.33%, respectively, and the
caller-candidate oracle leaves little recoverable error. This explains why
additional caller-output-only models did not improve reliably.

### 3. Guarded resolution of low-agreement loci

TwoThirdsGuardedCC has a precise and defensible interpretation: retain the
simple consensus when agreement is high, and invoke CC when the simple rule
would abstain. In the internal data it added correct calls without changing
accepted two-thirds consensus calls. This supports a claim about **increased
fixed-denominator call resolution**, not improved conditional accuracy on
already-called loci.

### 4. Transparent negative results

MetaConsensus, RefFormer, and EvidenceGatedCC were tested with explicit kill
criteria rather than being repeatedly tuned until a favourable result appeared.
Reporting these results as ablations strengthens the benchmark paper by showing
where complex ensemble methods fail and why strong simple voting is difficult
to beat.

### 5. WGS pipeline repair

The repaired full-CRAM process preserves mates and HLA alternate-contig reads.
For HG00096 and HG00097, all five callers completed and all 30 expected
caller-by-locus records passed checksum, native-output, normalisation, and
parser round-trip checks. This confirms the repaired technical workflow, but
two subjects cannot establish accuracy or significance.

## What did not work or cannot currently be claimed

- Raw CC does not currently beat MajorityVote in WES or RNA-seq.
- The historical WGS numbers are invalid and must not appear in the abstract,
  results headline, or performance figures.
- MV-floor has no WES/RNA advantage over MajorityVote and has no valid WGS
  evidence.
- MetaConsensus did not transfer to NCI-60 without regression.
- RefFormer did not pass either the overall improvement gate or the IMGT-value
  ablation gate.
- EvidenceGatedCC did not recover the observed MajorityVote errors.
- Significance versus the two-thirds baseline cannot be presented as
  superiority over MajorityVote because the two-thirds baseline abstains.
- Existing 1000G and NCI-60 analyses are development or previously examined
  evidence; they are not untouched independent confirmation.
- External confirmation is incomplete. The frozen WGS and WES rosters are
  feasible, but only 99 strictly eligible unseen RNA donors were identified,
  below the locked minimum of 130.

## Numerical and naming problems requiring correction

Before manuscript editing, the following must be reconciled in one locked
table:

1. The raw-CC WES count is 364/390 in the baseline-ladder artifact but 365/390
   in the two-thirds significance artifact.
2. MetaConsensus WGS is registered as 220/411 in one frozen artifact and
   230/411 in the baseline ladder. Neither value can support a claim until WGS
   is rerun from valid full-CRAM inputs.
3. The pooled raw-CC comparison reports 669/711, whereas the later Guarded CC
   totals 675/711. These are different methods and must not both be called
   simply “Champion–Challenger.”
4. The existing manuscript describes MV-floor as the operating method, while
   the recovery plan freezes TwoThirdsGuardedCC as the conditional method.
5. “MajorityVote,” “plurality,” “SimpleConsensusBaseline,” “two-thirds
   consensus,” raw CC, and Guarded CC must be defined separately and used
   consistently in text, tables, and figures.

Until these points are resolved, no numerical method claim should be placed in
the abstract.

## Recommended manuscript strategy

### Default: benchmark/software manuscript

The most defensible paper should focus on:

- multi-caller workflow integration and reproducibility;
- canonical HLA call harmonisation and auditability;
- modality-specific caller behaviour;
- MajorityVote saturation and candidate ceilings;
- why learned and read-evidence extensions did not generalise;
- WGS extraction and failure-propagation lessons;
- transparent reporting of both positive and negative ablations.

Guarded CC can be described as an internally promising mechanism for resolving
low-agreement loci, but not as a universally superior method.

### Conditional: method-focused manuscript

This route should be used only if the unchanged Guarded CC rule succeeds on
frozen unseen cohorts. The main results must still display MajorityVote,
raw CC, MV-floor, best-single-tool, MetaConsensus, and individual callers.
Success only against an abstaining baseline is insufficient for a broad claim
that ChampHLA improves HLA typing accuracy.

## Remaining work in priority order

1. Reconcile the conflicting raw-CC and MetaConsensus counts and regenerate one
   authoritative result table.
2. Obtain named human sign-off for the 30 technically reviewed WGS pilot
   records and add 20 valid records to meet the 50-record production audit.
3. Replace the permissive legacy Nextflow failure policy with a fail-closed,
   tested confirmation runner.
4. Run truth-blind predictions for the frozen 120-subject WGS and 89-subject
   WES rosters, checksum them, and only then join the locked 2014 laboratory
   truth once.
5. Run HPRC WGS as independent validation after the internal WGS pipeline gate
   passes.
6. Treat the current 99-subject unseen RNA cohort as infeasible for the locked
   three-modality claim unless a legitimate orthogonal truth cohort is found;
   do not lower the target after viewing results.
7. Update the benchmark manuscript from the authoritative registry. Remove the
   invalid WGS headline, distinguish fixed-denominator from called-only
   accuracy, and move FIMM, LOH, HED/survival, class II, and long-read analyses
   outside the main story.

## Recommended conclusion to communicate

ChampHLA has produced a useful and potentially publishable multi-caller
benchmark and software framework. The evidence does not currently support a
general claim that raw Champion–Challenger beats MajorityVote. The most
promising CC result is Guarded CC's ability to resolve loci where a conservative
two-thirds consensus abstains while preserving high-agreement calls. That
result is internally significant but still needs frozen confirmation. The
failed MetaConsensus, RefFormer, and EvidenceGatedCC experiments should be
reported transparently as informative ablations rather than omitted.

## Evidence files

This report was prepared from the isolated publication-candidate repository:

- `result_registry.tsv`
- `artifacts/source_results/primary_modality_results.tsv`
- `artifacts/source_results/publication_assessment.json`
- `artifacts/source_results/final_internal_decision.json`
- `artifacts/verification/discovery_reproduction.json`
- `cohorts/frozen_initial_rosters.json`
- `artifacts/wgs_pilot_review/REVIEW_STATUS.md`
- `manuscripts/claim_audit.tsv`
- `decisions/20260901_publication_recovery.json`

All values in this report are development or secondary evidence unless
explicitly described otherwise. No result is currently authorised as a
confirmed abstract-level performance claim.
