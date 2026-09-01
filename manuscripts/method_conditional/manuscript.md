# ChampHLA Guarded Champion–Challenger for resolving low-agreement HLA consensus calls

## Status

This manuscript is a conditional working copy and is not submission-ready. It
cannot become the active manuscript until the frozen decision in
`decisions/20260901_publication_recovery.json` passes without policy changes.

## Proposed claim

The method preserves every call made by a fixed two-thirds consensus and uses a
frozen Champion–Challenger policy only at lower-agreement loci. The intended
claim is specifically that it resolves consensus failures at a fixed
denominator. It is not a claim of better conditional accuracy on already-called
loci and is not a claim that every always-call comparator is inferior.

## Required evidence

The manuscript will be activated only if subject-unseen 1000G WGS, WES, and
RNA-seq comparisons are each positive, Holm-significant, and have
subject-clustered confidence intervals above zero; per-gene harm and plurality
non-regression gates must pass; and independent HPRC WGS must exclude a
regression worse than two percentage points. No development number will be
copied into the abstract as confirmation evidence.

## Methods

The method, endpoint, comparators, and truth firewall are defined in
`manuscripts/shared/shared_methods.md`. The policy, caller panels, thresholds,
normalisation, references, and statistical family are frozen before truth is
joined.

## Results

No confirmatory result is available. Result prose and tables will be rendered
from `result_registry.tsv` only after the immutable prediction bundle is joined
to independently locked truth.

## Decision

If all gates pass, this copy becomes the method manuscript. If one or two
modalities pass, their results move to the benchmark manuscript with restricted
claims. If no modality passes, this copy is archived and baseline searching or
post-result retuning stops.

