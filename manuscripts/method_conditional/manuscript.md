# ChampHLA Guarded Champion–Challenger for resolving low-agreement HLA consensus calls

## Status: superseded archive

This manuscript is preserved only for provenance. The signed-state-independent
supersession chain in `decisions/publication_route.json` makes the
plurality-centered benchmark the active manuscript. Guarded CC remains an
ablation and comparator; this draft cannot become the primary manuscript by a
favorable future result.

## Proposed claim

The method preserves every call made by a fixed two-thirds consensus and uses a
frozen Champion–Challenger policy only at lower-agreement loci. The intended
claim is specifically that it resolves consensus failures at a fixed
denominator. It is not a claim of better conditional accuracy on already-called
loci and is not a claim that every always-call comparator is inferior.

## Historical activation rule

The earlier design would have required subject-unseen 1000G WGS, WES, and
RNA-seq comparisons are each positive, Holm-significant, and have
subject-clustered confidence intervals above zero; per-gene harm and plurality
non-regression gates must pass; and independent HPRC WGS must exclude a
regression worse than two percentage points. This rule is retained to explain
the archive and no longer controls manuscript selection.

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

This copy remains archived. Future results update the outcome-dependent wording
of the plurality-centered benchmark and do not reactivate this manuscript.
