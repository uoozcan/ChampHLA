# Shared verified methods

ChampHLA harmonises caller outputs into canonical unordered HLA-A, HLA-B, and
HLA-C allele pairs at two-field resolution. Each intended-use caller contributes
at most one callable pair per locus. Missing calls remain explicit.

`SimpleTwoThirdsConsensus` gives every callable caller one equal vote and emits
a pair only when its support is at least two thirds. `TwoThirdsGuardedCC`
preserves every such call, invokes the already frozen Champion–Challenger policy
below the threshold, and uses deterministic plurality only when the raw policy
is unavailable. The method always emits a top call. The truth-derived
homozygosity amendment is excluded.

The primary endpoint is exact unordered two-field correctness at the fixed
eligible-locus denominator. Only loci whose laboratory truth resolves uniquely
to a two-field pair enter that denominator; no representative allele is chosen
from a cross-two-field ambiguity. A missing call or `no_consensus` is incorrect.
Callability and accuracy among called loci are reported next to the primary
endpoint. Always-call plurality, raw Champion–Challenger, MV-floor,
development-frozen best-single-tool, MetaConsensus, and intended-use callers
remain visible comparators.

All predictions are produced without truth-bearing columns and checksummed
before truth is joined once. Statistical inference clusters by subject, uses
paired sign-flip tests and subject bootstrap intervals, and applies Holm
correction to prespecified families. WGS, WES, and RNA-seq have separate input,
reference, caller-panel, and QC records.
