# Champion–Challenger versus simple consensus: publication note

## Result

A truth-free family of four simple rules was evaluated: always-call plurality, strict
majority, two-thirds consensus, and unanimity. Each caller contributes one unweighted vote
for its exact canonical unordered two-field pair. No confidence, reliability, calibration,
gene-specific tuning, or modality-specific tuning is used.

The two-thirds rule provides the most useful predeclared operating point for external
validation:

| Modality | Two-thirds consensus | Callability | Champion–Challenger | Delta | Subject-cluster p | Holm p |
|---|---:|---:|---:|---:|---:|---:|
| WES | 328/390 (84.10%) | 339/390 (86.92%) | 365/390 (93.59%) | +9.49 points | 3.73e-9 | 2.61e-8 |
| RNA-seq | 296/321 (92.21%) | 308/321 (95.95%) | 304/321 (94.70%) | +2.49 points | 0.00781 | 0.0391 |

The tests are exact subject-cluster sign-flip tests. Holm correction covers all eight
policy-by-modality comparisons.

## What the result means

Champion–Challenger makes exactly the same calls as two-thirds consensus wherever that
baseline emits a call: zero call differences across 339 WES and 308 RNA-seq loci. It then
calls every unresolved locus and is correct at 37/51 WES and 8/13 RNA-seq loci.

This supports the claim that Champion–Challenger **preserves high-agreement consensus and
resolves otherwise uncalled discordant loci**. It does not show better conditional accuracy
on high-agreement loci: two-thirds consensus has 96.76% WES and 96.10% RNA-seq called-only
accuracy because it abstains on difficult loci.

## Manuscript use

Use fixed-denominator accuracy as the primary endpoint and place callability and called-only
accuracy immediately beside it. Show the complete agreement-threshold curve rather than only
the two-thirds row. Keep MajorityVote, MV-floor, best-single-tool, and MetaConsensus in the
main comparison table.

Suggested result sentence:

> Relative to a fixed two-thirds caller-agreement consensus, Champion–Challenger preserved
> every accepted consensus call while resolving additional discordant loci, increasing
> fixed-denominator exact accuracy by 9.49 percentage points in WES and 2.49 points in
> RNA-seq (Holm-adjusted subject-cluster permutation p=2.61e-8 and 0.039, respectively).

This is an internal discovery analysis. Because the threshold family was introduced during
method development, the two-thirds rule, endpoint, software hashes, and multiplicity policy
are now frozen. The sentence should become a primary manuscript claim only if the direction
replicates on an independent cohort without changing the rule.
