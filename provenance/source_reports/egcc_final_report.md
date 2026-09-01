# Evidence-Gated Champion-Challenger pilot: final internal decision

Decision date: 2026-08-31

## Outcome

The pilot does **not** support EvidenceGatedCC as a new accuracy-method headline. The
prespecified development gate was not met, so the sealed 38-subject holdout remains unopened,
external validation was not started, and no production ChampHLA interface or manuscript was
changed.

The defensible publication route is the benchmark/software paper. EvidenceGatedCC belongs in
that paper as a transparent negative ablation, while the modality/gene-weighted consensus can
be reported as a small exploratory WES result.

## Decisive corrected results

| Modality | Development loci | MajorityVote | Best weighted consensus | Rescue/harm | Delta | Read-evidence rescues among MV errors |
|---|---:|---:|---:|---:|---:|---:|
| WES | 276 | 253 | 255 | 2 / 0 | +0.725 points | 0 / 23 |
| RNA-seq | 231 | 222 | 222 | 0 / 0 | 0.000 points | 0 / 9 |

The WES weighted result has McNemar p=0.50. It is harm-free but below the two-point gate and
not statistically significant. In RNA-seq, the caller-pair oracle equals MajorityVote at
222/231, so caller reweighting has no development-set headroom.

Candidate generation was not the limiting component: corrected expanded candidate oracles
were 275/276 WES and 230/231 RNA-seq. The blocker was ranking. With corrected assembly-specific
read extraction and the best truth-blind 31-mer configuration, cross-modality top-call
concordance was 0.2982 versus 0.971 for MajorityVote; the ranker rescued none of the 32
development loci on which MajorityVote was wrong.

## Important audit findings

- The original RNA extraction queried GRCh38 coordinates against GRCh37 Geuvadis BAMs. Those
  RNA conclusions were superseded after assembly-aware extraction and header verification.
- Candidate generation omitted homozygous pairs for caller alleles outside the read top-N.
  The correction is test-covered and explicitly recorded as a post-registration amendment.
- The first read scorer had singleton-group certainty bias, inert IDF weighting, and
  hash-order nondeterminism. These defects were corrected and preserved in the audit trail.
- The current WGS table has systematic homozygous over-calling (58.4% called versus 13.4%
  truth) and must not be interpreted as WGS accuracy.

## Verification and preservation

- 29/29 isolated tests pass.
- Corrected Roihu shards: 130 WES and 106 RNA-seq.
- Large read/index artifacts remain on Roihu; compact results are synchronized to D.
- The sealed holdout remains sealed because no development method met its gate.
- Existing ChampHLA, RefFormer, MetaConsensus, manuscripts, defaults, and source projects were
  not changed.

## Manuscript recommendation

Frame the paper around reproducible eight-caller integration, normalization, modality-specific
failure modes, leakage-controlled benchmarking, auditability, and the empirical saturation of
MajorityVote. Report the two WES rescues as exploratory and EvidenceGatedCC as a negative
read-evidence ablation. Do not claim superiority over MajorityVote or external generalization.

Before submission, repair or exclude the WGS arm, run the already-planned independent public
cohorts for the benchmark claims, complete cohort-flow/denominator/ethics reporting, pin all
software and reference versions, and archive the reproduction bundle at a DOI.

## Added simple consensus baseline

Two truth-free baselines were added without changing the existing methods:

- `SimplePluralityLex`: one unweighted vote per callable caller's exact unordered pair,
  with a lexicographic tie. It is identical to the locked MajorityVote on all 711 WES/RNA
  loci, demonstrating that the observed MajorityVote result is reproducible with a very
  simple rule.
- `SimpleStrictMajority`: emit only when one exact pair has more than half of callable-tool
  votes; otherwise emit `no_consensus`.

Against strict majority, Champion–Challenger raises fixed-denominator correctness from
359/390 to 365/390 in WES and from 296/321 to 304/321 in RNA-seq. The interpretation is
coverage, not better accuracy on agreed loci: at the strict baseline's coverage, both methods
have exactly 359/376 WES and 296/308 RNA-seq correct calls. Champion–Challenger additionally
gets 6/14 WES and 8/13 RNA-seq no-consensus loci correct.

This is a useful secondary story—Champion–Challenger resolves otherwise uncalled loci—but it
does not replace the strongest-baseline comparison. MajorityVote, MV-floor, best-single-tool,
and MetaConsensus must remain visible in the manuscript.

The subsequent complete agreement-threshold analysis identified a frozen two-thirds
consensus operating point. CC preserved every two-thirds consensus call and increased
fixed-denominator correctness by 9.49 points in WES and 2.49 points in RNA-seq. Exact
subject-cluster tests remained significant after Holm correction across all eight
policy-by-modality comparisons (adjusted p=2.61e-8 and 0.0391). This is an internal discovery
result and requires unchanged external validation before becoming the primary method claim.

An experimental `TwoThirdsGuardedCC` policy now enforces this property by construction:
protect every two-thirds consensus call and invoke raw CC only on lower-agreement loci. It
protects 647 loci and invokes CC on 64; its current outputs are identical to raw CC on all
711 loci. In the complete individual-caller comparison, CC is significantly better after
subject clustering and Holm correction than WES HLA-HD, SpecHLA, and T1K and RNA-seq T1K.

For a single primary baseline analysis, the frozen two-thirds rule is exposed as
`SimpleConsensusBaseline`. Pooled WES/RNA correctness is 624/711 for the baseline and 669/711
for CC (difference +6.33 points, subject-clustered 95% CI +4.47 to +8.37; exact clustered
sign-flip p=2.91e-11). CC preserves all 647 baseline calls and correctly resolves 45 of its
64 no-consensus loci. This remains internal discovery evidence until unchanged external
confirmation.
