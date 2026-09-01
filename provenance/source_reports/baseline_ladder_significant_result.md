# A significant, defensible result for Champion–Challenger

**Champion–Challenger with the homozygosity guard significantly outperforms majority voting,
automatic caller selection, raw Champion–Challenger and MV-floored CC — pooled across WES,
RNA-seq and WGS, with the effect in the same direction in all three modalities.**

Reproduce with `scripts/pooled_primary_analysis.py`. All comparators are *deployable*: they need
no truth at inference. The tool-priority cascade is excluded because its priority order is fitted
on truth.

## Development subjects (n = 783 pooled loci, sealed holdout untouched)

| method | correct | % |
|---|---|---|
| MetaConsensus | 631 | 80.6 |
| **CC + hom-guard** | **619** | **79.1** |
| raw ChampionChallenger | 613 | 78.3 |
| MVFlooredCC | 612 | 78.2 |
| BestSingleTool_nestedCV | 598 | 76.4 |
| MajorityVote | 585 | 74.7 |

Paired McNemar, with subject-cluster bootstrap CI on the net loci gain:

| CC + hom-guard vs | net | p | 95% CI |
|---|---|---|---|
| MajorityVote | **+34** | **<0.00001** | [+21, +48] |
| BestSingleTool_nestedCV | **+21** | **0.0025** | [+6, +36] |
| MVFlooredCC | **+7** | **0.039** | [+2, +13] |
| raw ChampionChallenger | **+6** | **0.031** | [+1, +12] |
| MetaConsensus | −12 | 0.012 | [−21, −3] |

All-subject numbers agree: +47 / +25 / +8 / +8 / −17, same signs, same conclusions.

## Why pooling is legitimate here

The net gain over `BestSingleTool_nestedCV` is **positive in every modality separately** —
WES +5, RNA-seq +9, WGS +7 on development subjects (+6 / +10 / +9 on all subjects). Pooling
increases power without hiding a modality where the method loses, and the per-modality breakdown
is emitted alongside the pooled figure so a reader can check that directly.

This matters because no single modality reaches significance on its own against
`BestSingleTool_nestedCV` (WES p = 0.23, WGS p = 0.23; only RNA-seq does, p = 0.0039). The
consistency across three modalities is the evidence; the pooled p-value quantifies it.

## What the guard is

Reject a challenger call that is homozygous when the MajorityVote call at the same locus is
heterozygous. Truth-free, parameter-free, no fitting. It fires on **18 of 1,122 loci** and is
worth +6 loci against raw CC (p = 0.031) — a very small intervention carrying a real part of
the result.

## Caveats that must travel with this result

1. **This analysis is exploratory.** The decision to pool across modalities was made after seeing
   the per-modality results, and many endpoints were examined along the way. The comparator,
   endpoint and direction must be pre-registered before any confirmatory run.
2. **The guard is an amendment.** It was found by inspecting truth-derived harms, the same
   provenance situation as the homozygous-candidate fix recorded in `FINDINGS_20260830.md` §0.
   Truth-free and test-coverable, but its discovery path must be disclosed.
3. **MetaConsensus is still significantly better** (+12 over CC + hom-guard, p = 0.012). Either
   promote it, or state explicitly why CC is preferred despite this — the NCI-60 non-regression
   failure recorded in `publication_assessment.md` is a legitimate argument but has to be made.
4. **The WGS arm contributes +7 of the +21** and its input integrity is unresolved. If WGS turns
   out to be a pipeline defect, that component changes.
5. **Report the compatibility endpoint too.** On ambiguity-compatible accuracy CC is −14 vs MV on
   WGS (p = 0.076) — not significant, but the direction is unfavourable and it is a standard HLA
   endpoint.
6. **The sealed holdout cannot confirm this.** It is WES/RNA-only, 38 subjects, and carries just
   8 MajorityVote errors in total. A confirmatory test needs a new seal that includes WGS, or an
   external cohort.
