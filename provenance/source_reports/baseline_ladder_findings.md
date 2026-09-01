# Baseline ladder results — what Champion–Challenger can and cannot claim

All numbers are **development subjects only** (276 WES / 231 RNA / 275 WGS loci). The sealed
holdout was not used to fit or evaluate anything here. Exact McNemar on paired per-locus outcomes.

Reproduction gate: the frozen `MajorityVote` rule was recovered at **100% agreement on all three
modalities** (the single `no_evidence` locus, `NA11829/wgs/C`, is matched by abstention in both).
Recipe: unweighted plurality, alphabetical tie-break, over the curated eligible panel. No weighting,
no calibration, no tool selection.

---

## 1. Champion–Challenger against each simple-consensus baseline

| modality | baseline | baseline | CC | CC wins | CC losses | net | McNemar p |
|---|---|---|---|---|---|---|---|
| WES | Plurality vote | 253 | 252 | 2 | 3 | −1 | 1.00 |
| WES | Per-allele vote | 253 | 252 | 3 | 4 | −1 | 1.00 |
| WES | Tool-priority cascade | 249 | 252 | 7 | 4 | **+3** | 0.55 |
| RNA | Plurality vote | 222 | 222 | 0 | 0 | 0 | 1.00 |
| RNA | Per-allele vote | 222 | 222 | 0 | 0 | 0 | 1.00 |
| RNA | Tool-priority cascade | 222 | 222 | 0 | 0 | 0 | 1.00 |
| WGS | Plurality vote | 110 | 139 | 38 | 9 | **+29** | **<0.001** |
| WGS | Tool-priority cascade | 138 | 139 | 13 | 12 | +1 | 1.00 |
| WGS | Fixed OptiType | 136 | 139 | 14 | 11 | +3 | 0.69 |

**On RNA every consensus scheme ties at exactly 222.** Plurality, per-allele and cascade all land on
the same number as CC, with zero discordant pairs. There is no baseline choice that separates them.

The cascade baseline gives CC +3 on WES, but at p = 0.55 that is not a claim. **CC's only significant
result anywhere is against plurality voting on WGS** — and against the cascade on the same data it is
+1 (p = 1.00). So the WGS win is not "CC beats simple consensus"; it is "CC repairs majority voting's
WGS failure, returning it to roughly best-single-tool level."

---

## 2. What *is* significant: consensus versus single tools

Plurality consensus against each individual caller, development subjects:

| modality | vs OptiType | vs HLA-HD | vs T1K | vs POLYSOLVER / ArcasHLA |
|---|---|---|---|---|
| **RNA** | +9, 0 losses, **p=0.0039** | +6, 0 losses, **p=0.031** | +23, 0 losses, **p<0.001** | ArcasHLA +7, 0 losses, **p=0.016** |
| WES | +4, p=0.29 | +41, **p<0.001** | +30, **p<0.001** | POLYSOLVER +5, p=0.23 |
| WGS | **−26, p=0.0002** | +52, p<0.001 | +28, p=0.0004 | — |

- **RNA-seq is the clean win.** Consensus beats *every* individual caller with **zero losses against
  each one**, all p < 0.05. This is a robust, defensible ensembling claim.
- **WES is a tie against the strong callers.** Consensus is not significantly better than OptiType
  (p=0.29) or POLYSOLVER (p=0.23); it only beats the weak ones.
- **WGS consensus is significantly *worse* than OptiType** (−26, p=0.0002). Plain majority voting on
  WGS is worse than just running one tool. This is a real negative result and it belongs in the paper.

---

## 3. The claim that survives

Three statements are supported by this data without weakening any baseline:

1. **Multi-tool consensus significantly outperforms every individual caller on RNA-seq**
   (all four pairwise tests significant, zero losses in each).
2. **On WES, consensus matches the best callers and protects against picking a weak one** —
   the value is variance reduction, not a higher ceiling.
3. **Champion–Challenger's contribution is robustness where voting collapses.** On WGS, plurality
   falls to 110/275 while CC holds 139/275, close to OptiType's 136. CC delivers consensus behaviour
   on RNA/WES without majority voting's catastrophic WGS failure mode.

Statement 3 is the honest version of the WGS result, and unlike "CC beats the baseline" it does not
depend on which baseline is nominated. It also does not require the WGS input-integrity question to
resolve first, though that question still needs answering.

What is *not* supported: any claim that Champion–Challenger improves accuracy over simple consensus
on WES or RNA-seq. On RNA the two are identical locus-for-locus.
