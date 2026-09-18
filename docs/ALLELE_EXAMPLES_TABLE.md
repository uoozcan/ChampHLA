# Allele-Level Call Examples — Understanding HLA Typing Agreement and Errors

## Background

Each row below shows a real sample + HLA locus from the 1000 Genomes benchmark. The "Truth" column is the validated ground truth allele pair. Each tool column shows what that tool called. The final row (★ MajorityVote) shows the consensus call and how many tools agreed on it.

**Accuracy metric:** Exact two-field match — both alleles must match exactly (e.g., `B*08:01 / B*44:03` must not be called as `B*08:01 / B*44:156`).

---

## Case 1 — All Tools Agree and Are Correct (Ideal Case)
**Sample: HG00114 | Gene: HLA-C | WGS**

| Tool | Allele 1 | Allele 2 | Correct? |
|---|---|---|---|
| ArcasHLA | C\*07:01 | C\*07:02 | ✅ |
| HLA-HD | C\*07:01 | C\*07:02 | ✅ |
| Kourami | C\*07:01 | C\*07:02 | ✅ |
| OptiType | C\*07:01 | C\*07:02 | ✅ |
| T1K | C\*07:01 | C\*07:02 | ✅ |
| **★ MajorityVote** | **C\*07:01** | **C\*07:02** | ✅ **(5/5 agree)** |
| **Ground truth** | C\*07:01 | C\*07:02 | — |

**Interpretation:** When all tools independently converge on the same allele pair, the call is almost certainly correct. Unanimous agreement is the highest-confidence scenario. MajorityVote correctly selects C\*07:01 / C\*07:02 and matches the truth exactly.

---

## Case 2 — Majority Correct, One Tool Wrong (Error Filtered by Consensus)
**Sample: HG00097 | Gene: HLA-B | WGS**

| Tool | Allele 1 | Allele 2 | Correct? |
|---|---|---|---|
| ArcasHLA | B\*07:06 | B\*07:44 | ❌ |
| HLA-HD | B\*07:02 | B\*07:02 | ✅ |
| Kourami | B\*07:02 | B\*07:02 | ✅ |
| OptiType | B\*07:02 | B\*07:02 | ✅ |
| SpecHLA | B\*07:02 | B\*07:02 | ✅ |
| T1K | B\*07:02 | B\*07:02 | ✅ |
| **★ MajorityVote** | **B\*07:02** | **B\*07:02** | ✅ **(5/6 agree)** |
| **Ground truth** | B\*07:02 | B\*07:02 | — |

**Interpretation:** ArcasHLA miscalls B\*07:06/B\*07:44 (likely misalignment of RNA reads to a similar allele). Five other tools correctly identify B\*07:02 homozygous. MajorityVote filters out ArcasHLA's error by outvoting it 5:1 — this is the key advantage of ensemble methods. The correct call is preserved.

---

## Case 3 — Majority Wrong (Systematic Allele Confusion)
**Sample: HG00096 | Gene: HLA-B | WGS**

| Tool | Allele 1 | Allele 2 | Correct? |
|---|---|---|---|
| ArcasHLA | B\*07:44 | B\*08:01 | ❌ |
| HLA-HD | B\*08:01 | B\*44:156 | ❌ (B\*44:156 ≠ B\*44:03) |
| OptiType | B\*08:01 | B\*44:03 | ✅ |
| SpecHLA | B\*08:01 | B\*08:01 | ❌ (misses B\*44:03) |
| T1K | B\*07:386 | B\*08:01 | ❌ |
| **★ MajorityVote** | **B\*07:386** | **B\*08:01** | ❌ **(1/5 agree)** |
| **Ground truth** | B\*08:01 | B\*44:03 | — |

**Interpretation:** This case illustrates the core WGS limitation. The sample is heterozygous B\*08:01 / B\*44:03. All tools correctly find B\*08:01 (the more common allele) but struggle with B\*44:03 — each tool proposes a different allele for the second position. No single wrong call reaches majority, so MajorityVote selects a minority call (T1K's B\*07:386 paired with the widely-agreed B\*08:01). Only OptiType gets the exact truth pair. The difficulty: B\*44:03 is closely related to B\*44:156 and B\*07:44, and short WGS reads cannot always disambiguate at the 3rd/4th field level. **This is why WGS majority voting (38%) lags behind OptiType (50%).**

---

## Case 4 — Correct Tool Outvoted (The WGS Ensemble Problem)
**Sample: HG00099 | Gene: HLA-B | WGS**

| Tool | Allele 1 | Allele 2 | Correct? |
|---|---|---|---|
| ArcasHLA | B\*07:44 | B\*08:01 | ❌ |
| HLA-HD | B\*08:01 | B\*08:01 | ❌ (misses B\*44:02) |
| Kourami | B\*08:01 | B\*08:134 | ❌ |
| OptiType | B\*08:01 | B\*44:02 | ✅ |
| SpecHLA | B\*08:01 | B\*08:01 | ❌ (misses B\*44:02) |
| T1K | B\*08:01 | B\*44:02 | ✅ |
| **★ MajorityVote** | **B\*08:01** | **B\*08:01** | ❌ **(2/6 agree on correct)** |
| **Ground truth** | B\*08:01 | B\*44:02 | — |

**Interpretation:** OptiType and T1K correctly identify the heterozygous pair B\*08:01 / B\*44:02. However, HLA-HD and SpecHLA both call B\*08:01 / B\*08:01 (homozygous — they miss the second allele). These two tools together create a "homozygous" majority that outvotes the correct heterozygous call 2:2:1:1. MajorityVote selects the wrong homozygous call. **This directly illustrates why quality-filtered majority voting (excluding low-accuracy tools) would fix the WGS gap** — if only OptiType + T1K + HLA-HD participated, the 2-tool correct answer (B\*08:01 / B\*44:02) would win over the 1-tool alternatives.

---

## Case 5 — Perfect WES Agreement (Why MajorityVote Beats Best Tool in WES)
**Sample: HG00103 | Gene: HLA-A | WES**

| Tool | Allele 1 | Allele 2 | Correct? |
|---|---|---|---|
| HLA-HD | A\*01:01 | A\*01:01 | ✅ |
| Kourami | A\*01:01 | A\*01:01 | ✅ |
| OptiType | A\*01:01 | A\*01:01 | ✅ |
| POLYSOLVER | A\*01:01 | A\*01:01 | ✅ |
| SpecHLA | A\*01:01 | A\*01:01 | ✅ |
| T1K | A\*01:01 | A\*01:01 | ✅ |
| **★ MajorityVote** | **A\*01:01** | **A\*01:01** | ✅ **(6/6 agree)** |
| **Ground truth** | A\*01:01 | A\*01:01 | — |

**Interpretation:** In WES, the higher HLA-region coverage from exome capture means all tools converge on the correct answer simultaneously. When this happens across most samples, MajorityVote achieves slightly higher accuracy than any individual tool (92.9% vs OptiType 92.0%) because the few individual tool errors are filtered out by the majority agreement.

---

## Summary: When Does Majority Voting Help vs Hurt?

| Scenario | Agreement | MV result | Why |
|---|---|---|---|
| All tools agree + correct | High (5-6) | ✅ Correct | Independent confirmation — very high confidence |
| Majority correct, minority wrong | Medium (4-5) | ✅ Correct | Errors are filtered by the majority |
| All tools confused, no consensus | Low (1-2) | ❌ Usually wrong | No tool has clear signal; each proposes a different wrong answer |
| Best tool correct, others wrong | Low (1-2) | ❌ Outvoted | High-accuracy tool overruled by agreeing low-accuracy tools |
| All tools systematically wrong | Any | ❌ Wrong | Ensemble cannot recover from systematic bias |

**Key principle:** MajorityVote is most valuable when tool errors are *random and independent* (Case 2 — errors cancel out). It fails when errors are *systematic and shared* (Case 3/4 — tools make the same wrong call, or one tool is clearly better but outvoted).

**Modality comparison:**
- **WES/RNA-seq:** Most tools agree on most calls → errors are independent and get filtered → MV beats best individual tool ✅
- **WGS:** Tools have very different accuracy levels (8–50%) and make systematic shared errors on difficult heterozygous calls → MV gets dragged down by low-accuracy tools ⚠️
