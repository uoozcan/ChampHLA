# ChampHLA — Responses to Manuscript Comments

Point-by-point responses to the co-author/reviewer comments on `CHAMPHLA_MANUSCRIPT_V3.md`. Each item is
tagged **[ANSWER]** (conceptual question answered here), **[EDIT]** (applied to the manuscript this pass),
**[FIGURE]** (figure task, deferred to the figure pass), or **[DECISION]** (needs a choice; noted).

Last updated: 2026-07-16.

---

## 1. "…for cohorts lacking gold typing…" — jargon in the abstract
**[EDIT]** Agreed — "gold typing" is undefined jargon and the abstract has no room to define it. Reworded to
**"cohorts without reference-quality HLA genotypes"** (and, where the fuller phrase appears, "cohorts lacking
established/gold-standard HLA reference genotypes"). "Gold-standard" is retained only where it is contrasted
explicitly with the silver-standard candidate labels, because there the contrast carries the meaning.

## 2. "…CSC Puhti supercomputer…" — Finland-specific
**[EDIT]** Agreed. Generalised to **"HPC clusters"** with CSC/SLURM given only as *an example* of a provided
profile: "…including profiles for Docker-based local execution and SLURM-based HPC clusters (a ready-made
CSC/Puhti profile is bundled as an example)." This keeps the concrete artefact for reproducibility without
implying the tool is Finland-only.

## 3. "…most polymorphic locus in the entire human genome…" — needs citation
**[EDIT/ANSWER]** Added citation **Trowsdale & Knight (2013), *Annu. Rev. Genomics Hum. Genet.* 14:301–323**
(MHC genomics and human disease — a standard, correct source for HLA being the most polymorphic region).
Robinson et al. (2020, IMGT/HLA) is also already cited nearby for the allele count.

## 4. "…matched for at least HLA-A, -B, -C, and -DRB1 (8/8 match) or preferably also HLA-DQB1 (10/10)…" — explain the counts
**[ANSWER + EDIT]** The numbers are **allele-match counts**: each locus contributes two alleles (one per
chromosome). Matching at four loci — HLA-A, -B, -C, -DRB1 — gives **4 loci × 2 alleles = 8 possible matches
("8/8")**; adding HLA-DQB1 as a fifth locus gives **5 × 2 = 10 ("10/10")**. This is standard haematopoietic
stem-cell-transplant donor-matching nomenclature. Added a short parenthetical to the sentence: *"(each locus
contributes two alleles, so four loci give eight possible matches and five loci give ten)."*

## 5. "…the HLA locus poses unique computational challenges…" — needs citation
**[EDIT/ANSWER]** Added **Kiyotani, Mai & Nakamura (2017), *J. Hum. Genet.* 62:397–405** (comparison of
NGS-based HLA-typing tools; documents the alignment/paralogue challenges). *Please verify the exact
volume/pages against your reference manager;* an equally suitable alternative is a recent HLA-caller benchmark
review if you prefer.

## 6. "…allele-specific expression … miss the minor allele in low-expression conditions…" — needs citation
**[EDIT/ANSWER]** Added **Aguiar, Tishkoff, Meyer et al. (2019), *PLoS Genet.* 15:e1008091** ("Expression
estimation and eQTL mapping for HLA genes…"), which directly documents allele-specific expression at HLA loci.
Orenbuch et al. (2019, arcasHLA) is also already cited and discusses expression-driven RNA calling.

## 7. "SpecHLA combines variant-aware alignment with haplotype phasing…" — citation for SpecHLA
**[EDIT — needs your confirmation]** SpecHLA is **not currently in the reference list**. I have added an inline
citation and a References stub, but I could not verify the exact reference from the repository, so **please
insert the correct SpecHLA citation** (the tool is `deepomicslab/SpecHLA`; the accompanying paper). I have
marked the References entry `SpecHLA — [citation to be confirmed by author]` so it is impossible to miss.

## 8. The three Methods paragraphs (champion-challenger fallback + nested CV + reproducibility/sweep)
**[EDIT — reviewed and revised]** Your proposed text is accurate and clearer than the current version, and the
grid values match the code exactly (`SUPPORT_GRID = [0.20, 0.35, 0.50, 0.65]`,
`MARGIN_GRID = [0.0, 0.05, 0.10, 0.20]`, `TOOLS_GRID = [1, 2, 3]`). Three corrections were made while
integrating it:
- **"four decision thresholds are learned"** → **three** thresholds are grid-searched (support fraction,
  weight margin, minimum supporting tools); the fourth gate (challenger allele must be unambiguous) is a fixed
  requirement, not a searched parameter. Reworded accordingly.
- The decision-trace identifiers (`champion_retained`, etc.) were converted to plain language per comment 10.
- Kept your excellent in-sample-sweep caveat, and made the weighting description consistent with the corrected
  Methods §3 (base-reliability weights are a **fixed full-cohort prior**, not re-estimated per fold — only the
  champions and thresholds are out-of-fold; see the earlier weight-provenance correction).

## 9. "…normalized to [0, 1]…" — which normalization?
**[ANSWER + EDIT]** It is **neither min–max nor z-score**. Each tool's native confidence proxy is first
**clipped to the unit interval [0,1]** (`clip_unit` in `normalize_confidence_values`), then mapped to a
**calibrated probability by per-tool Platt scaling** (logistic calibration fitted on the benchmark,
`fit_platt_calibrator`), and finally gated by the Brier/ECE guardrail. The manuscript now states this
explicitly rather than saying "normalized to [0,1]".

## 10. "…pass the guardrail, effective_confidence is the…" — no `_` identifiers in the main text
**[EDIT]** Agreed. Decision-trace and weight identifiers written in `snake_case` (e.g. *champion_retained,
champion_missing, challenger_override, effective_confidence, base_reliability, final_weight*) are converted to
plain-language terms in the running text ("champion retained", "effective confidence", etc.). Genuine artefact
and file names (`benchmark_metadata.json`, `tool_weights_{wgs,wes,rna}.json`, the `.tsv`/`.json` sweep files,
`bin/*.py`) are kept in monospace because they are concrete reproducibility references, not jargon.

## 11. "final_weight = 0.7 × base reliability + 0.3 × effective confidence" — proper equation + rationale
**[ANSWER + EDIT]** Rendered as a numbered display equation:

> w_final = 0.7 · r_base + 0.3 · c_eff

**Rationale (added to the text):** the base reliability r_base is a *directly measured* quantity (each tool's
benchmark correct-call rate), so it is the primary, more trustworthy signal and receives the larger weight;
the effective confidence c_eff is a *derived, noisier, and often miscalibrated* proxy (which is exactly why it
is passed through the Brier/ECE guardrail), so it is down-weighted. The 0.7/0.3 split is a **fixed a-priori
default** (`weight_alpha`/`weight_beta`), deliberately *not* tuned, to avoid adding another optimised
hyperparameter; the routing-vs-weighting ablation and the new weighting-sensitivity control (Table 5) show the
reported accuracy is essentially insensitive to this choice (the significant WGS result is weight-independent).

## 12. "A 50-sample 'wave2' subset…" — reviewers will ask where wave1 is
**[ANSWER + EDIT]** Valid concern: **wave1 does exist** in the pipeline (`slurm_benchmark_wgs_wave1.sh`,
`conf/benchmark_1000g_wgs_wave1.yaml`) as an earlier WGS batch, but it is not part of the manuscript's
narrative, so "wave2" invites the question. Renamed throughout to a neutral, self-explanatory term: **"an
initial 50-sample WGS calibration subset."** It is used only for the confidence-calibration-guardrail
characterisation (Table 5), which the text already states.

## 13. "CEU (n=12), FIN (n=7), GBR (n=6), TSI (n=6), YRI (n=7)…" — expand abbreviations
**[EDIT + ANSWER]** These are 1000 Genomes population codes; expanded on first use: **CEU** (Utah residents of
Northern & Western European ancestry), **FIN** (Finnish in Finland), **GBR** (British in England & Scotland),
**TSI** (Toscani/Tuscans in Italy), **YRI** (Yoruba in Ibadan, Nigeria).

## 14. "A supplementary matched-subject trimodal robustness cohort of 106 samples" — where from?
**[ANSWER]** It is defined in Methods §6 (and Table 2, "Matched intersection"): the **intersection of samples
that have WGS, WES, and RNA-seq outputs and a ground-truth HLA call** — i.e. the subset of subjects typed on
all three modalities. It was only unexplained in the abstract; a three-word signpost
("WGS∩WES∩RNA intersection") was added there pointing to Methods §6.

## 15. "…propagated through `benchmark_metadata.json`…" — reword "propagated"
**[EDIT]** Changed to **"recorded in `benchmark_metadata.json`"**.

## 16. "Per-tool wall-clock time" — what is this?
**[ANSWER + EDIT]** **Wall-clock time** is the elapsed real time from the start to the end of a tool's run
("time on the wall"), as opposed to CPU time (which sums time across cores and can exceed wall-clock for
multi-threaded tools). It is read from the Nextflow execution trace `realtime` field. Added a one-line
definition on first use in Methods §7.

## 17. "Use only 2 digits after decimal for all numbers."
**[DECISION — recorded, applied in the figure/rounding pass]** You chose to **keep the key held-out
correct-call rates and their confidence intervals at 3–4 digits** (so e.g. WES 0.9436 vs 0.9359 stay distinct)
and round the *other* numbers (timings, RAM, CPU%, percentages) to 2 decimals. This global rounding is grouped
with the figure pass, not applied in this text-edit pass.

## 18. "Please use Table numbers for all tables."
**[EDIT]** The three panels I had appended under "Table 3" (statistical summary, weighting sensitivity,
override audit) are given **their own numbers — Tables 4, 5, 6** — and the subsequent tables renumbered
(old Tables 4–10 → 7–13), with every in-text reference updated.

## 19. "Put all mentioned figures into the docx document."
**[FIGURE — clarification needed]** All 14 supplementary figures **are already embedded** in the separate
`CHAMPHLA_MANUSCRIPT_V3_SUPPLEMENTARY.docx`, and the 7 main figures in `CHAMPHLA_MANUSCRIPT_V3.docx`. If you
instead want a **single combined `.docx`** containing main + supplementary figures, tell me and I will build
it — that is a small change to the docx builder. (Deferred with the figure pass.)

## 20. Supplementary Figure S11 (computational resources) — legend inside 3rd subplot, larger fonts, left-aligned titles
**[FIGURE — deferred]** Will regenerate `figure_8_computational_performance` with: legend moved inside the
third subplot (upper-right), larger tick/label fonts, and left-aligned panel titles rendered as
"1) Peak RAM (GB)", "2) CPU utilization (%)", "3) Wall-clock time (h)".

## 21. "For figure 13, bar sizes should be consistent in all subplots."
**[FIGURE — clarification needed]** After the figure consolidation there is no main Figure 13. The source stem
`figure_13` is the **NCI-60 external-validation figure, now Supplementary Figure S14** (panels A/B/C of bars).
I will make the bar widths consistent across its subplots — **please confirm** this is the figure you meant.

## 22. "Make all figures with white background and Arial fonts."
**[FIGURE — deferred]** Will set a global matplotlib style (figure & axes `facecolor='white'`,
`font.family='Arial'`/sans-serif fallback) in the shared plot style and regenerate all figures.

---

### Summary of what changed this pass
Applied as text edits (comments 1–16, 18): accessibility wording, four citations (one flagged for you to
confirm — SpecHLA), the 8/8–10/10 explanation, the revised three Methods paragraphs with a proper weight
equation and rationale, the clip-to-[0,1]+Platt normalization statement, removal of `snake_case` jargon from
prose, the wave2 rename, expanded population codes, the 106-cohort signpost, "propagated"→"recorded", the
wall-clock definition, and full table renumbering.

Deferred to the figure/rounding pass (comments 17, 19–22): the 2-decimal policy, Suppl. Fig S11 and NCI-60
figure edits, white-background/Arial across all figures, and (if wanted) a single combined docx. Two items
need a one-word confirmation from you: the SpecHLA citation (comment 7) and the "figure 13" identity
(comment 21).
