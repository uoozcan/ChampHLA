#!/usr/bin/env python3
"""
Inserts informal 1-2 paragraph explanations after each figure caption in the HTML report.
Two caption types handled:
  1. <div class="fig-caption">...</div>  (main figures 1-14)
  2. <figcaption>...</figcaption>        (trimodal + VENEX figures)

Usage:
  python3 add_figure_explanations.py \
      --input  /scratch/project_2008084/mvhla_report_v13.html \
      --output /scratch/project_2008084/mvhla_report_v14.html
"""
import argparse, re
from pathlib import Path

# ── Explanation CSS ───────────────────────────────────────────────────────────
EXPLAIN_CSS = """
  .fig-explain {
    background: #f8fbff;
    border-left: 3px solid #4fc3f7;
    padding: 10px 16px;
    margin: 8px 0 20px 0;
    font-size: 0.93em;
    color: #2a3a50;
    line-height: 1.65;
    border-radius: 0 6px 6px 0;
  }
  .fig-explain p { margin: 0 0 8px 0; }
  .fig-explain p:last-child { margin-bottom: 0; }
"""

# ── Explanations keyed by unique substring from figure title ─────────────────
EXPLANATIONS = {

"Figure 1 — Tool Accuracy Overview": """
<p>Think of this chart as a report card for every HLA typing tool in the benchmark, shown
side by side. Each dot on the horizontal bar marks how often a tool returned the correct
two-allele answer for a patient sample — a higher dot means better performance. The thin
lines through each dot show the range of uncertainty (if the line is long, the tool was
evaluated on fewer samples).</p>
<p>The most striking pattern is the gap between WGS (top panel) and the other two
modalities. In WGS, the best tool reaches only about 50% — meaning it gets the right
answer for roughly 1 in 2 samples. In WES and RNA-seq, most tools cluster above 90%. This
is not a software quality issue: WGS simply provides a harder input for HLA typing because
reads come from the whole genome and are harder to assign to the correct allele.</p>
""",

"Figure 9 — Cross-Modality Tool Performance Profile": """
<p>This figure traces each tool like a thread running from left (WGS) through the middle
(WES) to the right (RNA-seq). A rising line means the tool improves as the input gets more
targeted; a falling line means the tool is best suited to the modality on the left. The
dramatic rise of ArcasHLA — barely above zero in WGS but shooting up to >90% in RNA-seq —
makes visible what the text explains: this tool was designed specifically for RNA data.</p>
<p>The key practical message is that tool choice should always match the sequencing type.
Running an RNA-seq specialist tool on WGS data produces valid-looking output but incorrect
results, and vice versa. This figure is a diagnostic you can apply to any new HLA typing
tool: run it across modalities and watch whether the accuracy line makes biological sense.</p>
""",

"Figure 2 — Ensemble Methods vs. Best Single Tool": """
<p>This bar chart answers the question: does combining multiple tools actually help? Each
cluster of three bars represents one modality. The light-blue bar is the best tool in that
modality on its own. The red bar is what happens when you give every tool an equal vote
(MajorityVote). The dark-red bar is what happens when you weight votes by each tool's
proven reliability (WeightedConsensus).</p>
<p>In WGS, the red MajorityVote bar is <em>shorter</em> than the best-single-tool bar —
combining tools actually hurt. This is because five poorly-performing tools outvote the one
good tool. The dark-red WeightedConsensus bar recovers some accuracy by down-weighting
those bad tools, but it also chooses to skip the hardest loci entirely (which is why its
overall bar is the shortest — it abstains on 38% of WGS loci). In WES and RNA-seq, all
three bars are nearly the same height, confirming that when tools already agree well,
ensemble voting adds little extra accuracy but does add a quality-control safety net.</p>
""",

"Figure 10 — Callable Rate vs. Accuracy Trade-off": """
<p>Each bubble in this chart is one method in one modality. The further right a bubble is,
the more often that method produces a call (rather than saying "I'm not sure"). The higher
a bubble is, the more accurate those calls are. A perfect method would sit in the top-right
corner: 100% of loci called, 100% correctly. The bubble size reflects sample count.</p>
<p>The WGS bubbles cluster in the lower-left — low callable rate, lower accuracy — while
WES and RNA-seq bubbles cluster top-right. The split between MajorityVote (always calls,
lower accuracy) and WeightedConsensus (sometimes abstains, higher accuracy) is visible
within each modality as two bubbles at the same height but different horizontal positions.
This chart is useful for deciding which method to use: if you need every sample typed
regardless of confidence, move right; if you need the highest accuracy on those samples
you do type, move up.</p>
""",

"Figure 3 — Per-Gene Accuracy (HLA-A, -B, -C)": """
<p>HLA typing accuracy is not uniform across the three class I genes — and this chart shows
exactly how much it varies. Each cluster of bars covers HLA-A, HLA-B, and HLA-C side by
side for a given modality and method. If the bars were the same height, the genes would be
equally easy to type; the height differences tell you where the hard problems are.</p>
<p>In WGS, HLA-C is consistently the shortest bar: it is the hardest gene to type correctly
from whole-genome reads because of its proximity to pseudogenes that share similar sequences.
HLA-B tends to be easiest in WGS. In WES and RNA-seq the gene-to-gene differences are much
smaller, because the enrichment and expression steps already filter out many of the confusing
reads. If your study relies on accurate HLA-C typing from WGS, this figure explains why
you should expect more errors and plan additional validation accordingly.</p>
""",

"Figure 5 — Abstention–Accuracy Trade-off": """
<p>Imagine tuning a dial that controls how confident the ensemble needs to be before it
commits to an answer. Turn the dial to the left (low threshold): the method calls almost
every locus, but some of those calls are wrong. Turn it to the right (high threshold): the
method skips more loci but is more accurate on the ones it does call. This figure plots
both effects simultaneously for all three modalities.</p>
<p>The blue line (left axis) shows the callable rate falling as the threshold rises — fewer
loci get a committed answer. The coloured lines (right axis) show accuracy-among-called
rising as the threshold rises — the loci that do get called are increasingly trustworthy.
The default threshold of 0.55 sits near the "knee" of the curve: the point where accuracy
starts to plateau while callable rate is still reasonably high. The curves also show that
in RNA-seq, even very aggressive thresholds keep accuracy above 95%, while in WGS, accuracy
never exceeds ~65% even at the strictest settings — a hard ceiling imposed by tool quality.</p>
""",

"Figure 11 — WeightedConsensus vs. MajorityVote Accuracy Delta": """
<p>This diverging bar chart directly answers: does switching from equal-weight voting to
reliability-weighted voting help or hurt for each gene and modality? Bars pointing right
mean WeightedConsensus outperforms MajorityVote; bars pointing left mean it underperforms.
Note that WeightedConsensus abstains on some loci — so "outperforms" here means it is more
accurate on the loci it <em>does</em> commit to.</p>
<p>In WES and RNA-seq, the differences are small in either direction — the two methods are
nearly interchangeable because tools already agree well and the weighting has little effect.
In WGS, WeightedConsensus shows the most positive bars: by down-weighting the weak tools
and abstaining on conflicted loci, it recovers meaningful accuracy over the naive majority
vote. The negative bars in some WGS genes reflect loci where abstention reduces the overall
callable count without fully compensating in accuracy — a known limitation when the entire
tool pool performs poorly.</p>
""",

"Figure 12 — Tool Agreement Predicts Call Accuracy": """
<p>This chart asks: when tools agree with each other, does that agreement actually mean the
answer is more likely to be correct? The x-axis shows the "support fraction" — how much of
the total tool weight was behind the winning allele pair. The y-axis shows the fraction of
such calls that turned out to be correct. If tool agreement were meaningless, all the dots
would be at the same height regardless of x-position.</p>
<p>Instead, the dots rise from left to right: loci where only a slim majority of tools
agreed are correct less often, while loci where almost all tool weight lined up behind one
answer are nearly always correct. This validates the central premise of WeightedConsensus —
the support fraction is a genuine confidence signal, not just a voting artefact. Practically,
this means you can use the support fraction in the output files as a per-locus quality score:
calls with support fraction > 0.8 are extremely reliable; calls near the 0.55 threshold
deserve more scrutiny.</p>
""",

"Figure 13 — Per-Gene Accuracy Gain (WeightedConsensus vs. Baselines)": """
<p>This figure zooms into the per-gene improvement that WeightedConsensus achieves over two
baselines: MajorityVote (equal-weight vote of all tools) and the single best tool available
in each modality. Positive bars mean WeightedConsensus wins; values near zero mean it ties.
The figure is organised by gene (A, B, C) and coloured by modality.</p>
<p>The pattern reinforces the overall story: WGS benefits most from the weighted approach
(especially for HLA-B, where down-weighting weak tools rescues meaningful accuracy), while
WES and RNA-seq show modest gains because the best single tools and MajorityVote are already
near-optimal. For a practitioner deciding whether to run WeightedConsensus instead of just
using the best available tool, this figure provides the per-gene evidence: in WGS it clearly
pays off; in RNA-seq the benefit is marginal but never harmful.</p>
""",

"Figure 4 — Bimodal WES+RNA Joint Consensus": """
<p>This figure addresses a practical question: if a laboratory has both WES and RNA-seq data
for the same patient, should it combine the two for HLA typing? Six methods are compared
side by side for each of the three class I genes. The lighter bars are WES-only and RNA-only
baselines; the bold-outlined bars are what happens when both data types are pooled into a
single consensus vote.</p>
<p>The bimodal bars sit consistently at or above the best single-modality bars for every
gene. The gain is modest — about one percentage point per gene — but it is consistent and
free of charge: no additional sequencing is needed, just pooling the calls that already
exist. The practical upshot is clear: whenever both WES and RNA-seq data are available (as
they routinely are in cancer genomics workflows), running the bimodal consensus should be
the default choice over either modality alone.</p>
""",

"Figure 8 — LOH / Homozygosity False-Positive Rate": """
<p>A key quality metric for any HLA typing tool is whether it can tell the difference
between a patient who genuinely carries two copies of the same allele (true homozygosity)
and a patient who carries two different alleles but the tool only found one (allele dropout).
This figure measures that by computing the "false duplicate rate": out of all calls where a
tool reported both alleles as identical, what fraction of those cases was actually
heterozygous in ground truth?</p>
<p>A rate near 100% means the tool almost never correctly identifies true homozygosity — it
calls duplicates only when something went wrong technically. Rates near 0% mean the tool's
homozygous calls are trustworthy. The dramatic modality effect visible in this figure — WGS
tools clustering at very high false-duplicate rates while WES and RNA-seq tools cluster
lower — explains why allele dropout is considered a WGS-specific problem. If you need to
detect HLA Loss of Heterozygosity in tumour samples, this chart tells you which tool-modality
combinations can be trusted for that task.</p>
""",

"Figure 14 — LOH Status Distribution in VENEX-WGS Cohort": """
<p>This figure applies the LOH quality-control framework to the VENEX cancer cohort — a
real clinical dataset rather than the 1000 Genomes reference samples. Each horizontal bar
represents one HLA gene (A, B, or C) and shows how the 80 SpecHLA-typed samples were
classified. Green means the allele frequencies looked balanced and both alleles were
confidently detected. Red means only one allele was called — a potential LOH flag. Purple
means there were not enough heterozygous markers in the reads to make a reliable call
either way.</p>
<p>The dominance of purple ("insufficient evidence") is the most important result here.
It tells us that even at standard 30× WGS coverage, the HLA region does not consistently
provide enough signal for confident allele-balance QC in most samples. This motivates two
practical recommendations: (1) for high-stakes LOH calls (e.g., before neoantigen vaccine
design), increase WGS depth or add a secondary platform; (2) treat the red "homozygous
needs review" flags not as confirmed LOH but as candidates for orthogonal validation.</p>
""",

"Figure 6 — Failure Mode Profile by Tool and Modality": """
<p>When a tool gets an HLA call wrong, it does not always fail in the same way. This chart
breaks down incorrect calls into categories: wrong at the second field (the tool found the
right allele family but got the specific variant wrong), wrong at the first field (the tool
called a completely different allele family), abstained (the tool produced no call), or
ambiguous (multiple alleles scored equally). Each horizontal bar is one tool in one modality;
the colour segments show the mix of failure types.</p>
<p>The failure mode profile matters because different error types have different consequences
in clinical practice. A "wrong field 2" error (e.g., A*02:01 when truth is A*02:06) may not
affect transplant compatibility if the alleles share epitopes; a "wrong field 1" error
(e.g., A*02 when truth is A*24) is a serious mismatch. Tools with a high proportion of
field-1 errors in WGS are fundamentally misidentifying allele families — a problem no
ensemble strategy can fully correct. This chart helps prioritise which tools to include in
the ensemble and which to flag as unreliable for specific modalities.</p>
""",

"Figure 7 — Discordance Event Taxonomy": """
<p>Even when a consensus call has been made, the tools contributing to it may have
disagreed in structured ways. This figure counts those disagreements by category across all
three modalities. "Technical conflict" means tools disagreed without a clear biological
explanation — the most common category in WGS because tool accuracy is so variable.
"Possible expression bias" flags cases in RNA-seq where one allele may have been
transcribed more than the other, causing the expressed allele to dominate tool calls.
"DNA-RNA discordance" marks the rarest but most interesting events: DNA-based and
RNA-based tools systematically disagree for the same patient.</p>
<p>The 12-fold decrease in low-evidence conflicts from WGS (176) to RNA-seq (14) is the
most striking quantitative pattern. It reflects the underlying tool reliability: when
individual tools each score 90%+ accuracy, they naturally converge on the same answer.
When tools score 12–50%, each independently chooses a different wrong allele, maximising
disagreement. For users reading discordance tags in the PIHLA output files: a
"technical_conflict" tag in WGS is almost expected and does not necessarily mean the
consensus call is wrong; a "dna_rna_discordance" tag anywhere should prompt immediate
manual review.</p>
""",

# ── Trimodal (figcaption) ─────────────────────────────────────────────────────
"Figure 15. Trimodal": """
<p>This chart extends the bimodal comparison by adding WGS data to the mix. The lighter bars
repeat the per-modality WGS, WES, and RNA-seq baselines restricted to the 106 samples that
have all three data types available. The bold-outlined bars show what happens when all three
modalities are pooled into a single vote. The purple bar (Trimodal WeightedConsensus) is
the ultimate combination: every tool from every sequencing type contributing to one answer.</p>
<p>The key takeaway is that the purple bar barely rises above the dark-red bimodal bar for
any of the three genes. Adding WGS data to an already well-calibrated WES+RNA ensemble
provides essentially no benefit — because the WGS tools introduce noise at a level that
the weighting system can partially suppress but not eliminate. This is a practical finding
with direct cost implications: if you already have WES and RNA-seq, you do not need to
generate WGS data for the purpose of improving HLA typing accuracy.</p>
""",

# ── VENEX figures (figcaption) ────────────────────────────────────────────────
"VENEX Figure 1. Tool call coverage": """
<p>Before comparing tool accuracy, it is important to know which tools actually produced
usable results. This chart shows the coverage gaps: for the 87 VENEX WGS samples, only
SpecHLA returned calls for more than 90% of samples, while OptiType processed fewer than
40%. These coverage differences arise from practical issues — memory limits, convergence
failures, input format requirements — rather than tool design. In a production workflow,
missing calls mean missing data for those patients, so coverage rate is as important a
metric as accuracy.</p>
<p>Notice also that only SpecHLA, HLA-HD, and T1K provide class II calls (DRB1, DQB1)
for this cohort; OptiType was designed for class I only. If class II typing is required for
clinical decisions (transplantation, autoimmune disease risk), the choice of tool is
effectively made for you by this coverage constraint.</p>
""",

"VENEX Figure 2. Inter-tool allele-pair agreement": """
<p>Without external ground truth, the best internal quality check is whether multiple
independent tools agree. This tile grid shows, for each pair of tools that typed the same
samples, what fraction of those shared samples received the same allele pair from both
tools. A green tile (high agreement) means two tools independently reached the same
conclusion; a red tile means they frequently disagree.</p>
<p>High agreement between tools does not guarantee correctness — two tools can consistently
make the same mistake — but it does increase confidence. In practice, samples where three
or more tools agree deserve high confidence; samples where two tools calling from different
algorithmic approaches agree are also generally reliable. The off-diagonal tiles that remain
red or yellow identify the tool pairs most likely to generate conflicting calls, and flagging
those conflicts in the output is exactly what the PIHLA discordance taxonomy is designed
to do.</p>
""",

"VENEX Figure 3. Homozygous call rate": """
<p>When an HLA typing tool reports that both copies of a gene carry the exact same allele,
that can mean one of two things: the patient genuinely inherited the same allele from both
parents (true homozygosity, which is uncommon at highly polymorphic HLA loci), or the tool
failed to find the second allele and defaulted to reporting the first one twice (allele
dropout). This chart shows how often each tool makes homozygous calls for each gene in the
VENEX cohort.</p>
<p>The 1000 Genomes benchmark (where we know the truth) tells us that 61–91% of WGS
homozygous calls are dropout artefacts. Applying that calibration to the VENEX numbers
means that most of the homozygous calls visible here are likely technical failures rather
than genuine biology. The dashed red line at 10% marks an approximate expected rate if
dropout were absent. Any tool-gene combination exceeding this line should be treated with
caution, particularly in clinical settings where a false homozygous call could mask a
medically important second allele.</p>
""",

"VENEX Figure 4. LOH status distribution": """
<p>This stacked bar chart summarises the result of asking, for each patient and each HLA
gene: does the sequencing data support the presence of both alleles, or does it look like
one allele is missing? The green segment (balanced heterozygous) is the cleanest result —
both alleles detected at roughly equal read frequencies, suggesting no LOH. The red segment
(homozygous needs review) flags cases where the call was homozygous and the allele-balance
QC could not rule out dropout or true LOH. The purple segment (insufficient evidence) is
the honest "I don't know" category.</p>
<p>The dominance of purple across all three genes is itself informative: it tells you that
standard 30× WGS does not reliably provide enough heterozygous SNPs in the HLA region to
power allele-balance classification for most samples. This is not unique to PIHLA — it is
a general limitation of short-read WGS at the MHC. The 15 red-flagged (sample, gene) pairs
represent the highest-priority candidates for follow-up, but even these should be confirmed
with orthogonal methods (SNP arrays, deeper targeted sequencing) before being reported
clinically as LOH-positive.</p>
""",

"VENEX Figure 5. Number of unique 2-field alleles": """
<p>This chart measures something different from accuracy: how diverse are the alleles that
each tool reports across all patients? A tool that calls many different alleles is navigating
the full breadth of HLA diversity; a tool that repeatedly calls the same few alleles may
have a limited reference database or a tendency to collapse rare alleles into common ones.</p>
<p>Allele diversity also serves as a sanity check on the cohort itself. If the VENEX cohort
were highly homogeneous (e.g., all from the same ethnic group), you would expect to see
fewer unique alleles even if tools were working perfectly. The diversity numbers here —
spanning dozens of unique alleles per gene — suggest that the cohort covers a reasonable
range of HLA variation, and that the tools are not systematically collapsing rare alleles.
SpecHLA shows the highest diversity mainly because it typed the most samples, not because
it discovers more alleles per sample than other tools.</p>
""",

}

def find_explanation(title_snippet):
    for key, exp in EXPLANATIONS.items():
        if key in title_snippet:
            return exp.strip()
    return None

def build_explain_div(text):
    paragraphs = [p.strip() for p in text.strip().split("\n\n") if p.strip()]
    inner = "\n".join(f"<p>{p}</p>" for p in paragraphs)
    return f'\n<div class="fig-explain">\n{inner}\n</div>\n'

def patch(html):
    # 1. Inject CSS into <style> block
    html = html.replace("</style>", EXPLAIN_CSS + "\n</style>", 1)

    # 2. Patch <div class="fig-caption"> blocks (figures 1-14)
    def replace_fig_caption(m):
        full = m.group(0)
        # Extract title from <strong>...</strong>
        title_m = re.search(r'<strong>(Figure[^<]*)</strong>', full)
        if not title_m:
            return full
        title = title_m.group(1)
        exp = find_explanation(title)
        if not exp:
            return full
        return full + build_explain_div(exp)

    html = re.sub(
        r'<div class="fig-caption">.*?</div>',
        replace_fig_caption,
        html,
        flags=re.DOTALL
    )

    # 3. Patch <figcaption> blocks (trimodal + VENEX figures)
    def replace_figcaption(m):
        full = m.group(0)
        inner_m = re.search(r'<figcaption[^>]*>(.*?)</figcaption>', full, re.DOTALL)
        if not inner_m:
            return full
        inner_text = re.sub(r'<[^>]+>', '', inner_m.group(1))
        # Find matching explanation by checking start of cleaned text
        exp = None
        for key, text in EXPLANATIONS.items():
            # check key against inner_text
            key_clean = re.sub(r'[^\w\s]', '', key).lower()
            inner_clean = re.sub(r'[^\w\s]', '', inner_text[:len(key)+20]).lower()
            if key_clean[:30] in inner_clean:
                exp = text.strip()
                break
        if not exp:
            return full
        return full + build_explain_div(exp)

    html = re.sub(
        r'<figure[^>]*>.*?</figure>',
        replace_figcaption,
        html,
        flags=re.DOTALL
    )

    return html

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input",  default="/scratch/project_2008084/mvhla_report_v13.html")
    p.add_argument("--output", default="/scratch/project_2008084/mvhla_report_v14.html")
    args = p.parse_args()

    print(f"Reading {args.input} ...")
    html = Path(args.input).read_text(encoding="utf-8")
    print(f"  {len(html):,} chars")

    print("Inserting figure explanations...")
    patched = patch(html)

    # Count how many explanations were inserted
    n = patched.count('class="fig-explain"')
    print(f"  Inserted {n} explanation blocks")

    Path(args.output).write_text(patched, encoding="utf-8")
    size_kb = Path(args.output).stat().st_size // 1024
    print(f"  Written → {args.output}  ({size_kb} KB)")

if __name__ == "__main__":
    main()
