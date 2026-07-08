# PIHLA Operator Checklist: Claude <-> Codex Manuscript Workflow

Use this checklist to run the full PIHLA manuscript-framing loop without guessing who should do what next.

## 1. Confirm starting state
- Work from `/scratch/project_2008084/pihla-publish`.
- Treat `/scratch/project_2008084/pihla-publish` as runtime truth for execution and scientific reporting.
- Confirm the repo-local AI files are present:
  - `.ai/context.md`
  - `.ai/runbooks.md`
  - `.ai/handoffs/codex_framing_questions.md`
  - `.ai/handoffs/codex_response_framing_2026-04-18.md`
- Confirm Claude is the owner of implementation and manuscript edits.
- Confirm Codex is being used only for second-opinion framing, critique, and decision support.

## 2. Prepare Claude's framing request to Codex
- Open `ai-workflows/templates/pihla-codex-response-prompt.md`.
- Confirm the prompt points Codex to:
  - `.ai/handoffs/codex_framing_questions.md`
  - `.ai/context.md`
  - `docs/ANALYSIS_2026-04-01.md`
  - `docs/BENCHMARK_FIGURES.md`
  - `docs/MANUSCRIPT_DRAFT_V1.md`
  - `CHANGELOG.md`
- Confirm the prompt tells Codex:
  - do not implement code
  - do not rewrite the manuscript directly
  - answer all six questions
  - keep the response advisory only

## 3. Run the Codex consultation
- Give Codex the prompt from `ai-workflows/templates/pihla-codex-response-prompt.md`.
- Tell Codex to write the response into `.ai/handoffs/codex_response_framing_2026-04-18.md`.
- Keep the response in the required structure:
  - six question sections
  - editorial bottom line
  - blocking vs non-blocking

## 4. Check Codex output quality
- Confirm every question has:
  - `Primary recommendation`
  - `Reasoning`
  - `Confidence`
  - `Suggested follow-up for Claude`
- Confirm confidence values are only `high`, `medium`, or `low`.
- Confirm Codex did not drift into implementation instructions beyond advisory next steps.
- Confirm the response distinguishes:
  - safe current claims
  - stronger but riskier claims
  - caveats
  - blocking vs non-blocking items

## 5. Run Claude's synthesis step
- Open `ai-workflows/templates/pihla-claude-post-codex-synthesis-prompt.md`.
- Give Claude both handoff files:
  - `.ai/handoffs/codex_framing_questions.md`
  - `.ai/handoffs/codex_response_framing_2026-04-18.md`
- Confirm Claude is instructed to:
  - evaluate each Codex answer
  - mark each as `agree`, `partially agree`, or `disagree`
  - produce one final manuscript-framing decision
  - identify human-only decision points

## 6. Review Claude's synthesis output
- Confirm Claude gives a decision for all six questions.
- Confirm Claude separates:
  - safe claims
  - claims to soften or remove
  - mandatory limitations/caveats
  - next manuscript edits
- Confirm Claude keeps ownership of the final conclusion rather than deferring back to Codex.

## 7. Decide what happens next
- If Claude identifies only non-blocking issues:
  - proceed with manuscript revisions
- If Claude identifies blocking scientific/framing issues:
  - capture them in `.ai/handoffs/`
  - decide whether they need:
    - more analysis
    - more validation
    - wording changes only
    - a human decision before submission

## 8. Update project state
- Update any relevant manuscript plan or status files.
- Log important decisions with `ai-log` if you are tracking them in the workflow hub.
- Keep the final Codex response and Claude synthesis in `.ai/handoffs/` so the reasoning trail stays visible.

## 9. Stop conditions
- Stop the loop when Claude has produced:
  - a final manuscript framing recommendation
  - a clear set of next edits
  - a short list of remaining human decision points
- Do not ask Codex to implement manuscript changes.
- Do not let both tools edit the same manuscript task simultaneously.

## Quick Decision Rule
- Claude builds and decides.
- Codex critiques and advises.
- Claude synthesizes and edits.
