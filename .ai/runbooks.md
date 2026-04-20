# PIHLA Runbooks

## Runtime path sanity check
- Treat `/scratch/project_2008084/pihla-publish` as the execution truth.
- If the workspace copy differs, use the runtime path for actual runs and annotate the discrepancy in the task log.

## Benchmark triage and fixes
1. Confirm whether the issue is fixture-only or real-data.
2. Check `docs/ANALYSIS_2026-04-01.md` and `docs/BENCHMARK_FIGURES.md` for scientific reporting policy.
3. Claude Code investigates and implements the fix by default.
4. If you want option analysis, result evaluation, or scientific critique before changing code, create a handoff for Codex.
5. After Codex responds, Claude Code performs the implementation follow-up.

## Review workflow
1. Run `ai-review <target> --project pihla --base-branch main` when you want Codex review prep.
2. Open the generated handoff note in `.ai/handoffs/`.
3. Use Codex `/review` guidance against that target.
4. Let Claude Code apply any follow-up code changes unless a human explicitly reassigns ownership.
5. Log the result with `ai-log`.

## Manuscript framing consultation
1. Use `pihla-publish/.ai/handoffs/codex_framing_questions.md` as the active Claude-to-Codex framing handoff.
2. Use `ai-workflows/templates/pihla-codex-response-prompt.md` as the matching prompt when Claude asks Codex for a response.
3. Capture Codex's answer in `pihla-publish/.ai/handoffs/codex_response_framing_2026-04-18.md` using the existing handoff naming convention.
4. Before any Claude synthesis step, verify the response file is no longer `pending` and that all six question sections plus `Editorial Bottom Line` and `Blocking vs Non-Blocking` contain real content.
5. If the Codex response is incomplete, record the missing sections in `.ai/handoffs/` and wait rather than synthesizing partial guidance.
6. If the Codex response is complete, use `ai-workflows/templates/pihla-claude-post-codex-synthesis-prompt.md` immediately so Claude can assess the response and make the final framing recommendation.
7. Codex should return advisory framing recommendations only; Claude remains the owner of manuscript revisions and follow-up actions.
8. If Codex returns new blocking concerns, capture them in `.ai/handoffs/` before Claude edits the manuscript.

## Figure and manuscript consistency
- Claude Code owns interpretation, wording, and direct edits.
- Codex may critique framing, compare options, or verify evidence before Claude edits the text or figures.
