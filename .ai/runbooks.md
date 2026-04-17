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

## Figure and manuscript consistency
- Claude Code owns interpretation, wording, and direct edits.
- Codex may critique framing, compare options, or verify evidence before Claude edits the text or figures.
