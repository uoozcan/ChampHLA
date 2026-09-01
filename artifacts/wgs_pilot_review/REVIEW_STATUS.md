# Two-subject full-CRAM WGS review status

The bounded HG00096/HG00097 pilot is technically complete. The collected
bundle contains 30/30 expected callable records: two subjects, three class-I
genes, and five intended-use callers. Its freeze validation has no failures.

On 2026-09-01, Codex independently inspected the native A/B/C summary lines for
HLA-HD, Kourami, OptiType, SpecHLA, and T1K for both subjects. The inspection
confirmed all 30 canonical unordered two-field pairs recorded in
`harmonized_calls.tsv`. It specifically checked high-field reduction,
homozygous single-allele/dash output, gene prefixes, and allele-order changes.
The separate automated audit likewise reports 30 passes and zero mismatches.

This is an assistant-supported technical review, not named human adjudication.
Accordingly, the generated manual-review packet remains `pending`; it has not
been relabelled as passed. A supervisor or designated analyst must sign the 30
rows. The production rule also requires ten reviewed records per caller (50
total), so 20 additional records from valid full-CRAM runs are required before
the WGS expansion is authorized.

No truth was used for caller collection, parsing, or this audit.
