# ChampHLA-RefFormer publication assessment

Decision: **benchmark_software_paper; report_refformer_and_metaconsensus_as_ablations**.

## Evidence

- Candidate reranker prespecified gate passed: False.
- Pooled IMGT delta versus no-IMGT: -0.0053.
- Logistic MetaConsensus internal gate passed: True.
- Population-stratified MetaConsensus WGS: 0.5353 versus 0.5012; delta +0.0341, cluster 95% CI [+0.0097, +0.0584].
- Synthetic alignment-aware held-out accuracy: 0.0702 (random expectation 0.1250).
- Full fold-0 imgt real pilot: top-1 0.0078, top-5 0.0233, missing tensors 0.
- Full fold-0 no_imgt real pilot: top-1 0.0000, top-5 0.0078, missing tensors 0.
- Supplemental frozen MetaConsensus nci60_rna: 0.7927 versus 0.8171; delta -0.0244, 2-point non-regression passed: False.
- Supplemental frozen MetaConsensus nci60_wes: 0.8171 versus 0.8415; delta -0.0244, 2-point non-regression passed: False.
- Qualifying external validation complete: False.
- The planned MetaConsensus external gate is incomplete: HPRC WGS n=1, verified IHWG WES availability n=2, and E-MTAB-197 has only 35 development-nonoverlapping subjects rather than the target 50.
- Separate frozen-CC benchmark evidence is more substantial in NCI-60 (n=42 RNA and n=42 WES), with small GIAB, HPRC, and IHWG pilots. It supports a benchmark/software paper but is not external validation of the newly frozen MetaConsensus model.
- Original source integrity: 243/243 recorded files unchanged.

## Manuscript consequence

Neither RefFormer nor MetaConsensus should be the method headline. Report both transparently as ablations: MetaConsensus passed internally in WGS but missed the two-point non-regression criterion in both available independent NCI-60 modalities. Center the publishable paper on reproducible multi-caller harmonization, auditability, modality-specific benchmark behavior, candidate-set ceilings, runtime, and the finding that learned selection did not reliably improve saturated WES/RNA consensus.

Do not claim completed external validation from the existing pilot-sized folders, and do not spend the planned 10x5 nested compute on the standalone model after this feasibility failure.
