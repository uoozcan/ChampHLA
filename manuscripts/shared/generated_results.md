# Registered ChampHLA results

| Result | Cohort | Modality | Method | Evidence role | Validity | Correct / loci | Accuracy | Abstract |
|---|---|---|---|---|---|---:|---:|---|
| DEV_WES_CONSENSUS | 1000G-development | wes | SimplePluralityLex | development | valid | 365 / 390 | 0.935897 | no |
| DEV_RNA_CONSENSUS | 1000G-development | rnaseq | SimplePluralityLex | development | valid | 306 / 321 | 0.953271 | no |
| DEV_WGS_CONSENSUS | 1000G-development | wgs | SimplePluralityLex | invalid | invalid | 164 / 411 | 0.399027 | no |
| DEV_TT_WES_BASE | 1000G-development | wes | SimpleTwoThirdsConsensus | development | valid | 328 / 390 | 0.841026 | no |
| DEV_TT_WES_GCC | 1000G-development | wes | TwoThirdsGuardedCC | development | valid | 369 / 390 | 0.946154 | no |
| DEV_TT_RNA_BASE | 1000G-development | rnaseq | SimpleTwoThirdsConsensus | development | valid | 296 / 321 | 0.922118 | no |
| DEV_TT_RNA_GCC | 1000G-development | rnaseq | TwoThirdsGuardedCC | development | valid | 306 / 321 | 0.953271 | no |
| DEV_WGS_LEGACY_BASE | 1000G-development | wgs | SimpleTwoThirdsConsensus | invalid | invalid | 79 / 411 | 0.192214 | no |
| DEV_WGS_LEGACY_GCC | 1000G-development | wgs | TwoThirdsGuardedCC | invalid | invalid | 205 / 411 | 0.498783 | no |
| META_WGS_INTERNAL | 1000G-development | wgs | MetaConsensus | invalid | invalid | 220 / 411 | 0.535280 | no |
| META_NCI60_WES | NCI-60 | wes | MetaConsensus | exploratory | valid | 67 / 82 | 0.817073 | no |
| META_NCI60_RNA | NCI-60 | rnaseq | MetaConsensus | exploratory | valid | 65 / 82 | 0.792683 | no |
| REFFORMER_GATE | 1000G-development | all | RefFormer | exploratory | valid | — | — | no |
| EGCC_GATE | 1000G-development | wes+rnaseq | EvidenceGatedCC | exploratory | valid | — | — | no |

## Generated result sentences

- SimplePluralityLex produced 365 correct genotypes among 390 eligible wes loci in 1000G-development (development) [RESULT:DEV_WES_CONSENSUS].
- SimplePluralityLex produced 306 correct genotypes among 321 eligible rnaseq loci in 1000G-development (development) [RESULT:DEV_RNA_CONSENSUS].
- SimpleTwoThirdsConsensus produced 328 correct genotypes among 390 eligible wes loci in 1000G-development (development) [RESULT:DEV_TT_WES_BASE].
- TwoThirdsGuardedCC produced 369 correct genotypes among 390 eligible wes loci in 1000G-development (development) [RESULT:DEV_TT_WES_GCC].
- SimpleTwoThirdsConsensus produced 296 correct genotypes among 321 eligible rnaseq loci in 1000G-development (development) [RESULT:DEV_TT_RNA_BASE].
- TwoThirdsGuardedCC produced 306 correct genotypes among 321 eligible rnaseq loci in 1000G-development (development) [RESULT:DEV_TT_RNA_GCC].
- MetaConsensus produced 67 correct genotypes among 82 eligible wes loci in NCI-60 (exploratory) [RESULT:META_NCI60_WES].
- MetaConsensus produced 65 correct genotypes among 82 eligible rnaseq loci in NCI-60 (exploratory) [RESULT:META_NCI60_RNA].

## Invalid diagnostic records (supplement only)

- INVALID DIAGNOSTIC: SimplePluralityLex recorded 164 of 411 loci; this cannot be used as performance evidence [RESULT:DEV_WGS_CONSENSUS].
- INVALID DIAGNOSTIC: SimpleTwoThirdsConsensus recorded 79 of 411 loci; this cannot be used as performance evidence [RESULT:DEV_WGS_LEGACY_BASE].
- INVALID DIAGNOSTIC: TwoThirdsGuardedCC recorded 205 of 411 loci; this cannot be used as performance evidence [RESULT:DEV_WGS_LEGACY_GCC].
- INVALID DIAGNOSTIC: MetaConsensus recorded 220 of 411 loci; this cannot be used as performance evidence [RESULT:META_WGS_INTERNAL].
