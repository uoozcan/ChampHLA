# Compact plurality head-to-head

The current comparator manifest contains 7 integration or historical-ablation methods.

| Modality | Plurality correct / loci | Plurality call rate | Comparator role | Comparator | Comparator correct | Comparator − plurality | Simultaneous 95% CI | Holm-adjusted p | Status |
|---|---:|---:|---|---|---:|---:|---:|---:|---|
| wgs | — | — | — | — | — | — | — | — | invalid source excluded; awaiting corrected rerun |
| wes | 365 / 390 | 100.0% | best observed caller (descriptive) | Caller:OptiType | 360 | -1.28 | [-3.85, 1.03] | 1 | valid development result |
| wes | 365 / 390 | 100.0% | strongest integration | TwoThirdsGuardedCC | 365 | 0.00 | [-1.54, 1.28] | 1 | valid development result |
| rnaseq | 306 / 321 | 100.0% | best observed caller (descriptive) | Caller:HLA-HD | 297 | -2.80 | [-7.79, 0.00] | 0.625 | valid development result |
| rnaseq | 306 / 321 | 100.0% | strongest integration | MVFlooredCC | 306 | 0.00 | [0.00, 0.00] | 1 | valid development result |
