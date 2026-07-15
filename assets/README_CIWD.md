# CIWD 3.0.0 catalogue — vendored lookup tables

**What this is.** A normalised, machine-readable copy of the **Common, Intermediate and
Well-Documented (CIWD) HLA alleles, version 3.0.0** catalogue, used by the champHLA benchmark as an
*allele-commonness annotation layer* — **not** as a ground-truth cohort. CIWD contains no per-sample
genotypes and no sequencing reads (it is an allele-frequency classification derived from ~8M
donor-registry typings), so it cannot serve as benchmark truth; see
`../docs/NEW_HLA_BENCHMARK_DATASETS.md`.

**Citation.** Hurley CK, Kempenich J, Wadsworth K, et al. *Common, intermediate and well-documented
HLA alleles in world populations: CIWD version 3.0.0.* HLA. 2020;95(6):516-531.
doi:10.1111/tan.13811 (PMID 31970929; PMC7317522). Reference IPD-IMGT/HLA build **3.31.0**.

**Categories** (`ciwd_total` and the seven per-population columns):
`C` = common (≥1/10⁴), `I` = intermediate (≥1/10⁵), `WD` = well-documented (≥5 occurrences),
`o` = not-CIWD in 3.0.0 (listed because it was CWD in an earlier catalogue), blank = not assigned in
that population group. `bin/ciwd.py` maps these to `common / intermediate / well_documented /
not_ciwd`, and reports alleles absent from the table as `unknown` (novel / never observed).

**Population groups:** AFA (African), API (Asian/Pacific Islander), EURO (European),
MENA (Middle-East/North-African), HIS (Hispanic), NAM (Native American), UNK (unknown).

## Files

| File | Source xlsx | Rows |
|------|-------------|------|
| `ciwd_3.0.0.tsv` | `CIWD_Palleles-all_loci-2020320-ihws-website.xlsx` (P group, two-field) | 3250 alleles |
| `ciwd_3.0.0_ggroup.tsv` | `SupTbl5Rev-Ggroup-CIWD-20200323_ihws-website.xlsx` (G group) | 3758 alleles |

Columns: `allele, designation, ciwd_highest, ciwd_total, ciwd_AFA, ciwd_API, ciwd_EURO, ciwd_MENA,
ciwd_HIS, ciwd_NAM, ciwd_UNK, cwd2_category` (`cwd2_category` = the older ASHI CWD 2.0.0 designation,
for comparability with labs/tools keyed on the 2012 catalogue).

## Reproduce

Download date: 2026-07-09. Source: 18th IHIW website (AWS S3 mirror).

```bash
mkdir -p assets/ciwd_raw && cd assets/ciwd_raw
P=https://s3.eu-central-1.amazonaws.com/ihiw.website.data/CIWD-3.0
curl -sSfLO $P/CIWD_Palleles-all_loci-2020320-ihws-website.xlsx     # P group, two-field
curl -sSfLO $P/SupTbl5Rev-Ggroup-CIWD-20200323_ihws-website.xlsx    # G group
cd ../.. && python3 assets/build_ciwd_tsv.py                        # needs openpyxl (python-data/3.12)
```

The raw `.xlsx` files live in `assets/ciwd_raw/` and are the unmodified published supplements; the
`.tsv` files are the only artefacts the pipeline reads.
