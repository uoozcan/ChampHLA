#!/usr/bin/env python3.11
"""
Create pseudo-IDs (AML_1, AML_2, ...) for all FIMM sample IDs,
write the mapping table, and produce de-identified copies of every
input file.

Originals are never modified.
"""

from pathlib import Path
import pandas as pd

# ---------------------------------------------------------------------------
# Input files that contain a "sample_id" column with real FIMM IDs
# ---------------------------------------------------------------------------
INPUT_FILES = [
    Path("/scratch/project_2008084/pihla_local/fimm_results/fimm_hla_calls.tsv"),
    Path("/scratch/project_2008084/pihla_local/fimm_results/fimm_spechla_freq_rna.tsv"),
    Path("/scratch/project_2008084/pihla_local/fimm_results/fimm_spechla_freq_wes.tsv"),
    Path("/scratch/project_2008084/pihla_local/analysis/loh_analysis/fimm_loh_candidates_wes.tsv"),
    Path("/scratch/project_2008084/pihla-publish/fimm_results/fimm_hla_calls.tsv"),
    Path("/scratch/project_2008084/pihla-publish/fimm_results/fimm_scrna_matched.tsv"),
]

MAP_OUT_PATHS = [
    Path("/scratch/project_2008084/pihla_local/fimm_results/fimm_pseudo_id_map.tsv"),
    Path("/scratch/project_2008084/pihla-publish/fimm_results/fimm_pseudo_id_map.tsv"),
]

# ---------------------------------------------------------------------------
# 1. Collect all unique sample IDs across every input file
# ---------------------------------------------------------------------------
all_ids: set[str] = set()

for fpath in INPUT_FILES:
    if not fpath.exists():
        print(f"  SKIP (not found): {fpath}")
        continue
    df = pd.read_csv(fpath, sep="\t", dtype=str)
    if "sample_id" not in df.columns:
        print(f"  WARN: no 'sample_id' column in {fpath.name}")
        continue
    ids = df["sample_id"].dropna().unique().tolist()
    all_ids.update(ids)
    print(f"  found {len(ids):4d} unique IDs in {fpath.name}")

all_ids_sorted = sorted(all_ids)
print(f"\nTotal unique sample IDs across all files: {len(all_ids_sorted)}")

# ---------------------------------------------------------------------------
# 2. Build mapping:  real_sample_id -> AML_N  (1-indexed, alphabetical order)
# ---------------------------------------------------------------------------
mapping: dict[str, str] = {
    sid: f"AML_{i+1}"
    for i, sid in enumerate(all_ids_sorted)
}

# ---------------------------------------------------------------------------
# 3. Write mapping table
# ---------------------------------------------------------------------------
map_df = pd.DataFrame(
    [{"real_sample_id": k, "pseudo_id": v} for k, v in mapping.items()]
)

for out_path in MAP_OUT_PATHS:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    map_df.to_csv(out_path, sep="\t", index=False)
    print(f"  Wrote mapping -> {out_path}  ({len(map_df)} rows)")

# ---------------------------------------------------------------------------
# 4. Write de-identified copies  (original files untouched)
# ---------------------------------------------------------------------------
print("\nDe-identifying files:")
deident_paths = []

for fpath in INPUT_FILES:
    if not fpath.exists():
        continue
    df = pd.read_csv(fpath, sep="\t", dtype=str)
    if "sample_id" not in df.columns:
        continue

    unmapped = set(df["sample_id"].dropna().unique()) - set(mapping.keys())
    if unmapped:
        print(f"  WARN {fpath.name}: {len(unmapped)} IDs not in map — keeping as-is")

    df["sample_id"] = df["sample_id"].map(mapping).fillna(df["sample_id"])

    out_path = fpath.with_name(fpath.stem + "_deident.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    deident_paths.append(out_path)
    print(f"  {fpath.name:45s} -> {out_path.name}")

# ---------------------------------------------------------------------------
# 5. Verification: spot-check no real IDs remain in de-identified files
# ---------------------------------------------------------------------------
print("\nVerification (grep for known real ID prefixes):")
real_prefixes = ["FH_", "FHRB_", "FHRB-", "FHRB.", "VX_", "vcp-"]
all_clean = True
for dp in deident_paths:
    text = dp.read_text()
    hits = [p for p in real_prefixes if p in text]
    if hits:
        print(f"  FAIL {dp.name}: found real-ID-like strings {hits}")
        all_clean = False
    else:
        print(f"  OK   {dp.name}")

if all_clean:
    print("\nAll de-identified files: CLEAN (no real ID prefixes found).")
else:
    print("\nWARNING: Some files may still contain real IDs — review above.")

print(f"\nDone.  {len(all_ids_sorted)} IDs mapped.  {len(deident_paths)} files de-identified.")
