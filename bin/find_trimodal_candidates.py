#!/usr/bin/env python3
"""
find_trimodal_candidates.py — Enumerate all 1000G samples eligible for
trimodal analysis (WGS + WES + RNA-seq) and report processing gaps.

Outputs:
  - Summary table to stdout
  - sample_fastq_urls_batch5.tsv  (RNA index for Geuvadis+truth samples needing RNA)
  - wes_priority_manifest.txt     (sample list needing WES, ready for slurm_wes_batch10.sh)

Usage:
  python3 bin/find_trimodal_candidates.py [--dry-run]
"""
import argparse
import os
import sys
from collections import defaultdict

# ── Paths ─────────────────────────────────────────────────────────────────────
CALIB_ROOT    = "/scratch/project_2008084/hla_calibration"
PIHLA_ROOT    = "/scratch/project_2008084/pihla-publish"
ANALYSIS_ROOT = f"{PIHLA_ROOT}/analysis"

GEUVADIS_TABLE = f"{CALIB_ROOT}/rna/index/geuvadis_run_table.tsv"
TRUTH_FILE     = f"{PIHLA_ROOT}/analysis/1000g_realdata/full_cohort_inputs/truth_long.tsv"

WES_INDICES = [
    f"{CALIB_ROOT}/wes/index/sample_bam_urls_batch2.tsv",
    f"{CALIB_ROOT}/wes/index/sample_bam_urls_batch3.tsv",
    f"{CALIB_ROOT}/wes/index/sample_bam_urls_batch4.tsv",
    f"{CALIB_ROOT}/wes/index/sample_bam_urls_recovery_wave001.tsv",
    f"{CALIB_ROOT}/wes/index/sample_bam_urls.tsv",
]
WGS_CRAM_INDEX = f"{CALIB_ROOT}/wgs/index/sample_cram_urls.tsv"

WGS_BATCHES = f"{CALIB_ROOT}/wgs_batches"
WES_BATCHES = f"{CALIB_ROOT}/wes_batches"
RNA_BATCHES = f"{CALIB_ROOT}/rna_batches"

OUT_RNA_BATCH5  = f"{CALIB_ROOT}/rna/index/sample_fastq_urls_batch5.tsv"
OUT_WES_MANIFEST = f"{ANALYSIS_ROOT}/wes_priority_manifest.txt"


def load_geuvadis(path):
    """
    Returns:
      sample_runs: dict[sample_id -> list of (run_acc, url_r1, url_r2)]
    """
    sample_runs = defaultdict(list)
    with open(path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            run_acc     = parts[0]
            sample_alias = parts[1]          # e.g. "GEUV:HG00096"
            fastq_ftp   = parts[3]           # "url_r1;url_r2"
            sample_id = sample_alias.replace("GEUV:", "").strip()
            urls = fastq_ftp.split(";")
            if len(urls) >= 2:
                sample_runs[sample_id].append((run_acc, urls[0].strip(), urls[1].strip()))
    return sample_runs


def load_truth_samples(path):
    """Returns set of sample IDs present in the truth file."""
    samples = set()
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.split("\t")
            if parts:
                samples.add(parts[0].strip())
    return samples


def load_wes_index(paths):
    """Returns dict[sample_id -> primary_bam_url] from all WES index files (first match wins)."""
    index = {}
    for path in paths:
        if not os.path.isfile(path):
            continue
        with open(path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    sample = parts[0].strip()
                    url    = parts[1].strip()
                    if sample and url and sample not in index:
                        index[sample] = url
    return index


def load_wgs_cram_index(path):
    """Returns set of sample IDs in the EBI CRAM index."""
    samples = set()
    if not os.path.isfile(path):
        return samples
    with open(path) as f:
        for line in f:
            parts = line.split("\t")
            if parts:
                samples.add(parts[0].strip())
    return samples


def processed_samples(batch_dir):
    """Returns set of sample IDs that have a results directory under batch_dir."""
    if not os.path.isdir(batch_dir):
        return set()
    return {
        d for d in os.listdir(batch_dir)
        if os.path.isdir(os.path.join(batch_dir, d))
        and not d.startswith(".")
        and "SMOKE" not in d
    }


def wgs_processed_samples(batch_dir):
    """
    WGS batches contain per-sample directories. Returns those that have at
    least one tool result (optitype result file is a reliable marker).
    """
    done = set()
    if not os.path.isdir(batch_dir):
        return done
    for sample in os.listdir(batch_dir):
        sample_dir = os.path.join(batch_dir, sample)
        if not os.path.isdir(sample_dir):
            continue
        results = os.path.join(sample_dir, "results", sample, "optitype")
        if os.path.isdir(results) and any(
            f.endswith("_result.tsv") for f in os.listdir(results)
        ):
            done.add(sample)
        # Older per-tool layout: results/NA.../optitype/...
        elif any(
            "optitype" in root
            for root, dirs, files in os.walk(os.path.join(sample_dir, "results"))
        ):
            done.add(sample)
    return done


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be written without creating files")
    args = parser.parse_args()

    # ── Load data ──────────────────────────────────────────────────────────────
    print("Loading Geuvadis run table...", flush=True)
    geuvadis   = load_geuvadis(GEUVADIS_TABLE)          # 464 samples
    print(f"  {len(geuvadis)} unique Geuvadis samples")

    print("Loading HLA truth...", flush=True)
    truth      = load_truth_samples(TRUTH_FILE)         # 138 samples
    print(f"  {len(truth)} truth samples")

    print("Loading WES index...", flush=True)
    wes_index  = load_wes_index(WES_INDICES)
    print(f"  {len(wes_index)} samples in WES index")

    print("Loading WGS EBI CRAM index...", flush=True)
    wgs_cram_indexed = load_wgs_cram_index(WGS_CRAM_INDEX)
    print(f"  {len(wgs_cram_indexed)} samples in EBI CRAM index")

    print("Scanning processed batches...", flush=True)
    wgs_done = wgs_processed_samples(WGS_BATCHES)
    wes_done = processed_samples(WES_BATCHES)
    rna_done = processed_samples(RNA_BATCHES)
    print(f"  WGS done: {len(wgs_done)}")
    print(f"  WES done: {len(wes_done)}")
    print(f"  RNA done: {len(rna_done)}")

    # ── Candidate universe: Geuvadis ∩ truth ──────────────────────────────────
    universe = set(geuvadis.keys()) & truth
    print(f"\nGeuvadis+truth intersection (candidate universe): {len(universe)} samples")

    # ── Classify each candidate ───────────────────────────────────────────────
    trimodal        = []   # all 3 done
    needs_wes_only  = []   # WGS+RNA done, WES missing
    needs_rna_only  = []   # WGS+WES done, RNA missing
    needs_wes_rna   = []   # WGS done, WES+RNA missing
    needs_wgs       = []   # WES+RNA done, WGS missing (NYGC-blocked, may have EBI CRAM)
    needs_all_three = []   # none of the 3 done

    for s in sorted(universe):
        has_wgs = s in wgs_done
        has_wes = s in wes_done
        has_rna = s in rna_done

        if has_wgs and has_wes and has_rna:
            trimodal.append(s)
        elif has_wgs and has_rna and not has_wes:
            needs_wes_only.append(s)
        elif has_wgs and has_wes and not has_rna:
            needs_rna_only.append(s)
        elif has_wgs and not has_wes and not has_rna:
            needs_wes_rna.append(s)
        elif has_wes and has_rna and not has_wgs:
            needs_wgs.append(s)
        else:
            needs_all_three.append(s)

    total_actionable = (len(needs_wes_only) + len(needs_rna_only)
                        + len(needs_wes_rna) + len(needs_wgs))

    # ── Print summary ─────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("  TRIMODAL CANDIDATE SUMMARY")
    print("=" * 72)
    print(f"  Candidate universe (Geuvadis+truth): {len(universe):>4}")
    print(f"  Already trimodal (all 3 done):       {len(trimodal):>4}  ← current")
    print(f"  Need WES only (WGS+RNA done):        {len(needs_wes_only):>4}")
    print(f"  Need RNA only (WGS+WES done):        {len(needs_rna_only):>4}")
    print(f"  Need WES+RNA  (WGS done):            {len(needs_wes_rna):>4}")
    print(f"  WGS-blocked   (WES+RNA done):        {len(needs_wgs):>4}  ← EBI CRAM check below")
    print(f"  Need all three:                      {len(needs_all_three):>4}")
    print(f"  Total actionable (non-WGS-blocked):  {total_actionable:>4}")
    trimodal_after = len(trimodal) + len(needs_wes_only) + len(needs_rna_only) + len(needs_wes_rna)
    print(f"  Trimodal after completing above:     {trimodal_after:>4}")
    print("=" * 72)

    def show_list(label, samples):
        if samples:
            print(f"\n  {label} ({len(samples)}):")
            for s in samples:
                print(f"    {s}")

    show_list("Need WES only",  needs_wes_only)
    show_list("Need RNA only",  needs_rna_only)
    show_list("Need WES+RNA",   needs_wes_rna)
    show_list("WGS-blocked — WES+RNA done (check EBI CRAM)", needs_wgs)
    if needs_all_three:
        show_list("Need all three (lowest priority)", needs_all_three)

    # ── EBI CRAM check for WGS-blocked ────────────────────────────────────────
    wgs_recoverable = [s for s in needs_wgs if s in wgs_cram_indexed]
    if wgs_recoverable:
        print(f"\n  WGS-blocked samples IN EBI CRAM index ({len(wgs_recoverable)}):")
        for s in wgs_recoverable:
            print(f"    {s}  → can recover via slurm_download_wgs_cram.sh")

    # ── Truth-only samples (not in Geuvadis) ─────────────────────────────────
    truth_no_rna = truth - set(geuvadis.keys())
    print(f"\n  Truth samples with NO Geuvadis RNA data: {len(truth_no_rna)}")
    print("  (These cannot reach trimodal via Geuvadis — no action needed)")
    print()

    # ── Emit batch5 RNA index ─────────────────────────────────────────────────
    rna_needed = sorted(needs_rna_only + needs_wes_rna)
    if rna_needed:
        print(f"RNA batch5: {len(rna_needed)} samples need RNA index entries")
        rows = []
        missing_in_geuvadis = []
        for s in rna_needed:
            runs = geuvadis.get(s, [])
            if not runs:
                missing_in_geuvadis.append(s)
                continue
            # Use first run; URL format matches existing index (no ftp:// prefix)
            _, url_r1, url_r2 = runs[0]
            rows.append((s, url_r1, url_r2))

        if missing_in_geuvadis:
            print(f"  WARNING: not found in Geuvadis table: {missing_in_geuvadis}")

        if rows:
            if args.dry_run:
                print(f"  [dry-run] would write {len(rows)} rows to {OUT_RNA_BATCH5}")
                for s, r1, r2 in rows:
                    print(f"    {s}\t{r1}\t{r2}")
            else:
                os.makedirs(os.path.dirname(OUT_RNA_BATCH5), exist_ok=True)
                with open(OUT_RNA_BATCH5, "w") as f:
                    for s, r1, r2 in rows:
                        f.write(f"{s}\t{r1}\t{r2}\n")
                print(f"  Written: {OUT_RNA_BATCH5}")
                for s, r1, r2 in rows:
                    print(f"    {s}\t{r1}\t{r2}")
    else:
        print("No samples need a new RNA index entry.")

    # ── Emit WES priority manifest ─────────────────────────────────────────────
    wes_needed = sorted(needs_wes_only + needs_wes_rna)
    # Check which are actually in the WES index
    wes_indexable  = [s for s in wes_needed if s in wes_index]
    wes_not_indexed = [s for s in wes_needed if s not in wes_index]

    if wes_needed:
        print(f"\nWES priority manifest: {len(wes_needed)} samples need WES")
        print(f"  In WES index (actionable): {len(wes_indexable)}")
        if wes_not_indexed:
            print(f"  NOT in any WES index (blocked): {wes_not_indexed}")

        if args.dry_run:
            print(f"  [dry-run] would write {len(wes_indexable)} lines to {OUT_WES_MANIFEST}")
        else:
            os.makedirs(os.path.dirname(OUT_WES_MANIFEST), exist_ok=True)
            with open(OUT_WES_MANIFEST, "w") as f:
                for s in wes_indexable:
                    f.write(f"{s}\n")
            print(f"  Written: {OUT_WES_MANIFEST}")
        for s in wes_indexable:
            print(f"    {s}")

    print()
    print("Next steps:")
    if rna_needed:
        print(f"  1. bash {PIHLA_ROOT}/submit_rna_6trimodal.sh")
    if wes_indexable:
        print(f"  2. sbatch {PIHLA_ROOT}/slurm_wes_batch10.sh \\")
        print(f"       {' '.join(wes_indexable)}")
    print(f"  3. sbatch {PIHLA_ROOT}/slurm_benchmark_trimodal.sh  (after jobs complete)")


if __name__ == "__main__":
    main()
