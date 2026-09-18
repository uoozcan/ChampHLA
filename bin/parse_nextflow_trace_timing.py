"""
parse_nextflow_trace_timing.py — Aggregate per-tool timing from Nextflow execution traces.

Reads execution_trace.txt files produced by the ChampHLA Nextflow pipeline (trace
enabled in nextflow.config) and computes per-tool × per-modality timing statistics
(median wall-clock time, peak RAM, CPU utilisation) across n≥30 samples.

Usage:
    python3.11 bin/parse_nextflow_trace_timing.py \
        --wes-traces  /scratch/.../hla_calibration/wes_batches/*/results/pipeline_info/execution_trace.txt \
        --rna-traces  /scratch/.../hla_calibration/rna_batches/*/results/pipeline_info/execution_trace.txt \
        --wgs-traces  analysis/performance_benchmark_n30/runs/wgs/*/pipeline_info/execution_trace.txt \
        --n-per-modality 30 \
        --out-dir analysis/performance_benchmark_n30/tables/

Outputs:
    timing_summary_n30.tsv   — one row per tool × modality with median/P25/P75

Trace fields used (columns from execution_trace.txt):
    name      — process name, e.g. OPTITYPE (NA12345)
    status    — COMPLETED | FAILED | ABORTED
    realtime  — wall-clock duration string, e.g. "1m 30s", "2h 15m 3s"
    %cpu      — average CPU utilisation, e.g. "143.8%"
    peak_rss  — peak resident set size, e.g. "1.4 GB" or "340.8 MB"
"""

import argparse
import csv
import glob
import re
import statistics
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Process name → tool name mapping
# ---------------------------------------------------------------------------
PROCESS_TOOL_MAP = {
    "OPTITYPE":     "OptiType",
    "ARCASHLA":     "ArcasHLA",
    "HLAHD":        "HLA-HD",
    "T1K":          "T1K",
    "SPECHLA":      "SpecHLA",
    "KOURAMI":      "Kourami",
    "POLYSOLVER":   "POLYSOLVER",
    "SEQ2HLA":      "Seq2HLA",
}

# Processes to skip (preprocessing steps, not tool typing)
SKIP_PREFIXES = {
    "BAM_TO_FASTQ", "EXTRACT_HLA", "CRAM_TO_BAM", "AGGREGATE",
    "MAJORITY_VOTING", "HARMONISE", "CONSENSUS",
}


def parse_duration_to_hours(s: str) -> float | None:
    """
    Parse Nextflow human-readable duration string to decimal hours.
    Examples:
        "4.7s"       → 0.001306
        "1m 30s"     → 0.025
        "2h 15m 3s"  → 2.25083
        "1h 2m 9s"   → 1.0358
    """
    if not s or s.strip() in ("-", ""):
        return None
    s = s.strip()
    hours = 0.0
    # hours
    m = re.search(r"(\d+(?:\.\d+)?)\s*h", s)
    if m:
        hours += float(m.group(1))
    # minutes
    m = re.search(r"(\d+(?:\.\d+)?)\s*m(?!s)", s)
    if m:
        hours += float(m.group(1)) / 60.0
    # seconds
    m = re.search(r"(\d+(?:\.\d+)?)\s*(?:ms|s)", s)
    if m:
        unit_str = s[m.start():]
        if "ms" in unit_str[:4]:
            hours += float(m.group(1)) / 3_600_000.0
        else:
            hours += float(m.group(1)) / 3_600.0
    return hours if hours > 0 else None


def parse_memory_to_gb(s: str) -> float | None:
    """
    Parse Nextflow memory string to GB.
    Examples: "340.8 MB" → 0.333, "10.7 GB" → 10.7, "1.4 TB" → 1433.6
    """
    if not s or s.strip() in ("-", ""):
        return None
    s = s.strip()
    m = re.match(r"([\d.]+)\s*(KB|MB|GB|TB)", s, re.IGNORECASE)
    if not m:
        return None
    val, unit = float(m.group(1)), m.group(2).upper()
    return {
        "KB": val / 1_048_576,
        "MB": val / 1_024,
        "GB": val,
        "TB": val * 1_024,
    }[unit]


def parse_cpu_pct(s: str) -> float | None:
    """Parse "%cpu" field like "143.8%" or "143.8" → 143.8."""
    if not s or s.strip() in ("-", ""):
        return None
    return float(s.strip().rstrip("%"))


def process_name_to_tool(name: str) -> str | None:
    """
    Map a Nextflow process name (e.g. 'OPTITYPE (NA12345)') to a tool label.
    Returns None if the process is a preprocessing step.
    """
    proc = name.split("(")[0].strip().upper()
    for skip in SKIP_PREFIXES:
        if proc.startswith(skip):
            return None
    for key, label in PROCESS_TOOL_MAP.items():
        if proc.startswith(key):
            return label
    return None


def read_trace(path: str) -> list[dict]:
    """Read a single execution_trace.txt and return list of row dicts."""
    rows = []
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rows.append(row)
    except Exception as e:
        print(f"  WARN: could not read {path}: {e}", file=sys.stderr)
    return rows


def extract_sample_from_name(name: str) -> str:
    """Extract sample ID from process name like 'OPTITYPE (NA12345)'."""
    m = re.search(r"\(([^)]+)\)", name)
    return m.group(1).strip() if m else "unknown"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_traces(trace_paths: list[str], modality: str, n_max: int) -> list[dict]:
    """
    Parse a list of trace files for one modality.
    Returns list of per-measurement dicts with keys:
        tool, modality, sample, realtime_h, peak_rss_gb, cpu_pct
    """
    records = []
    seen_samples: set[str] = set()

    # Map sample → list of tool rows (to check completeness)
    sample_tool_rows: dict[str, list[dict]] = {}

    for path in trace_paths:
        rows = read_trace(path)
        for row in rows:
            if row.get("status", "").strip() != "COMPLETED":
                continue
            tool = process_name_to_tool(row.get("name", ""))
            if tool is None:
                continue
            sample = extract_sample_from_name(row.get("name", ""))
            sample_tool_rows.setdefault(sample, []).append({
                "tool": tool,
                "realtime_h": parse_duration_to_hours(row.get("realtime", "")),
                "peak_rss_gb": parse_memory_to_gb(row.get("peak_rss", "")),
                "cpu_pct": parse_cpu_pct(row.get("%cpu", "")),
                "sample": sample,
                "modality": modality,
            })

    # Select samples with the most tools completed (prefer fully complete ones)
    sorted_samples = sorted(
        sample_tool_rows.items(),
        key=lambda x: len(x[1]),
        reverse=True,
    )

    selected = 0
    for sample, tool_rows in sorted_samples:
        if selected >= n_max:
            break
        if sample in seen_samples:
            continue
        seen_samples.add(sample)
        records.extend(tool_rows)
        selected += 1

    print(
        f"  {modality}: selected {selected} samples "
        f"(from {len(sample_tool_rows)} with data); "
        f"{len(records)} tool-sample records",
        file=sys.stderr,
    )
    return records


def aggregate(records: list[dict]) -> list[dict]:
    """
    Group by tool × modality and compute median / P25 / P75.
    Returns list of summary row dicts.
    """
    from collections import defaultdict
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in records:
        groups[(r["tool"], r["modality"])].append(r)

    summary = []
    for (tool, modality), rows in sorted(groups.items()):
        runtimes = [r["realtime_h"] for r in rows if r["realtime_h"] is not None]
        rams = [r["peak_rss_gb"] for r in rows if r["peak_rss_gb"] is not None]
        cpus = [r["cpu_pct"] for r in rows if r["cpu_pct"] is not None]

        def med(vals):
            return round(statistics.median(vals), 4) if vals else None

        def p25(vals):
            if not vals:
                return None
            s = sorted(vals)
            return round(s[int(len(s) * 0.25)], 4)

        def p75(vals):
            if not vals:
                return None
            s = sorted(vals)
            return round(s[int(len(s) * 0.75)], 4)

        summary.append({
            "tool": tool,
            "modality": modality,
            "n_samples": len(rows),
            "median_runtime_hours": med(runtimes),
            "p25_runtime_hours": p25(runtimes),
            "p75_runtime_hours": p75(runtimes),
            "median_peak_ram_gb": med(rams),
            "p25_peak_ram_gb": p25(rams),
            "p75_peak_ram_gb": p75(rams),
            "median_cpu_pct": med(cpus),
            "p25_cpu_pct": p25(cpus),
            "p75_cpu_pct": p75(cpus),
        })
    return summary


def write_tsv(rows: list[dict], path: Path) -> None:
    if not rows:
        print(f"WARNING: no data to write to {path}", file=sys.stderr)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wes-traces", nargs="*", default=[], metavar="GLOB")
    parser.add_argument("--rna-traces", nargs="*", default=[], metavar="GLOB")
    parser.add_argument("--wgs-traces", nargs="*", default=[], metavar="GLOB")
    parser.add_argument("--n-per-modality", type=int, default=30)
    parser.add_argument(
        "--out-dir",
        default="analysis/performance_benchmark_n30/tables",
        metavar="DIR",
    )
    args = parser.parse_args()

    all_records = []
    n = args.n_per_modality

    for modality, patterns in [
        ("WES",  args.wes_traces),
        ("RNA",  args.rna_traces),
        ("WGS",  args.wgs_traces),
    ]:
        # Expand globs (shell may not expand wildcards passed as string args)
        paths = []
        for pat in patterns:
            expanded = glob.glob(pat)
            paths.extend(expanded if expanded else [pat])

        if not paths:
            print(f"  {modality}: no trace files provided, skipping", file=sys.stderr)
            continue

        print(f"\nProcessing {modality} ({len(paths)} trace files)...", file=sys.stderr)
        records = parse_traces(paths, modality, n)
        all_records.extend(records)

    if not all_records:
        print("ERROR: no records parsed. Check trace file paths.", file=sys.stderr)
        sys.exit(1)

    summary = aggregate(all_records)
    out = Path(args.out_dir) / "timing_summary_n30.tsv"
    write_tsv(summary, out)
    print(f"\nDone. {len(summary)} tool×modality combinations in {out}")


if __name__ == "__main__":
    main()
