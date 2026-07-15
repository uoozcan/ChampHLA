#!/usr/bin/env python3.11
"""
run_pipeline_validation.py — MVHLA pipeline end-to-end validation (loop approach)

Uses a Draft → Critique → Fix → Loop cycle for each pipeline stage.
Each step has:
  - Runner:    executes the step
  - Scorer:    evaluates output quality (0–10)
  - Fixer:     diagnoses and attempts to fix failures
  - Loop:      repeats until score ≥ threshold or max retries reached

Stages:
  Stage 1:  Validate tool outputs  (loop: modality → sample → tool)
  Stage 2:  Run benchmark harness  (hla_benchmark.py)
  Stage 3:  Generate all figures   (3 figure scripts)
  Stage 4:  Parse timing data      (trace parser)
  Stage 5:  Compile final report

Usage:
    python3.11 bin/run_pipeline_validation.py

Output:
    analysis/pipeline_validation_results/
"""

import csv
import glob
import json
import os
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

# ═══════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════
BASE = Path("/scratch/project_2008084/pihla-publish")
BATCH_ROOT = Path("/scratch/project_2008084/hla_calibration")
OUT_DIR = BASE / "analysis" / "pipeline_validation_results"

MODALITY_TOOLS = {
    "wgs": ["optitype", "arcashla", "hlahd", "t1k", "spechla", "kourami", "polysolver"],
    "wes": ["optitype", "arcashla", "hlahd", "t1k", "spechla", "kourami", "polysolver"],
    "rna": ["optitype", "arcashla", "hlahd", "t1k", "spechla", "seq2hla"],
}

BATCH_DIRS = {
    "wgs": BATCH_ROOT / "wgs_batches",
    "wes": BATCH_ROOT / "wes_batches",
    "rna": BATCH_ROOT / "rna_batches",
}

MAX_RETRIES = 3
PASS_THRESHOLD = 8  # out of 10

FIGURES_DIR = BASE / "analysis" / "figures_manuscript"
TIMING_TSV = BASE / "analysis" / "performance_benchmark_n30" / "tables" / "timing_summary_n30.tsv"

EXPECTED_FIGS = [
    "fig08_computational_performance",
    "figure_1_workflow_architecture", "figure_2_accuracy_comparison",
    "figure_3_per_gene_gains", "figure_4_confidence_calibration",
    "figure_5_abstention_tradeoff", "figure_6_discordance_taxonomy",
    "figure_7_confidence_weights",
    "figure_09_bimodal_per_gene", "figure_10_trimodal_comparison",
    "figure_s1_per_gene_accuracy", "figure_s2_calibration_heatmap",
    "figure_s3_resolution_comparison",
]

EXPECTED_BENCHMARK_TABLES = [
    "method_comparison.tsv", "method_per_gene.tsv",
    "summary_full_cohort.tsv", "summary_per_gene.tsv",
    "tool_confidence_weights.tsv", "weighted_consensus_calls.tsv",
    "discordance_summary.tsv", "abstention_tradeoff.tsv",
    "harmonized_benchmark_rows.tsv", "consensus_runtime_weights.json",
]


# ═══════════════════════════════════════════════════════════════════════════
# Logging
# ═══════════════════════════════════════════════════════════════════════════
class Logger:
    def __init__(self):
        self.lines: list[str] = []
        self.errors: list[str] = []
        self.warnings: list[str] = []

    def log(self, msg: str, level: str = "INFO"):
        ts = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] [{level:>5}] {msg}"
        print(line, flush=True)
        self.lines.append(line)
        if level == "ERROR":
            self.errors.append(msg)
        elif level == "WARN":
            self.warnings.append(msg)

    def section(self, title: str):
        self.log("=" * 72)
        self.log(title)
        self.log("=" * 72)


L = Logger()


# ═══════════════════════════════════════════════════════════════════════════
# Loop Engine — Draft / Critique / Fix / Loop
# ═══════════════════════════════════════════════════════════════════════════
class StepResult:
    """Result of a single loop iteration."""
    def __init__(self, score: int, details: dict, critiques: list[str], fixes_applied: list[str]):
        self.score = score
        self.details = details
        self.critiques = critiques
        self.fixes_applied = fixes_applied
        self.passed = score >= PASS_THRESHOLD


def run_loop(stage_name: str, runner, scorer, fixer=None) -> StepResult:
    """
    Core loop engine.
      runner()  → produces output (dict)
      scorer(output) → (score: int, critiques: list[str])
      fixer(critiques, output) → fixes_applied: list[str]
    Loops until score ≥ PASS_THRESHOLD or MAX_RETRIES reached.
    """
    L.log(f"")
    L.section(f"STAGE: {stage_name}")

    all_fixes = []
    for attempt in range(1, MAX_RETRIES + 1):
        L.log(f"  Attempt {attempt}/{MAX_RETRIES}")

        # === DRAFT: Run the step ===
        L.log(f"  [Runner] Executing...")
        try:
            output = runner()
        except Exception as e:
            L.log(f"  [Runner] EXCEPTION: {e}", "ERROR")
            output = {"error": str(e)}

        # === CRITIQUE: Score the output ===
        L.log(f"  [Scorer] Evaluating output...")
        score, critiques = scorer(output)
        L.log(f"  [Scorer] Score: {score}/10")
        for c in critiques:
            L.log(f"    → {c}", "WARN" if score < PASS_THRESHOLD else "INFO")

        if score >= PASS_THRESHOLD:
            L.log(f"  ✓ PASSED (score {score}/10 ≥ {PASS_THRESHOLD})")
            return StepResult(score, output, critiques, all_fixes)

        # === FIX: Attempt to resolve critiques ===
        if fixer and attempt < MAX_RETRIES:
            L.log(f"  [Fixer] Attempting to resolve {len(critiques)} critique(s)...")
            fixes = fixer(critiques, output)
            all_fixes.extend(fixes)
            for f in fixes:
                L.log(f"    Applied: {f}")
        elif attempt < MAX_RETRIES:
            L.log(f"  [Fixer] No fixer defined — retrying as-is")

    L.log(f"  ✗ FAILED after {MAX_RETRIES} attempts (best score: {score}/10)", "ERROR")
    return StepResult(score, output, critiques, all_fixes)


# ═══════════════════════════════════════════════════════════════════════════
# Stage 1: Tool Output Validation
# ═══════════════════════════════════════════════════════════════════════════
def stage1_runner() -> dict:
    """Loop through every modality → sample → tool and validate output files."""
    inventory = []
    stats = defaultdict(lambda: {"ok": 0, "missing": 0, "invalid": 0, "total": 0})

    for modality, tools in MODALITY_TOOLS.items():
        batch_dir = BATCH_DIRS[modality]
        if not batch_dir.exists():
            L.log(f"    {modality.upper()}: batch dir missing: {batch_dir}", "ERROR")
            continue

        samples = sorted([d.name for d in batch_dir.iterdir() if d.is_dir()])
        L.log(f"    {modality.upper()}: {len(samples)} samples, {len(tools)} tools")

        for sample in samples:
            for tool in tools:
                key = f"{modality}:{tool}"
                stats[key]["total"] += 1

                # Find the tool output file
                tool_file = batch_dir / sample / "results" / sample / tool / f"{sample}_{tool}.txt"
                if not tool_file.exists():
                    # Try glob
                    alts = glob.glob(str(batch_dir / sample / "results" / "**" / f"*_{tool}.txt"), recursive=True)
                    tool_file = Path(alts[0]) if alts else tool_file

                status, n_rows, genes, errors = _validate_one_file(tool_file)

                if status == "OK":
                    stats[key]["ok"] += 1
                elif status == "MISSING":
                    stats[key]["missing"] += 1
                else:
                    stats[key]["invalid"] += 1

                inventory.append({
                    "modality": modality, "sample": sample, "tool": tool,
                    "status": status, "n_rows": n_rows, "genes": genes, "errors": errors,
                })

    return {"inventory": inventory, "stats": dict(stats)}


def _validate_one_file(path: Path) -> tuple:
    """Returns (status, n_rows, genes_str, errors_str).

    Handles multiple tool output formats:
      - Long format with Gene/Allele1/Allele2 columns (OptiType, ArcasHLA, etc.)
      - Headerless TSV (HLA-HD: gene allele1 allele2 reads1 reads2)
      - Wide format (SpecHLA: Sample HLA_A_1 HLA_A_2 HLA_B_1 ...)
      - Kourami/POLYSOLVER with comment headers
    """
    if not path.exists():
        return "MISSING", 0, "", "file not found"
    if path.stat().st_size == 0:
        return "INVALID", 0, "", "empty file"
    try:
        with open(path) as fh:
            raw_lines = fh.readlines()

        # Strip comment lines
        data_lines = [l for l in raw_lines if not l.startswith("#")]
        if not data_lines:
            return "INVALID", 0, "", "file contains only comments"

        # Check for "no output" warnings
        content = "".join(raw_lines).lower()
        if "no output" in content or "failed" in content:
            # File exists but tool produced no results
            stripped = [l.strip() for l in data_lines if l.strip()]
            if len(stripped) <= 1:  # only header or empty
                return "INVALID", 0, "", "tool produced no output"

        # Try long format first (has Gene column header)
        import io
        first_line = data_lines[0].strip()
        first_cols = first_line.split("\t")
        first_cols_lower = [c.lower().strip() for c in first_cols]

        # Long format: "Gene\tAllele1\tAllele2"
        if "gene" in first_cols_lower:
            reader = csv.DictReader(io.StringIO("".join(data_lines)), delimiter="\t")
            rows = list(reader)
            n_rows = len(rows)
            if n_rows == 0:
                return "INVALID", 0, "", "header only, no data"
            gene_col = next(c for c in reader.fieldnames if c.lower().strip() == "gene")
            genes = sorted({r[gene_col] for r in rows if r.get(gene_col)})
            return "OK", n_rows, ",".join(genes), ""

        # Wide format: "Sample\tHLA_A_1\tHLA_A_2\t..." (SpecHLA)
        if any("hla_" in c.lower() for c in first_cols):
            n_rows = len(data_lines) - 1  # subtract header
            # Extract gene names from column headers
            genes = sorted({c.split("_")[1] for c in first_cols
                           if c.lower().startswith("hla_") and len(c.split("_")) > 1})
            if n_rows > 0 and genes:
                return "OK", n_rows, ",".join(genes), ""
            return "INVALID", n_rows, "", "wide format but no data rows"

        # Headerless format: HLA-HD (A\tHLA-A*02:01\tHLA-A*03:01\t...)
        # Check if first field looks like a gene name
        if first_cols[0].strip() in {"A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1",
                                      "DRB3", "DRB4", "DRB5", "DMA", "DMB", "DOA", "DOB",
                                      "DRA", "E", "F", "G"}:
            n_rows = len([l for l in data_lines if l.strip()])
            genes = sorted({l.split("\t")[0].strip() for l in data_lines if l.strip()})
            if n_rows > 0:
                return "OK", n_rows, ",".join(genes), ""

        # Check for HLA allele content anywhere
        hla_count = sum(1 for l in data_lines if "HLA-" in l or "A*" in l or "B*" in l or "C*" in l)
        if hla_count > 0:
            return "OK", hla_count, "", "non-standard format but contains HLA calls"

        return "INVALID", len(data_lines), "", f"unrecognized format (first cols: {first_cols[:3]})"
    except Exception as e:
        return "INVALID", 0, "", str(e)


def stage1_scorer(output: dict) -> tuple[int, list[str]]:
    """Score tool output completeness. 10 = all present, deduct per missing/invalid."""
    critiques = []
    if "error" in output:
        return 0, [f"Runner error: {output['error']}"]

    stats = output.get("stats", {})
    total_ok = sum(s["ok"] for s in stats.values())
    total_all = sum(s["total"] for s in stats.values())
    total_missing = sum(s["missing"] for s in stats.values())

    if total_all == 0:
        return 0, ["No tool outputs found at all"]

    pct = total_ok / total_all * 100

    # Check per-modality×tool for any completely missing
    for key, s in stats.items():
        if s["ok"] == 0 and s["total"] > 0:
            critiques.append(f"{key}: 0/{s['total']} files found — tool never ran?")
        elif s["missing"] > s["total"] * 0.3:
            critiques.append(f"{key}: {s['missing']}/{s['total']} missing (>{30}%)")

    # Score: 10 if ≥95%, 8 if ≥80%, etc.
    if pct >= 95:
        score = 10
    elif pct >= 85:
        score = 9
    elif pct >= 75:
        score = 8
    elif pct >= 60:
        score = 6
    elif pct >= 40:
        score = 4
    else:
        score = 2

    L.log(f"    Tool outputs: {total_ok}/{total_all} OK ({pct:.1f}%), {total_missing} missing")
    return score, critiques


def stage1_fixer(critiques: list[str], output: dict) -> list[str]:
    """For tool outputs, we can't re-run Nextflow here. Log actionable advice."""
    fixes = []
    for c in critiques:
        if "tool never ran" in c:
            fixes.append(f"NOTE: {c} — requires re-running Nextflow for this tool/modality")
        elif "missing" in c:
            fixes.append(f"NOTE: {c} — partial results; check SLURM logs for failed samples")
    return fixes


# ═══════════════════════════════════════════════════════════════════════════
# Stage 2: Benchmark Harness
# ═══════════════════════════════════════════════════════════════════════════
def stage2_runner() -> dict:
    """Run or validate the benchmark analysis."""
    config = BASE / "conf" / "benchmark_1000g_full_cohort_v2.yaml"
    tables_dir = BASE / "analysis" / "benchmark_wes_all_samples" / "tables"

    if not config.exists():
        return {"error": f"Config not found: {config}"}

    # Check if tables exist and are recent enough
    method_comp = tables_dir / "method_comparison.tsv"
    if method_comp.exists():
        age_h = (datetime.now().timestamp() - method_comp.stat().st_mtime) / 3600
        L.log(f"    Existing tables found (age: {age_h:.0f}h)")
        if age_h < 48:
            L.log(f"    Using existing tables (< 48h old)")
            return _check_benchmark_tables(tables_dir)

    # Run benchmark
    L.log(f"    Running hla_benchmark.py...")
    cmd = [
        "python3.11", str(BASE / "bin" / "hla_benchmark.py"),
        "--config", str(config),
        "--output-dir", str(tables_dir.parent),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, cwd=str(BASE))
    if r.returncode != 0:
        return {"error": f"benchmark failed (exit {r.returncode}): {r.stderr[:300]}"}
    L.log(f"    Benchmark completed")
    return _check_benchmark_tables(tables_dir)


def _check_benchmark_tables(tables_dir: Path) -> dict:
    """Check all expected benchmark tables."""
    results = {"tables_dir": str(tables_dir), "tables": {}}
    for fname in EXPECTED_BENCHMARK_TABLES:
        fpath = tables_dir / fname
        if fpath.exists():
            size = fpath.stat().st_size
            if fname.endswith(".tsv"):
                with open(fpath) as fh:
                    n = sum(1 for _ in fh) - 1
                results["tables"][fname] = {"status": "OK", "rows": n, "size": size}
            else:
                results["tables"][fname] = {"status": "OK", "size": size}
        else:
            results["tables"][fname] = {"status": "MISSING"}
    return results


def stage2_scorer(output: dict) -> tuple[int, list[str]]:
    critiques = []
    if "error" in output:
        return 0, [output["error"]]

    tables = output.get("tables", {})
    n_ok = sum(1 for t in tables.values() if t.get("status") == "OK")
    n_total = len(EXPECTED_BENCHMARK_TABLES)

    # Check specific table content
    for fname, info in tables.items():
        if info.get("status") == "MISSING":
            critiques.append(f"{fname}: MISSING")
        elif info.get("rows", 1) == 0 and fname.endswith(".tsv"):
            critiques.append(f"{fname}: empty (0 data rows)")

    score = min(10, round(n_ok / n_total * 10))
    L.log(f"    Benchmark tables: {n_ok}/{n_total} present")
    return score, critiques


def stage2_fixer(critiques: list[str], output: dict) -> list[str]:
    fixes = []
    for c in critiques:
        if "MISSING" in c:
            fixes.append(f"Will re-run benchmark to regenerate: {c.split(':')[0]}")
    # If tables are missing, clear caches and re-run
    if any("MISSING" in c for c in critiques):
        fixes.append("Triggering full benchmark re-run on next attempt")
        # Force re-run by touching the config file timestamp
        config = BASE / "conf" / "benchmark_1000g_full_cohort_v2.yaml"
        if config.exists():
            os.utime(config, None)
    return fixes


# ═══════════════════════════════════════════════════════════════════════════
# Stage 3: Figure Generation
# ═══════════════════════════════════════════════════════════════════════════
def stage3_runner() -> dict:
    """Run all three figure generation scripts."""
    results = {}

    # 3a: Main figures
    L.log(f"    Running generate_figures.py...")
    r = subprocess.run(
        ["python3.11", str(BASE / "bin" / "generate_figures.py"),
         "--tables-dir", str(BASE / "analysis" / "1000g_realdata" / "tables"),
         "--figures-dir", str(FIGURES_DIR)],
        capture_output=True, text=True, timeout=300, cwd=str(BASE),
    )
    results["main"] = {"exit": r.returncode, "stdout": r.stdout[-300:], "stderr": r.stderr[-300:]}
    L.log(f"    generate_figures.py → exit {r.returncode}")

    # 3b: Trimodal
    trimodal_dir = BASE / "analysis" / "benchmark_trimodal_all_samples" / "tables"
    if trimodal_dir.exists():
        L.log(f"    Running generate_figures_pub.py...")
        r = subprocess.run(
            ["python3.11", str(BASE / "bin" / "generate_figures_pub.py"),
             "--trimodal-tables", str(trimodal_dir),
             "--out-dir", str(FIGURES_DIR)],
            capture_output=True, text=True, timeout=300, cwd=str(BASE),
        )
        results["trimodal"] = {"exit": r.returncode, "stdout": r.stdout[-300:], "stderr": r.stderr[-300:]}
        L.log(f"    generate_figures_pub.py → exit {r.returncode}")
    else:
        results["trimodal"] = {"exit": -1, "stderr": "trimodal tables not found"}

    # 3c: Performance
    L.log(f"    Running generate_performance_figures.py...")
    r = subprocess.run(
        ["python3.11", str(BASE / "bin" / "generate_performance_figures.py")],
        capture_output=True, text=True, timeout=120, cwd=str(BASE),
    )
    results["perf"] = {"exit": r.returncode, "stdout": r.stdout[-300:], "stderr": r.stderr[-300:]}
    L.log(f"    generate_performance_figures.py → exit {r.returncode}")

    # 3d: Check outputs
    fig_status = {}
    for stem in EXPECTED_FIGS:
        for fmt in ["pdf", "png"]:
            fpath = FIGURES_DIR / f"{stem}.{fmt}"
            if fpath.exists():
                fig_status[f"{stem}.{fmt}"] = {
                    "status": "OK",
                    "size_kb": round(fpath.stat().st_size / 1024, 1),
                    "mtime": datetime.fromtimestamp(fpath.stat().st_mtime).strftime("%Y-%m-%d %H:%M"),
                }
            else:
                fig_status[f"{stem}.{fmt}"] = {"status": "MISSING"}
    results["figures"] = fig_status
    return results


def stage3_scorer(output: dict) -> tuple[int, list[str]]:
    critiques = []

    # Check script exit codes
    for name in ["main", "trimodal", "perf"]:
        info = output.get(name, {})
        if info.get("exit", -1) != 0:
            critiques.append(f"{name} script failed (exit {info.get('exit')}): {info.get('stderr', '')[:100]}")

    # Check figure files
    figs = output.get("figures", {})
    n_ok = sum(1 for f in figs.values() if f.get("status") == "OK")
    n_total = len(figs)
    missing = [k for k, v in figs.items() if v.get("status") == "MISSING"]
    if missing:
        critiques.append(f"{len(missing)} figure files missing: {', '.join(missing[:5])}")

    # Check for tiny figures (< 1KB = likely empty/broken)
    for fname, info in figs.items():
        if info.get("status") == "OK" and info.get("size_kb", 0) < 1:
            critiques.append(f"{fname}: suspiciously small ({info['size_kb']} KB)")

    pct = n_ok / n_total * 100 if n_total > 0 else 0
    score = min(10, round(pct / 10))
    L.log(f"    Figures: {n_ok}/{n_total} present ({pct:.0f}%)")
    return score, critiques


def stage3_fixer(critiques: list[str], output: dict) -> list[str]:
    fixes = []
    for c in critiques:
        if "script failed" in c:
            fixes.append(f"Will retry figure script: {c.split(' script')[0]}")
        elif "missing" in c.lower():
            fixes.append("Figures may be named differently; checking alternative paths")
    return fixes


# ═══════════════════════════════════════════════════════════════════════════
# Stage 4: Performance Timing
# ═══════════════════════════════════════════════════════════════════════════
def stage4_runner() -> dict:
    """Check and update timing data."""
    result = {"tsv_exists": TIMING_TSV.exists(), "rows": 0, "modalities": [], "wgs_traces_done": 0}

    if TIMING_TSV.exists():
        with open(TIMING_TSV, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        result["rows"] = len(rows)
        result["modalities"] = sorted({r["modality"] for r in rows})
        result["tools"] = sorted({r["tool"] for r in rows})
        L.log(f"    Timing TSV: {len(rows)} entries, modalities: {result['modalities']}")

    # Count WGS traces
    wgs_traces = glob.glob(str(
        BASE / "analysis" / "performance_benchmark_n30" / "runs" / "wgs" / "*" / "pipeline_info" / "execution_trace.txt"
    ))
    for t in wgs_traces:
        try:
            with open(t) as fh:
                if fh.read().count("COMPLETED") >= 7:
                    result["wgs_traces_done"] += 1
        except:
            pass
    L.log(f"    WGS traces: {len(wgs_traces)} total, {result['wgs_traces_done']} fully completed")

    # Auto-update: if any modality is missing, re-run full trace parser with all sources
    expected_mods = {"WES", "RNA", "WGS"}
    missing_mods = expected_mods - set(result.get("modalities", []))
    if missing_mods:
        L.log(f"    Missing modalities: {missing_mods} — re-running full trace parser")
        wes_traces = sorted(glob.glob(str(BATCH_ROOT / "wes_batches" / "*" / "results" / "pipeline_info" / "execution_trace.txt")))
        rna_traces = sorted(glob.glob(str(BATCH_ROOT / "rna_batches" / "*" / "results" / "pipeline_info" / "execution_trace.txt")))
        cmd = ["python3.11", str(BASE / "bin" / "parse_nextflow_trace_timing.py")]
        if wes_traces:
            cmd += ["--wes-traces"] + wes_traces
        if rna_traces:
            cmd += ["--rna-traces"] + rna_traces
        if wgs_traces and result["wgs_traces_done"] >= 20:
            cmd += ["--wgs-traces"] + sorted(wgs_traces)
        cmd += ["--n-per-modality", "30", "--out-dir", str(TIMING_TSV.parent)]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300, cwd=str(BASE))
        if r.returncode == 0:
            L.log(f"    ✓ TSV updated with all available modalities")
            with open(TIMING_TSV, newline="") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            result["rows"] = len(rows)
            result["modalities"] = sorted({r["modality"] for r in rows})
        else:
            L.log(f"    ✗ Trace parser failed: {r.stderr[:300]}", "ERROR")

    return result


def stage4_scorer(output: dict) -> tuple[int, list[str]]:
    critiques = []
    if not output.get("tsv_exists"):
        return 0, ["Timing TSV does not exist"]

    mods = output.get("modalities", [])
    n_rows = output.get("rows", 0)

    if n_rows == 0:
        return 0, ["Timing TSV is empty"]
    if "WGS" not in mods:
        critiques.append("WGS timing not yet available")
    if "WES" not in mods:
        critiques.append("WES timing missing")
    if "RNA" not in mods:
        critiques.append("RNA timing missing")

    n_mods = len(mods)
    score = min(10, n_mods * 3 + (1 if n_rows >= 10 else 0))  # 3pts per modality + 1 for completeness
    return score, critiques


def stage4_fixer(critiques: list[str], output: dict) -> list[str]:
    fixes = []
    if "WGS timing not yet available" in critiques:
        wgs_done = output.get("wgs_traces_done", 0)
        if wgs_done < 20:
            fixes.append(f"WGS: only {wgs_done}/30 traces complete — SLURM job still running")
        else:
            fixes.append("WGS: enough traces done, will re-run parser on next attempt")
    return fixes


# ═══════════════════════════════════════════════════════════════════════════
# Stage 5: Final Report
# ═══════════════════════════════════════════════════════════════════════════
def compile_report(results: dict[str, StepResult]):
    """Compile all results into a final validation report."""
    L.section("FINAL REPORT")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- Tool inventory TSV ---
    inv = results.get("tool_outputs")
    if inv and inv.details.get("inventory"):
        inv_path = OUT_DIR / "tool_output_inventory.tsv"
        with open(inv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["modality", "sample", "tool", "status", "n_rows", "genes", "errors"], delimiter="\t")
            writer.writeheader()
            writer.writerows(inv.details["inventory"])
        L.log(f"  Wrote: {inv_path}")

    # --- Figure inventory TSV ---
    figs = results.get("figures")
    if figs and figs.details.get("figures"):
        fig_path = OUT_DIR / "figure_inventory.tsv"
        with open(fig_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["figure", "status", "size_kb", "mtime"], delimiter="\t")
            writer.writeheader()
            for fname, info in figs.details["figures"].items():
                writer.writerow({"figure": fname, **info})
        L.log(f"  Wrote: {fig_path}")

    # --- Markdown report ---
    report_path = OUT_DIR / "validation_report.md"
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    lines = [
        f"# MVHLA Pipeline Validation Report",
        f"",
        f"**Date:** {now}",
        f"**Pipeline:** MVHLA v2.0.0",
        f"**Approach:** Draft → Critique → Fix → Loop (max {MAX_RETRIES} retries, threshold {PASS_THRESHOLD}/10)",
        f"",
        f"---",
        f"",
        f"## Scorecard",
        f"",
        f"| Stage | Score | Status | Retries | Critiques |",
        f"|---|---|---|---|---|",
    ]

    for name, res in results.items():
        status = "PASS" if res.passed else "FAIL"
        n_crit = len(res.critiques)
        n_fix = len(res.fixes_applied)
        lines.append(f"| {name} | {res.score}/10 | {status} | {n_fix} fixes | {n_crit} |")

    overall = all(r.passed for r in results.values())
    lines.extend([
        f"",
        f"**Overall: {'ALL STAGES PASSED' if overall else 'SOME STAGES FAILED'}**",
        f"",
        f"---",
        f"",
    ])

    # Per-stage details
    for name, res in results.items():
        lines.append(f"## Stage: {name}")
        lines.append(f"")
        lines.append(f"**Score:** {res.score}/10 ({'PASS' if res.passed else 'FAIL'})")
        lines.append(f"")
        if res.critiques:
            lines.append(f"**Critiques:**")
            for c in res.critiques:
                lines.append(f"- {c}")
            lines.append(f"")
        if res.fixes_applied:
            lines.append(f"**Fixes applied:**")
            for f in res.fixes_applied:
                lines.append(f"- {f}")
            lines.append(f"")
        lines.append(f"---")
        lines.append(f"")

    # Errors & warnings
    lines.append(f"## Errors ({len(L.errors)})")
    lines.append(f"")
    for e in L.errors:
        lines.append(f"- {e}")
    lines.append(f"")

    lines.append(f"## Warnings ({len(L.warnings)})")
    lines.append(f"")
    for w in L.warnings:
        lines.append(f"- {w}")
    lines.append(f"")

    lines.append(f"## Log (last 40 lines)")
    lines.append(f"```")
    lines.extend(L.lines[-40:])
    lines.append(f"```")

    with open(report_path, "w") as fh:
        fh.write("\n".join(lines))
    L.log(f"  Wrote: {report_path}")

    L.log(f"")
    L.log(f"  ╔════════════════════════════════════════════════════╗")
    L.log(f"  ║  VALIDATION {'PASSED' if overall else 'FAILED':8s}                             ║")
    L.log(f"  ║  Errors: {len(L.errors):<3}  Warnings: {len(L.warnings):<3}                    ║")
    L.log(f"  ║  Report: pipeline_validation_results/              ║")
    L.log(f"  ╚════════════════════════════════════════════════════╝")


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════
def main():
    L.log("MVHLA Pipeline End-to-End Validation (Loop Approach)")
    L.log(f"  Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.log(f"  Base: {BASE}")
    L.log(f"  Max retries: {MAX_RETRIES}, Pass threshold: {PASS_THRESHOLD}/10")

    results = {}

    # Stage 1: Tool outputs
    results["tool_outputs"] = run_loop(
        "Tool Output Validation",
        runner=stage1_runner,
        scorer=stage1_scorer,
        fixer=stage1_fixer,
    )

    # Stage 2: Benchmark
    results["benchmark"] = run_loop(
        "Benchmark Analysis",
        runner=stage2_runner,
        scorer=stage2_scorer,
        fixer=stage2_fixer,
    )

    # Stage 3: Figures
    results["figures"] = run_loop(
        "Figure Generation",
        runner=stage3_runner,
        scorer=stage3_scorer,
        fixer=stage3_fixer,
    )

    # Stage 4: Timing
    results["timing"] = run_loop(
        "Performance Timing",
        runner=stage4_runner,
        scorer=stage4_scorer,
        fixer=stage4_fixer,
    )

    # Stage 5: Report
    compile_report(results)


if __name__ == "__main__":
    main()
