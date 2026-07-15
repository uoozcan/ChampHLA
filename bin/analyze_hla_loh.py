#!/usr/bin/env python3.11
"""Summarize HLA LOH-like patterns from live /scratch results with single-allele QC."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


CLASS_I_GENES = ("A", "B", "C")
CLASSICAL_GENES = {
    "A",
    "B",
    "C",
    "DRB1",
    "DRB3",
    "DRB4",
    "DRB5",
    "DQA1",
    "DQB1",
    "DPA1",
    "DPB1",
}
DEFAULT_BENCHMARK_TABLES = Path(
    "/scratch/project_2008084/pihla-publish/analysis/benchmark_trimodal_all_samples/tables"
)
DEFAULT_VENEX_RESULTS = Path("/scratch/project_2008084/hla_calibration/venex/results/wgs")
DEFAULT_VENEX_WORK = Path("/scratch/project_2008084/hla_calibration/venex/work/wgs")
DEFAULT_BENCHMARK_ROOTS = [
    Path("/scratch/project_2008084/hla_calibration/wgs/results"),
    Path("/scratch/project_2008084/hla_calibration/wgs_batches"),
    Path("/scratch/project_2008084/hla_calibration/wes/results"),
    Path("/scratch/project_2008084/hla_calibration/wes_batches"),
    Path("/scratch/project_2008084/hla_calibration/rna/results"),
    Path("/scratch/project_2008084/hla_calibration/rna_3sample/results"),
]
DEFAULT_OUTPUT_DIR = Path("/scratch/project_2008084/pihla-publish/analysis/loh_analysis")
LOH_FREQ_MAJOR_THRESHOLD = 0.80
LOH_FREQ_BALANCED_MIN = 0.35
LOH_FREQ_BALANCED_MAX = 0.65
LOH_HET_VARIANT_MIN = 20
FALSE_DUPLICATE_RATE_CAUTION = 0.50


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-tables", type=Path, default=DEFAULT_BENCHMARK_TABLES)
    parser.add_argument("--venex-results", type=Path, default=DEFAULT_VENEX_RESULTS)
    parser.add_argument("--venex-work", type=Path, default=DEFAULT_VENEX_WORK)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_comment_aware_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    if not lines:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def write_tsv(path: Path, fieldnames: List[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def normalize_gene_name(value: str) -> str:
    gene = value.strip()
    if gene.startswith("HLA-"):
        gene = gene[4:]
    elif gene.startswith("HLA_"):
        gene = gene[4:]
    if gene.startswith("HLA"):
        gene = gene[3:]
    return gene.strip("_")


def normalize_sample_from_path(source_file: str) -> Optional[str]:
    parts = Path(source_file).parts
    for idx, part in enumerate(parts):
        if part == "results" and idx + 1 < len(parts):
            return parts[idx + 1]
    return None


def find_arcashla_json_for_source(source_file: str) -> Optional[Path]:
    source = Path(source_file)
    sample = normalize_sample_from_path(source_file)
    if not sample:
        return None
    direct = source.with_suffix(".json")
    if direct.exists():
        return direct
    parent = source.parent
    candidate = parent / f"{sample}_arcashla.json"
    if candidate.exists():
        return candidate
    matches = list(parent.glob("*.json"))
    return matches[0] if matches else None


def build_arcashla_raw_single_map(benchmark_rows: List[Dict[str, str]]) -> Dict[Tuple[str, str, str, str], bool]:
    mapping: Dict[Tuple[str, str, str, str], bool] = {}
    cache: Dict[Path, Dict[str, int]] = {}
    for row in benchmark_rows:
        if row.get("tool") != "ArcasHLA":
            continue
        json_path = find_arcashla_json_for_source(row.get("source_file", ""))
        if not json_path or not json_path.exists():
            continue
        if json_path not in cache:
            try:
                raw = json.loads(json_path.read_text())
            except Exception:
                cache[json_path] = {}
            else:
                gene_lengths = {}
                for gene_name, alleles in raw.items():
                    gene = normalize_gene_name(gene_name)
                    if gene in CLASSICAL_GENES and isinstance(alleles, list):
                        gene_lengths[gene] = len(alleles)
                cache[json_path] = gene_lengths
        sample = row["sample"]
        modality = row["modality"]
        gene = row["gene"]
        length = cache[json_path].get(gene)
        if length == 1:
            mapping[(sample, modality, "ArcasHLA", gene)] = True
        elif length is not None:
            mapping[(sample, modality, "ArcasHLA", gene)] = False
    return mapping


def classify_benchmark_row(
    row: Dict[str, str], raw_single_map: Dict[Tuple[str, str, str, str], bool]
) -> str:
    key = (row["sample"], row["modality"], row["tool"], row["gene"])
    if raw_single_map.get(key):
        return "single_allele_raw"
    allele1 = row.get("allele1", "").strip()
    allele2 = row.get("allele2", "").strip()
    truth1 = row.get("truth_allele1", "").strip()
    truth2 = row.get("truth_allele2", "").strip()
    if not allele1 and not allele2:
        return "no_call"
    if allele1 and (not allele2 or allele2 == "-"):
        return "single_allele_raw"
    if allele1 == allele2:
        if truth1 == truth2 and allele1 == truth1:
            return "true_homozygous_supported"
        return "possible_allele_dropout"
    return "heterozygous"


def benchmark_qc_rows(
    benchmark_rows: List[Dict[str, str]],
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], Dict[Tuple[str, str, str], Dict[str, float]]]:
    raw_single_map = build_arcashla_raw_single_map(benchmark_rows)
    detailed_rows: List[Dict[str, object]] = []
    summary_counter: Dict[Tuple[str, str, str], Counter] = defaultdict(Counter)
    duplicate_counter: Dict[Tuple[str, str, str], Counter] = defaultdict(Counter)

    for row in benchmark_rows:
        gene = row.get("gene")
        if gene not in CLASS_I_GENES:
            continue
        qc_class = classify_benchmark_row(row, raw_single_map)
        allele1 = row.get("allele1", "")
        allele2 = row.get("allele2", "")
        truth1 = row.get("truth_allele1", "")
        truth2 = row.get("truth_allele2", "")
        duplicate_call = bool(allele1 and allele2 and allele1 == allele2)
        truth_homozygous = bool(truth1 and truth2 and truth1 == truth2)

        detailed_rows.append(
            {
                "sample": row["sample"],
                "modality": row["modality"],
                "tool": row["tool"],
                "gene": gene,
                "truth_allele1": truth1,
                "truth_allele2": truth2,
                "allele1": allele1,
                "allele2": allele2,
                "call_status": row.get("call_status", ""),
                "correct_status": row.get("correct_status", ""),
                "qc_class": qc_class,
                "duplicate_call": int(duplicate_call),
                "truth_homozygous": int(truth_homozygous),
                "source_file": row.get("source_file", ""),
            }
        )

        key = (row["tool"], row["modality"], gene)
        summary_counter[key]["rows"] += 1
        summary_counter[key][qc_class] += 1
        if duplicate_call:
            duplicate_counter[key]["duplicate_rows"] += 1
            if qc_class == "true_homozygous_supported":
                duplicate_counter[key]["duplicate_matches_homo_truth"] += 1
            elif qc_class in {"possible_allele_dropout", "single_allele_raw"}:
                duplicate_counter[key]["duplicate_not_homo_truth"] += 1

    summary_rows: List[Dict[str, object]] = []
    caution_summary: Dict[Tuple[str, str, str], Dict[str, float]] = {}
    for key in sorted(summary_counter):
        tool, modality, gene = key
        counts = summary_counter[key]
        dup = duplicate_counter.get(key, Counter())
        duplicate_rows = dup.get("duplicate_rows", 0)
        duplicate_bad = dup.get("duplicate_not_homo_truth", 0)
        false_duplicate_rate = (
            duplicate_bad / duplicate_rows if duplicate_rows else 0.0
        )
        caution_summary[key] = {
            "duplicate_rows": float(duplicate_rows),
            "duplicate_not_homo_truth": float(duplicate_bad),
            "false_duplicate_rate": false_duplicate_rate,
        }
        summary_rows.append(
            {
                "tool": tool,
                "modality": modality,
                "gene": gene,
                "rows": counts.get("rows", 0),
                "heterozygous": counts.get("heterozygous", 0),
                "true_homozygous_supported": counts.get("true_homozygous_supported", 0),
                "single_allele_raw": counts.get("single_allele_raw", 0),
                "possible_allele_dropout": counts.get("possible_allele_dropout", 0),
                "no_call": counts.get("no_call", 0),
                "duplicate_rows": duplicate_rows,
                "duplicate_not_homo_truth": duplicate_bad,
                "false_duplicate_rate": f"{false_duplicate_rate:.4f}",
            }
        )
    return detailed_rows, summary_rows, caution_summary


def parse_spechla_summary(path: Path) -> Dict[str, Dict[str, str]]:
    rows = read_comment_aware_tsv(path)
    if not rows:
        return {}
    row = rows[0]
    result = {}
    for gene in CLASS_I_GENES:
        allele1 = row.get(f"HLA_{gene}_1", "").strip()
        allele2 = row.get(f"HLA_{gene}_2", "").strip()
        result[gene] = {"allele1": allele1, "allele2": allele2}
    return result


def parse_freq_file(path: Path) -> Dict[str, object]:
    allele_freqs: Dict[str, float] = {}
    het_variant_count: Optional[int] = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# The number of heterozygous variant is"):
                match = re.search(r"(\d+)$", line)
                if match:
                    het_variant_count = int(match.group(1))
                continue
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                allele_freqs[parts[0]] = float(parts[1])
    allele1_freq = allele_freqs.get(path.stem.replace("_freq", "").replace("HLA_", "hla.allele.1.HLA_") + ".fasta")
    allele2_freq = allele_freqs.get(path.stem.replace("_freq", "").replace("HLA_", "hla.allele.2.HLA_") + ".fasta")
    # If keys were not synthesized exactly, fall back to first two sorted entries.
    if allele1_freq is None or allele2_freq is None:
        ordered = [allele_freqs[key] for key in sorted(allele_freqs)]
        if len(ordered) >= 2:
            allele1_freq = ordered[0]
            allele2_freq = ordered[1]
        elif len(ordered) == 1:
            allele1_freq = ordered[0]
            allele2_freq = None
    return {
        "allele1_freq": allele1_freq,
        "allele2_freq": allele2_freq,
        "heterozygous_variant_count": het_variant_count,
    }


def parse_low_depth_genes(path: Path) -> set[str]:
    genes = set()
    if not path.exists():
        return genes
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) >= 1 and parts[0].startswith("HLA_"):
                genes.add(parts[0].replace("HLA_", ""))
    return genes


def parse_hla_result_details(path: Path) -> Dict[str, Dict[str, str]]:
    details: Dict[str, Dict[str, str]] = {}
    if not path.exists():
        return details
    for row in read_comment_aware_tsv(path):
        gene_field = row.get("Gene", "")
        match = re.match(r"HLA_([A-Z0-9]+)_(1|2)$", gene_field)
        if not match:
            continue
        gene = match.group(1)
        allele_idx = match.group(2)
        if gene not in CLASS_I_GENES:
            continue
        details.setdefault(gene, {})
        prefix = f"allele{allele_idx}"
        details[gene][f"{prefix}_best"] = row.get("G_best", "")
        detail_field = row.get("details:allele;Score;Caucasian;Black;Asian", "")
        parts = detail_field.split(";")
        if len(parts) >= 2:
            details[gene][f"{prefix}_score"] = parts[1]
    return details


def find_venex_work_sample_dir(sample: str, venex_work: Path) -> Optional[Path]:
    candidates = list(venex_work.glob(f"wgs_{sample}/**/{sample}/{sample}"))
    if not candidates:
        return None
    # Prefer the directory that actually contains HLA freq files (multiple Nextflow runs may exist)
    for candidate in candidates:
        if any(candidate.glob("HLA_*_freq.txt")):
            return candidate
    return candidates[0]


def parse_published_spechla_samples(venex_results: Path) -> List[Tuple[str, Path]]:
    pairs = []
    for summary in sorted(venex_results.glob("*/results/*/spechla/*_spechla.txt")):
        sample = summary.stem.replace("_spechla", "")
        pairs.append((sample, summary))
    return pairs


def classify_venex_gene(
    gene_call: Dict[str, str],
    freq_info: Dict[str, object],
    low_depth_genes: set[str],
    caution: Dict[str, float],
    gene: str,
) -> Tuple[str, str]:
    allele1 = gene_call.get("allele1", "")
    allele2 = gene_call.get("allele2", "")
    false_duplicate_rate = caution.get("false_duplicate_rate", 0.0)
    if not allele1 and not allele2:
        return "no_call", "insufficient_evidence"
    if allele1 and (not allele2 or allele2 == "-"):
        return "single_allele_raw", "insufficient_evidence"
    if allele1 == allele2:
        if false_duplicate_rate >= FALSE_DUPLICATE_RATE_CAUTION:
            return "possible_allele_dropout", "likely_allele_dropout"
        return "duplicated_call_supported", "putative_homozygous"

    allele1_freq = freq_info.get("allele1_freq")
    allele2_freq = freq_info.get("allele2_freq")
    het_var_count = freq_info.get("heterozygous_variant_count")
    if allele1_freq is None or allele2_freq is None:
        return "heterozygous", "insufficient_evidence"

    major = max(float(allele1_freq), float(allele2_freq))
    minor = min(float(allele1_freq), float(allele2_freq))
    gene_low_depth = gene in low_depth_genes

    if gene_low_depth:
        return "heterozygous", "technical_caution_low_depth"
    if (
        major >= LOH_FREQ_MAJOR_THRESHOLD
        and het_var_count is not None
        and int(het_var_count) >= LOH_HET_VARIANT_MIN
    ):
        return "heterozygous", "candidate_loh"
    if LOH_FREQ_BALANCED_MIN <= minor <= LOH_FREQ_BALANCED_MAX and LOH_FREQ_BALANCED_MIN <= major <= LOH_FREQ_BALANCED_MAX:
        return "heterozygous", "balanced_heterozygous"
    return "heterozygous", "allelic_imbalance_review"


def build_venex_rows(
    venex_results: Path,
    venex_work: Path,
    caution_summary: Dict[Tuple[str, str, str], Dict[str, float]],
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    output_rows: List[Dict[str, object]] = []
    counters = Counter()
    for sample, summary_path in parse_published_spechla_samples(venex_results):
        counters["published_spechla_samples"] += 1
        calls = parse_spechla_summary(summary_path)
        work_dir = find_venex_work_sample_dir(sample, venex_work)
        if work_dir:
            counters["samples_with_work_dir"] += 1
        low_depth_genes = parse_low_depth_genes(work_dir / "low_depth.bed") if work_dir else set()
        details = parse_hla_result_details(work_dir / "hla.result.details.txt") if work_dir else {}
        for gene in CLASS_I_GENES:
            call = calls.get(gene, {"allele1": "", "allele2": ""})
            freq_path = work_dir / f"HLA_{gene}_freq.txt" if work_dir else None
            has_freq = bool(freq_path and freq_path.exists())
            if has_freq:
                counters["gene_rows_with_freq"] += 1
            freq_info = parse_freq_file(freq_path) if has_freq else {
                "allele1_freq": None,
                "allele2_freq": None,
                "heterozygous_variant_count": None,
            }
            caution = caution_summary.get(("SpecHLA", "wgs", gene), {})
            qc_flag, loh_status = classify_venex_gene(call, freq_info, low_depth_genes, caution, gene)
            counters[f"loh_status::{loh_status}"] += 1
            output_rows.append(
                {
                    "cohort": "VENEX_WGS",
                    "sample": sample,
                    "gene": gene,
                    "tool": "SpecHLA",
                    "allele1": call.get("allele1", ""),
                    "allele2": call.get("allele2", ""),
                    "qc_flag": qc_flag,
                    "loh_status": loh_status,
                    "allele1_freq": (
                        f"{freq_info['allele1_freq']:.3f}"
                        if freq_info.get("allele1_freq") is not None
                        else ""
                    ),
                    "allele2_freq": (
                        f"{freq_info['allele2_freq']:.3f}"
                        if freq_info.get("allele2_freq") is not None
                        else ""
                    ),
                    "major_allele_fraction": (
                        f"{max(float(freq_info['allele1_freq']), float(freq_info['allele2_freq'])):.3f}"
                        if freq_info.get("allele1_freq") is not None and freq_info.get("allele2_freq") is not None
                        else ""
                    ),
                    "heterozygous_variant_count": freq_info.get("heterozygous_variant_count", ""),
                    "low_depth_flag": int(gene in low_depth_genes),
                    "details_allele1_score": details.get(gene, {}).get("allele1_score", ""),
                    "details_allele2_score": details.get(gene, {}).get("allele2_score", ""),
                    "work_dir": str(work_dir) if work_dir else "",
                    "published_result_file": str(summary_path),
                    "freq_file": str(freq_path) if has_freq else "",
                    "benchmark_false_duplicate_rate": (
                        f"{caution.get('false_duplicate_rate', 0.0):.4f}" if caution else ""
                    ),
                    "benchmark_duplicate_rows": int(caution.get("duplicate_rows", 0)) if caution else 0,
                }
            )
    return output_rows, dict(counters)


def write_metadata(path: Path, payload: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def main() -> int:
    args = parse_args()
    benchmark_rows = read_tsv(args.benchmark_tables / "harmonized_benchmark_rows.tsv")
    benchmark_qc, benchmark_summary, caution_summary = benchmark_qc_rows(benchmark_rows)
    venex_rows, venex_counts = build_venex_rows(
        args.venex_results,
        args.venex_work,
        caution_summary,
    )

    output_dir = args.output_dir
    write_tsv(
        output_dir / "tool_homozygosity_qc.tsv",
        [
            "sample",
            "modality",
            "tool",
            "gene",
            "truth_allele1",
            "truth_allele2",
            "allele1",
            "allele2",
            "call_status",
            "correct_status",
            "qc_class",
            "duplicate_call",
            "truth_homozygous",
            "source_file",
        ],
        benchmark_qc,
    )
    write_tsv(
        output_dir / "benchmark_homozygosity_error_summary.tsv",
        [
            "tool",
            "modality",
            "gene",
            "rows",
            "heterozygous",
            "true_homozygous_supported",
            "single_allele_raw",
            "possible_allele_dropout",
            "no_call",
            "duplicate_rows",
            "duplicate_not_homo_truth",
            "false_duplicate_rate",
        ],
        benchmark_summary,
    )
    write_tsv(
        output_dir / "sample_loh_candidates.tsv",
        [
            "cohort",
            "sample",
            "gene",
            "tool",
            "allele1",
            "allele2",
            "qc_flag",
            "loh_status",
            "allele1_freq",
            "allele2_freq",
            "major_allele_fraction",
            "heterozygous_variant_count",
            "low_depth_flag",
            "details_allele1_score",
            "details_allele2_score",
            "work_dir",
            "published_result_file",
            "freq_file",
            "benchmark_false_duplicate_rate",
            "benchmark_duplicate_rows",
        ],
        venex_rows,
    )
    write_metadata(
        output_dir / "loh_analysis_metadata.json",
        {
            "authoritative_root": "/scratch/project_2008084",
            "benchmark_tables": str(args.benchmark_tables),
            "venex_results": str(args.venex_results),
            "venex_work": str(args.venex_work),
            "output_dir": str(output_dir),
            "class_i_genes": list(CLASS_I_GENES),
            "loh_thresholds": {
                "major_allele_fraction_candidate_loh": LOH_FREQ_MAJOR_THRESHOLD,
                "balanced_min": LOH_FREQ_BALANCED_MIN,
                "balanced_max": LOH_FREQ_BALANCED_MAX,
                "heterozygous_variant_min": LOH_HET_VARIANT_MIN,
                "false_duplicate_rate_caution": FALSE_DUPLICATE_RATE_CAUTION,
            },
            "venex_summary": venex_counts,
        },
    )

    print(output_dir)
    print(f"benchmark_qc_rows\t{len(benchmark_qc)}")
    print(f"benchmark_summary_rows\t{len(benchmark_summary)}")
    print(f"venex_gene_rows\t{len(venex_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
def read_comment_aware_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    if not lines:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))
