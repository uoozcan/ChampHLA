#!/usr/bin/env python3
"""Run the real-data 1000 Genomes HLA benchmark workflow."""

import argparse
import csv
import json
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path

import yaml

import hla_benchmark as hb
from build_1000g_benchmark_manifests import (
    build_cohort_manifest,
    load_sequencing_manifest_rows,
    load_truth_manifest_rows,
    parse_population_map,
    write_tsv,
)

REQUIRED_MODALITIES = ["wgs", "wes", "rnaseq"]


def string_value(value):
    return "" if value is None else str(value)


def parse_args():
    parser = argparse.ArgumentParser(description="Run the 1000 Genomes real-data HLA benchmark.")
    parser.add_argument("--config", required=True, help="YAML config describing manifests, truth, and result globs.")
    parser.add_argument("--output-dir", required=True, help="Final output directory for the real-data benchmark.")
    parser.add_argument("--weight-alpha", type=float, default=0.7,
                        help="Weight on base_reliability in final_weight formula (default: 0.7).")
    parser.add_argument("--weight-beta", type=float, default=0.3,
                        help="Weight on effective_confidence in final_weight formula (default: 0.3).")
    parser.add_argument("--benchmark-mode", default="",
                        help="Override benchmark.mode from config (legacy_heuristic, probabilistic_recalibrated, bayesian_shrinkage).")
    return parser.parse_args()


def load_yaml(path):
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def save_yaml(path, payload):
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, default_flow_style=False)


def load_rows(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_truth_supported_loci(config, cohort_rows):
    configured = hb.parse_gene_list(config.get("truth", {}).get("supported_loci") or config.get("benchmark", {}).get("genes"))
    if configured:
        return configured
    loci = []
    for row in cohort_rows:
        if row.get("include") == "1":
            for gene in hb.parse_gene_list(row.get("truth_supported_loci", "")):
                if gene not in loci:
                    loci.append(gene)
    return sorted(loci, key=hb.gene_sort_key)


def ensure_manifests(config, base_output_dir):
    manifests_cfg = config.setdefault("manifests", {})
    truth_manifest = manifests_cfg.get("truth_manifest")
    sequencing_manifest = manifests_cfg.get("sequencing_manifest")
    cohort_manifest = manifests_cfg.get("cohort_manifest")
    if truth_manifest and sequencing_manifest and cohort_manifest:
        return Path(truth_manifest), Path(sequencing_manifest), Path(cohort_manifest)

    sequencing_source = manifests_cfg.get("sequencing_source") or manifests_cfg.get("sequencing_manifest_source")
    if not sequencing_source:
        raise SystemExit("1000 Genomes benchmark requires manifests.truth_manifest/manifests.sequencing_manifest/manifests.cohort_manifest or manifests.sequencing_source")
    population_manifest = manifests_cfg.get("population_manifest", "")
    generated_dir = base_output_dir / "manifests"
    generated_dir.mkdir(parents=True, exist_ok=True)
    supported_loci = hb.parse_gene_list(config.get("truth", {}).get("supported_loci") or "A,B,C,DRB1,DQB1")
    population_map = parse_population_map(population_manifest)
    truth_rows = load_truth_manifest_rows(config["truth"]["path"], supported_loci, config["truth"].get("source", ""), string_value(config["truth"].get("acquisition_date", "")), population_map)
    sequencing_rows = load_sequencing_manifest_rows(sequencing_source, population_map)
    # Determine required modalities from config runs (if all three are present, use default behaviour)
    config_modalities = list({run["modality"] for run in config.get("runs", []) if "modality" in run})
    required_modalities = config_modalities if config_modalities else None
    cohort_rows = build_cohort_manifest(truth_rows, sequencing_rows, supported_loci, required_modalities=required_modalities)
    truth_manifest = generated_dir / "truth_manifest.tsv"
    sequencing_manifest = generated_dir / "sequencing_manifest.tsv"
    cohort_manifest = generated_dir / "cohort_manifest.tsv"
    write_tsv(truth_manifest, truth_rows, ["sample", "population", "truth_source", "acquisition_date", "truth_supported_loci", "truth_gene_count"])
    write_tsv(sequencing_manifest, sequencing_rows, ["sample", "population", "modality", "data_locator", "available"])
    cohort_fields = ["sample", "population", "include"]
    if any(hb.clean_token(row.get("split", "")) for row in cohort_rows):
        cohort_fields.append("split")
    cohort_fields.extend(["truth_supported_loci", "wgs_available", "wes_available", "rnaseq_available", "excluded_reason"])
    write_tsv(cohort_manifest, cohort_rows, cohort_fields)
    manifests_cfg["truth_manifest"] = str(truth_manifest)
    manifests_cfg["sequencing_manifest"] = str(sequencing_manifest)
    manifests_cfg["cohort_manifest"] = str(cohort_manifest)
    return truth_manifest, sequencing_manifest, cohort_manifest


def build_filtered_truth_rows(truth_map, samples, supported_loci):
    rows = []
    for sample in sorted(samples):
        for gene in sorted(truth_map.get(sample, {}), key=hb.gene_sort_key):
            if supported_loci and gene not in supported_loci:
                continue
            pair = truth_map[sample][gene]
            rows.append({"sample": sample, "gene": gene, "allele1": pair[0], "allele2": pair[1]})
    return rows


def write_filtered_truth(path, truth_rows):
    write_tsv(path, truth_rows, ["sample", "gene", "allele1", "allele2"])


def subset_config(config, truth_path, supported_loci):
    cfg = deepcopy(config)
    cfg.setdefault("truth", {})
    cfg["truth"]["path"] = str(truth_path)
    cfg["truth"]["supported_loci"] = supported_loci
    cfg.setdefault("benchmark", {})
    cfg["benchmark"].pop("runtime_weight_override", None)
    return cfg


def run_benchmark(config_payload, output_dir, weight_alpha=0.7, weight_beta=0.3, benchmark_mode="", population_manifest=""):
    with tempfile.TemporaryDirectory() as tmpdir:
        cfg_path = Path(tmpdir) / "benchmark.yaml"
        save_yaml(cfg_path, config_payload)
        cmd = [
            "python3",
            str(Path(__file__).with_name("hla_benchmark.py")),
            "--config",
            str(cfg_path),
            "--output-dir",
            str(output_dir),
            "--weight-alpha",
            str(weight_alpha),
            "--weight-beta",
            str(weight_beta),
        ]
        if benchmark_mode:
            cmd.extend(["--benchmark-mode", benchmark_mode])
        if population_manifest:
            cmd.extend(["--population-manifest", population_manifest])
        subprocess.run(cmd, check=True)


def copy_manifest_outputs(output_dir, truth_manifest, sequencing_manifest, cohort_manifest):
    tables_dir = Path(output_dir) / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(truth_manifest, tables_dir / "truth_manifest.tsv")
    shutil.copy2(sequencing_manifest, tables_dir / "sequencing_manifest.tsv")
    shutil.copy2(cohort_manifest, tables_dir / "cohort_manifest.tsv")


def build_population_counts(cohort_rows, included_only=True):
    counter = Counter()
    for row in cohort_rows:
        if included_only and row.get("include") != "1":
            continue
        counter[row.get("population", "") or "unknown"] += 1
    return dict(sorted(counter.items()))


def build_excluded_summary(cohort_rows):
    counter = Counter()
    for row in cohort_rows:
        if row.get("include") != "1":
            counter[row.get("excluded_reason", "unspecified") or "unspecified"] += 1
    return dict(sorted(counter.items()))


def build_modality_counts(sequencing_rows, cohort_rows):
    before = defaultdict(set)
    after = defaultdict(set)
    included = {row["sample"] for row in cohort_rows if row.get("include") == "1"}
    for row in sequencing_rows:
        if row.get("available") == "1":
            before[row["modality"]].add(row["sample"])
            if row["sample"] in included:
                after[row["modality"]].add(row["sample"])
    return {
        "before_filtering": {modality: len(before.get(modality, set())) for modality in REQUIRED_MODALITIES},
        "after_filtering": {modality: len(after.get(modality, set())) for modality in REQUIRED_MODALITIES},
    }


def expected_tools_by_modality(config):
    grouped = defaultdict(set)
    for run in config.get("runs", []):
        modality = hb.clean_token(run.get("modality", "")).lower()
        tool = hb.clean_token(run.get("tool", ""))
        if modality in REQUIRED_MODALITIES and tool:
            grouped[modality].add(tool)
    return {modality: sorted(grouped.get(modality, set())) for modality in REQUIRED_MODALITIES}


def build_tool_availability(config, harmonized_rows, included_samples):
    expected = expected_tools_by_modality(config)
    observed = Counter()
    for row in harmonized_rows:
        sample = hb.clean_token(row.get("sample", ""))
        modality = hb.clean_token(row.get("modality", "")).lower()
        tool = hb.clean_token(row.get("tool", ""))
        if sample in included_samples and modality in REQUIRED_MODALITIES and tool:
            observed[(sample, modality, tool)] += 1

    sample_rows = []
    summary_counter = Counter()
    available_tools = {modality: set() for modality in REQUIRED_MODALITIES}
    not_available_tools = {modality: set() for modality in REQUIRED_MODALITIES}
    for sample in sorted(included_samples):
        for modality in REQUIRED_MODALITIES:
            for tool in expected.get(modality, []):
                n_rows = observed.get((sample, modality, tool), 0)
                status = "available" if n_rows > 0 else "not_available"
                sample_rows.append(
                    {
                        "sample": sample,
                        "modality": modality,
                        "tool": tool,
                        "status": status,
                        "observed_rows": str(n_rows),
                    }
                )
                summary_counter[(modality, tool, status)] += 1
                if status == "available":
                    available_tools[modality].add(tool)
                else:
                    not_available_tools[modality].add(tool)

    summary_rows = []
    for key, count in sorted(summary_counter.items(), key=lambda kv: (hb.modality_sort_key(kv[0][0]), kv[0][1], kv[0][2])):
        summary_rows.append({"modality": key[0], "tool": key[1], "status": key[2], "n_samples": count})

    meta = {
        "analysis_tools_by_modality": expected,
        "available_tools_by_modality": {k: sorted(v) for k, v in available_tools.items()},
        "not_available_tools_by_modality": {k: sorted(v) for k, v in not_available_tools.items()},
    }
    return sample_rows, summary_rows, meta


def enrich_metadata(output_dir, config, cohort_rows, sequencing_rows, supported_loci, tool_meta):
    metadata_path = Path(output_dir) / "tables" / "benchmark_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.exists() else {}
    metadata["truth_acquisition_date"] = string_value(config.get("truth", {}).get("acquisition_date", ""))
    metadata["final_tri_modal_cohort_size"] = sum(1 for row in cohort_rows if row.get("include") == "1")
    metadata["supported_loci"] = supported_loci
    metadata["population_counts"] = build_population_counts(cohort_rows)
    metadata["excluded_sample_counts_by_reason"] = build_excluded_summary(cohort_rows)
    metadata["per_modality_sample_counts"] = build_modality_counts(sequencing_rows, cohort_rows)
    metadata["tool_coverage_policy"] = config.get("benchmark", {}).get("tool_coverage_policy", "phase_gated")
    metadata.update(tool_meta)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")


def update_reference_metadata(output_dir, config, supported_loci):
    path = Path(output_dir) / "tables" / "reference_metadata.tsv"
    rows = load_rows(path) if path.exists() else [{}]
    row = rows[0] if rows else {}
    row["truth_acquisition_date"] = string_value(config.get("truth", {}).get("acquisition_date", ""))
    row["supported_loci"] = ",".join(supported_loci)
    row["benchmark_mode"] = string_value(config.get("benchmark", {}).get("mode", ""))
    write_tsv(path, [row], ["truth_source", "truth_path", "imgt_hla_version", "primary_resolution", "secondary_resolutions", "truth_acquisition_date", "supported_loci", "benchmark_mode"])


def split_sample_sets(cohort_rows):
    grouped = defaultdict(list)
    for row in cohort_rows:
        if row.get("include") != "1":
            continue
        split = hb.clean_token(row.get("split", "")).lower()
        if split:
            grouped[split].append(row["sample"])
    required = {"training", "validation", "holdout"}
    if not required.issubset(grouped):
        return None
    return {key: sorted(values) for key, values in grouped.items()}


def tune_validation_support(validation_dir, config):
    path = Path(validation_dir) / "tables" / "abstention_tradeoff.tsv"
    if not path.exists():
        return config.get("benchmark", {}).get("consensus", {}).get("min_support", 0.55)
    rows = load_rows(path)
    if not rows:
        return config.get("benchmark", {}).get("consensus", {}).get("min_support", 0.55)
    best = max(
        rows,
        key=lambda row: (
            float(row.get("overall_correct_call_rate", 0.0) or 0.0),
            float(row.get("accuracy_among_called", 0.0) or 0.0),
            float(row.get("call_rate", 0.0) or 0.0),
        ),
    )
    return float(best.get("min_support", 0.55) or 0.55)


def copy_training_artifacts(training_dir, output_dir):
    source_tables = Path(training_dir) / "tables"
    dest_tables = Path(output_dir) / "tables"
    for name in [
        "tool_confidence_weights.tsv",
        "tool_confidence_weights_by_gene.tsv",
        "tool_confidence_weights_cv.tsv",
        "cross_validation_weight_summary.tsv",
        "confidence_calibration_summary.tsv",
        "consensus_runtime_weights.json",
        "tool_confidence_weights_by_population.tsv",
        "tool_population_diagnostics.tsv",
        "tool_population_calibration.tsv",
    ]:
        src = source_tables / name
        if src.exists():
            shutil.copy2(src, dest_tables / name)


def main():
    args = parse_args()
    config = load_yaml(Path(args.config))
    if args.benchmark_mode:
        config.setdefault("benchmark", {})
        config["benchmark"]["mode"] = args.benchmark_mode
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    truth_manifest, sequencing_manifest, cohort_manifest = ensure_manifests(config, output_dir)
    sequencing_rows = load_rows(sequencing_manifest)
    cohort_rows = load_rows(cohort_manifest)
    supported_loci = parse_truth_supported_loci(config, cohort_rows)
    included_rows = [row for row in cohort_rows if row.get("include") == "1"]
    if not included_rows:
        raise SystemExit("No tri-modal truth-backed samples were included in the cohort manifest.")

    truth_map = hb.load_truth(config["truth"])
    all_samples = [row["sample"] for row in included_rows]
    population_manifest_path = str(cohort_manifest)
    split_sets = split_sample_sets(cohort_rows)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        if split_sets:
            train_truth = tmpdir / "training_truth.tsv"
            validation_truth = tmpdir / "validation_truth.tsv"
            holdout_truth = tmpdir / "holdout_truth.tsv"
            write_filtered_truth(train_truth, build_filtered_truth_rows(truth_map, split_sets["training"], supported_loci))
            write_filtered_truth(validation_truth, build_filtered_truth_rows(truth_map, split_sets["validation"], supported_loci))
            write_filtered_truth(holdout_truth, build_filtered_truth_rows(truth_map, split_sets["holdout"], supported_loci))

            train_cfg = subset_config(config, train_truth, supported_loci)
            training_dir = output_dir / "training"
            validation_dir = output_dir / "validation"
            run_benchmark(
                train_cfg,
                training_dir,
                weight_alpha=args.weight_alpha,
                weight_beta=args.weight_beta,
                benchmark_mode=config.get("benchmark", {}).get("mode", ""),
                population_manifest=population_manifest_path,
            )

            validation_cfg = subset_config(config, validation_truth, supported_loci)
            validation_cfg.setdefault("benchmark", {})
            validation_cfg["benchmark"]["runtime_weight_override"] = str(training_dir / "tables" / "consensus_runtime_weights.json")
            run_benchmark(
                validation_cfg,
                validation_dir,
                weight_alpha=args.weight_alpha,
                weight_beta=args.weight_beta,
                benchmark_mode=config.get("benchmark", {}).get("mode", ""),
                population_manifest=population_manifest_path,
            )

            tuned_support = tune_validation_support(validation_dir, config)
            holdout_cfg = subset_config(config, holdout_truth, supported_loci)
            holdout_cfg.setdefault("benchmark", {})
            holdout_cfg["benchmark"]["runtime_weight_override"] = str(training_dir / "tables" / "consensus_runtime_weights.json")
            holdout_cfg.setdefault("benchmark", {}).setdefault("consensus", {})
            holdout_cfg["benchmark"]["consensus"]["min_support"] = tuned_support
            run_benchmark(
                holdout_cfg,
                output_dir,
                weight_alpha=args.weight_alpha,
                weight_beta=args.weight_beta,
                benchmark_mode=config.get("benchmark", {}).get("mode", ""),
                population_manifest=population_manifest_path,
            )
            copy_training_artifacts(training_dir, output_dir)
        else:
            all_truth = tmpdir / "all_truth.tsv"
            write_filtered_truth(all_truth, build_filtered_truth_rows(truth_map, all_samples, supported_loci))
            all_cfg = subset_config(config, all_truth, supported_loci)
            run_benchmark(
                all_cfg,
                output_dir,
                weight_alpha=args.weight_alpha,
                weight_beta=args.weight_beta,
                benchmark_mode=config.get("benchmark", {}).get("mode", ""),
                population_manifest=population_manifest_path,
            )

    copy_manifest_outputs(output_dir, truth_manifest, sequencing_manifest, cohort_manifest)

    harmonized_path = output_dir / "tables" / "harmonized_benchmark_rows.tsv"
    harmonized_rows = load_rows(harmonized_path) if harmonized_path.exists() else []
    included_samples = {row["sample"] for row in included_rows}
    sample_tool_rows, summary_tool_rows, tool_meta = build_tool_availability(config, harmonized_rows, included_samples)
    write_tsv(output_dir / "tables" / "tool_availability_by_sample.tsv", sample_tool_rows, ["sample", "modality", "tool", "status", "observed_rows"])
    write_tsv(output_dir / "tables" / "tool_availability_by_modality.tsv", summary_tool_rows, ["modality", "tool", "status", "n_samples"])

    enrich_metadata(output_dir, config, cohort_rows, sequencing_rows, supported_loci, tool_meta)
    if split_sets:
        metadata_path = output_dir / "tables" / "benchmark_metadata.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.exists() else {}
        metadata["split_counts"] = {name: len(samples) for name, samples in split_sets.items()}
        metadata["weight_learning_split"] = "training"
        metadata["threshold_tuning_split"] = "validation"
        metadata["final_evaluation_split"] = "holdout"
        metadata["population_manifest_path"] = population_manifest_path
        metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    update_reference_metadata(output_dir, config, supported_loci)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
