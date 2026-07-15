#!/usr/bin/env python3
"""Benchmark HLA typing runs, cohort-aware confidence weights, and figure-ready outputs."""

import argparse
import copy
import csv
import glob
import json
import math
import re
import statistics
import sys
from collections import defaultdict, namedtuple
from pathlib import Path

import yaml

GENE_ORDER = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5"]
MODALITY_ORDER = ["wes", "wgs", "rnaseq"]
MODALITY_STYLES = {
    "wes": {"label": "WES", "fill": "#1F6FEB"},
    "wgs": {"label": "WGS", "fill": "#F97316"},
    "rnaseq": {"label": "RNA-seq", "fill": "#0F766E"},
}
MISSING_TOKENS = {"", "na", "n/a", "none", "null", "nan", "-", ".", "failed", "fail", "no_call", "no result"}
DEFAULT_TARGET_READS = 50.0
WEIGHT_VERSION = "read_confidence_v2"
RNA_MIN_READ_SUPPORT = 10.0  # reads; below this threshold RNA calls are treated as unreliable
DEFAULT_BENCHMARK_MODE = "legacy_heuristic"
DEFAULT_PROBABILISTIC_CALIBRATION = {
    "method": "platt",
    "cv_strategy": "loo",
    "min_rows_for_isotonic": 20,
}
DEFAULT_BAYESIAN_SHRINKAGE = {
    "alpha": 1.0,
    "beta": 1.0,
}
DEFAULT_CONFIDENCE_GUARDRAIL = {
    "enabled": True,
    "max_expected_calibration_error_for_boost": 0.35,
    "max_brier_score_for_boost": 0.35,
    "min_confidence_coverage_for_boost": 0.5,
    "guardrail_steepness": 0,  # 0 = strict binary cutoff; >0 = soft exponential dampening (e.g. 2.0)
}
GROUP_SUFFIXES = {"G", "P"}


CallRecord = namedtuple("CallRecord", ["sample", "gene", "allele1", "allele2", "source_file"])
RunRecord = namedtuple(
    "RunRecord",
    [
        "tool",
        "modality",
        "result_glob",
        "parser",
        "runtime_glob",
        "sample_from_path_regex",
        "sample_from_path_group",
        "gene_filter",
        "confidence_glob",
        "confidence_parser",
        "confidence_sample_from_path_regex",
        "confidence_sample_from_path_group",
        "confidence_defaults",
    ],
)


def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark HLA typing runs.")
    parser.add_argument("--config", required=True, help="YAML configuration describing truth and run outputs.")
    parser.add_argument("--output-dir", required=True, help="Directory for harmonized tables, summaries, and SVG figures.")
    parser.add_argument("--weight-alpha", type=float, default=0.7,
                        help="Weight on base_reliability in final_weight formula (default: 0.7).")
    parser.add_argument("--weight-beta", type=float, default=0.3,
                        help="Weight on effective_confidence in final_weight formula (default: 0.3).")
    parser.add_argument("--cv-folds", type=int, default=0,
                        help="Number of cross-validation folds for weight uncertainty estimation (0 = disabled).")
    parser.add_argument("--population-manifest", default=None,
                        help="TSV with columns sample/superpopulation/population for ancestry-stratified weights.")
    parser.add_argument("--benchmark-mode", default="",
                        help="Override benchmark.mode from config (legacy_heuristic, probabilistic_recalibrated, bayesian_shrinkage).")
    return parser.parse_args()


def load_yaml(path):
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def detect_delimiter(path):
    return "," if path.suffix.lower() == ".csv" else "\t"


def clean_token(value):
    if value is None:
        return ""
    return str(value).strip()


def normalize_column_name(name):
    return re.sub(r"[^a-z0-9]+", "", clean_token(name).lower())


def first_existing(mapping, candidates):
    for candidate in candidates:
        if candidate in mapping:
            return candidate
    return None


def normalize_gene(raw):
    token = clean_token(raw).upper().replace("HLA_", "").replace("HLA-", "").replace("HLA", "")
    token = token.strip("_* -")
    if token.startswith("DRB"):
        match = re.match(r"(DRB\d)", token)
        return match.group(1) if match else token
    return token


def gene_sort_key(gene):
    try:
        return (GENE_ORDER.index(gene), gene)
    except ValueError:
        return (len(GENE_ORDER), gene)


def modality_sort_key(modality):
    try:
        return (MODALITY_ORDER.index(modality), modality)
    except ValueError:
        return (len(MODALITY_ORDER), modality)


def ordered_modalities(rows):
    return sorted({row["modality"] for row in rows}, key=modality_sort_key)


def modality_style(modality):
    return MODALITY_STYLES.get(modality, {"label": modality.upper(), "fill": "#6B7280"})


def benchmark_mode(config, override=""):
    token = clean_token(override or (config or {}).get("benchmark", {}).get("mode", "")).lower()
    if token in {"legacy_heuristic", "probabilistic_recalibrated", "bayesian_shrinkage"}:
        return token
    return DEFAULT_BENCHMARK_MODE


def probabilistic_calibration_settings(config):
    benchmark_cfg = (config or {}).get("benchmark", {}) if isinstance(config, dict) else {}
    raw = benchmark_cfg.get("probabilistic_calibration", {}) if isinstance(benchmark_cfg, dict) else {}
    settings = dict(DEFAULT_PROBABILISTIC_CALIBRATION)
    if isinstance(raw, dict):
        settings.update(raw)
    settings["method"] = clean_token(settings.get("method", DEFAULT_PROBABILISTIC_CALIBRATION["method"])).lower() or "platt"
    settings["cv_strategy"] = clean_token(settings.get("cv_strategy", DEFAULT_PROBABILISTIC_CALIBRATION["cv_strategy"])).lower() or "loo"
    settings["min_rows_for_isotonic"] = int(settings.get("min_rows_for_isotonic", DEFAULT_PROBABILISTIC_CALIBRATION["min_rows_for_isotonic"]))
    return settings


def bayesian_shrinkage_settings(config):
    benchmark_cfg = (config or {}).get("benchmark", {}) if isinstance(config, dict) else {}
    raw = benchmark_cfg.get("bayesian_shrinkage", {}) if isinstance(benchmark_cfg, dict) else {}
    settings = dict(DEFAULT_BAYESIAN_SHRINKAGE)
    if isinstance(raw, dict):
        settings.update(raw)
    settings["alpha"] = float(settings.get("alpha", DEFAULT_BAYESIAN_SHRINKAGE["alpha"]))
    settings["beta"] = float(settings.get("beta", DEFAULT_BAYESIAN_SHRINKAGE["beta"]))
    return settings


def benchmark_analysis_settings(config):
    benchmark_cfg = (config or {}).get("benchmark_analysis", {}) if isinstance(config, dict) else {}
    raw = benchmark_cfg.get("threshold_sweeps", {}) if isinstance(benchmark_cfg, dict) else {}
    if not isinstance(raw, dict):
        return {"enabled": False, "weighted_consensus": {}, "champion_challenger": {}}
    return {
        "enabled": bool(raw.get("enabled")),
        "weighted_consensus": raw.get("weighted_consensus", {}) if isinstance(raw.get("weighted_consensus", {}), dict) else {},
        "champion_challenger": raw.get("champion_challenger", {}) if isinstance(raw.get("champion_challenger", {}), dict) else {},
    }


def weighted_threshold_sweep_settings(config):
    analysis = benchmark_analysis_settings(config)
    raw = analysis.get("weighted_consensus", {})
    support_values = [float(value) for value in raw.get("min_support_values", []) or []]
    margin_values = [float(value) for value in raw.get("min_margin_values", []) or []]
    return {
        "enabled": analysis.get("enabled", False) and bool(support_values) and bool(margin_values),
        "min_support_values": support_values,
        "min_margin_values": margin_values,
    }


def champion_override_sweep_settings(config):
    analysis = benchmark_analysis_settings(config)
    raw = analysis.get("champion_challenger", {})
    return {
        "enabled": analysis.get("enabled", False),
        "champion_by_gene": {normalize_gene(gene): clean_token(tool) for gene, tool in (raw.get("champion_by_gene", {}) or {}).items() if clean_token(tool)},
        "min_challenger_support_fraction_values": [float(value) for value in raw.get("min_challenger_support_fraction_values", []) or []],
        "min_challenger_margin_values": [float(value) for value in raw.get("min_challenger_margin_values", []) or []],
        "min_supporting_tools_values": [int(value) for value in raw.get("min_supporting_tools_values", []) or []],
        "require_non_ambiguity_override": bool(raw.get("require_non_ambiguity_override", True)),
    }


def sweep_target_modality(rows):
    modalities = sorted({clean_token(row.get("modality")).lower() for row in rows if clean_token(row.get("modality"))}, key=modality_sort_key)
    if not modalities:
        return ""
    if len(modalities) == 1:
        return modalities[0]
    for preferred in ["wes", "rnaseq", "wgs"]:
        if preferred in modalities:
            return preferred
    return modalities[0]


def sweep_file_prefix(modality):
    token = clean_token(modality).lower()
    if token == "rnaseq":
        return "rna"
    return token


def spechla_forced_accuracy(config, modality):
    token = clean_token(modality).lower()
    if token != "rnaseq":
        return False
    for run in (config or {}).get("runs", []):
        if clean_token(run.get("tool")) == "SpecHLA" and clean_token(run.get("modality")).lower() == token:
            return not bool(run.get("coverage_only", False))
    return False


def is_missing(value):
    return clean_token(value).lower() in MISSING_TOKENS


def split_allele_components(raw):
    value = clean_token(raw)
    if is_missing(value):
        return "", [], ""
    value = value.split(";")[0].strip()
    value = value.replace("HLA-", "").replace("HLA_", "").replace("_", "*", 1)
    value = value.replace("_", ":")
    value = value.replace(" ", "")
    if "*" not in value:
        return value.upper(), [], ""
    gene, fields = value.split("*", 1)
    suffix = ""
    match = re.search(r"([A-Za-z])$", fields)
    if match:
        suffix = match.group(1).upper()
        fields = fields[:-1]
    parts = [part for part in fields.split(":") if part]
    return gene.upper(), parts, suffix


def normalize_allele(raw, resolution=2, preserve_group_suffix=False):
    gene, parts, suffix = split_allele_components(raw)
    if not gene:
        return ""
    if not parts:
        return gene if "*" not in clean_token(raw) else ""
    normalized = "%s*%s" % (gene, ":".join(parts[:resolution]))
    if preserve_group_suffix and suffix in GROUP_SUFFIXES:
        return normalized + suffix
    return normalized


def sort_alleles(alleles, resolution=2, preserve_group_suffix=False):
    normalized = []
    for allele in alleles:
        value = normalize_allele(allele, resolution=resolution, preserve_group_suffix=preserve_group_suffix)
        if value:
            normalized.append(value)
    normalized = sorted(normalized)
    if len(normalized) == 0:
        return ("", "")
    if len(normalized) == 1:
        return (normalized[0], normalized[0])
    return (normalized[0], normalized[1])


def normalize_group_allele(raw, suffix, resolution=3):
    normalized = normalize_allele(raw, resolution=resolution, preserve_group_suffix=True)
    return normalized if normalized.endswith(suffix) else ""


def sort_group_alleles(alleles, suffix, resolution=3):
    normalized = sorted([normalize_group_allele(allele, suffix, resolution=resolution) for allele in alleles if normalize_group_allele(allele, suffix, resolution=resolution)])
    if len(normalized) == 0:
        return ("", "")
    if len(normalized) == 1:
        return (normalized[0], normalized[0])
    return (normalized[0], normalized[1])


def match_pair(pair_a, pair_b):
    return all(pair_a) and all(pair_b) and pair_a == pair_b


def ambiguity_match_value(is_comparable, is_match):
    if not is_comparable:
        return ""
    return "1" if is_match else "0"


def allele_match_grade(raw_typed, raw_truth):
    typed_2field = sort_alleles(raw_typed, resolution=2)
    truth_2field = sort_alleles(raw_truth, resolution=2)
    typed_3field = sort_alleles(raw_typed, resolution=3)
    truth_3field = sort_alleles(raw_truth, resolution=3)
    typed_g_group = sort_group_alleles(raw_typed, "G", resolution=3)
    truth_g_group = sort_group_alleles(raw_truth, "G", resolution=3)
    typed_p_group = sort_group_alleles(raw_typed, "P", resolution=3)
    truth_p_group = sort_group_alleles(raw_truth, "P", resolution=3)
    if match_pair(typed_3field, truth_3field):
        return "exact_3field"
    if match_pair(typed_2field, truth_2field):
        return "exact_2field"
    if all(typed_g_group) and all(truth_g_group) and typed_g_group == truth_g_group:
        return "g_group"
    if all(typed_p_group) and all(truth_p_group) and typed_p_group == truth_p_group:
        return "p_group"
    return "incompatible"


def ambiguity_compatible(raw_typed, raw_truth):
    typed_tokens = [clean_token(value) for value in raw_typed if clean_token(value)]
    truth_tokens = [clean_token(value) for value in raw_truth if clean_token(value)]
    if len(typed_tokens) != 2 or len(truth_tokens) != 2:
        return False
    truth_2field = {normalize_allele(value, resolution=2) for value in truth_tokens if normalize_allele(value, resolution=2)}
    if not truth_2field:
        return False
    for token in typed_tokens:
        options = [part.strip() for part in token.split("/") if part.strip()]
        if not options:
            return False
        option_matches = {normalize_allele(option, resolution=2) for option in options if normalize_allele(option, resolution=2)}
        if not (option_matches & truth_2field):
            return False
    return True


def evaluate_ambiguity(raw_typed, raw_truth, configured_resolution):
    typed_configured = sort_alleles(raw_typed, resolution=configured_resolution)
    truth_configured = sort_alleles(raw_truth, resolution=configured_resolution)
    typed_2field = sort_alleles(raw_typed, resolution=2)
    truth_2field = sort_alleles(raw_truth, resolution=2)
    typed_3field = sort_alleles(raw_typed, resolution=3)
    truth_3field = sort_alleles(raw_truth, resolution=3)
    typed_g_group = sort_group_alleles(raw_typed, "G", resolution=3)
    truth_g_group = sort_group_alleles(raw_truth, "G", resolution=3)
    typed_p_group = sort_group_alleles(raw_typed, "P", resolution=3)
    truth_p_group = sort_group_alleles(raw_truth, "P", resolution=3)

    configured_match = match_pair(typed_configured, truth_configured)
    match_2field = match_pair(typed_2field, truth_2field)
    match_3field = match_pair(typed_3field, truth_3field)
    g_group_comparable = all(typed_g_group) and all(truth_g_group)
    p_group_comparable = all(typed_p_group) and all(truth_p_group)
    match_g_group = g_group_comparable and typed_g_group == truth_g_group
    match_p_group = p_group_comparable and typed_p_group == truth_p_group

    ambiguity_compatible_match = ambiguity_compatible(raw_typed, raw_truth)
    if not all(typed_configured):
        match_grade = "missing"
    elif match_3field:
        match_grade = "exact_3field"
    elif match_2field:
        match_grade = "exact_2field"
    elif match_g_group:
        match_grade = "g_group"
    elif match_p_group:
        match_grade = "p_group"
    elif ambiguity_compatible_match:
        match_grade = "ambiguity_compatible"
    else:
        match_grade = "mismatch"
    compatibility_grade = (
        "missing" if not all(typed_configured) else
        "exact_3field" if match_3field else
        "exact_2field" if match_2field else
        "g_group" if match_g_group else
        "p_group" if match_p_group else
        "ambiguity_compatible" if ambiguity_compatible_match else
        "incompatible"
    )

    return {
        "typed_configured": typed_configured,
        "truth_configured": truth_configured,
        "typed_2field": typed_2field,
        "truth_2field": truth_2field,
        "typed_3field": typed_3field,
        "truth_3field": truth_3field,
        "configured_match": configured_match,
        "match_2field": match_2field,
        "match_3field": match_3field,
        "g_group_comparable": g_group_comparable,
        "p_group_comparable": p_group_comparable,
        "match_g_group": match_g_group,
        "match_p_group": match_p_group,
        "match_grade": match_grade,
        "ambiguity_compatible": ambiguity_compatible_match,
        "is_resolution_compatible": match_2field or match_3field or match_g_group or match_p_group or ambiguity_compatible_match,
        "compatibility_grade": compatibility_grade,
    }


def coerce_float(value):
    token = clean_token(value)
    if not token:
        return None
    try:
        return float(token)
    except ValueError:
        match = re.search(r"[-+]?\d*\.?\d+", token)
        return float(match.group(0)) if match else None


def as_string_number(value):
    if value is None or value == "":
        return ""
    return "%.4f" % float(value)


def clip_unit(value):
    if value is None:
        return None
    return max(0.0, min(1.0, float(value)))


def read_table(path, delimiter=None):
    delim = delimiter or detect_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        filtered = [line for line in handle if line.strip() and not line.startswith("#")]
    if not filtered:
        return []
    reader = csv.DictReader(filtered, delimiter=delim)
    return [{clean_token(k): clean_token(v) for k, v in row.items() if k is not None} for row in reader]


def read_raw_lines(path):
    with path.open("r", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


def derive_sample_from_path(path, regex, group):
    if regex:
        match = re.search(regex, str(path))
        if match:
            return clean_token(match.group(group))
    return clean_token(path.parent.name or path.stem)


def extract_wide_gene_columns(columns):
    column_map = {}
    for col in columns:
        match = re.match(r"(?:HLA[_-])?([A-Za-z0-9]+)[_-]([12])$", clean_token(col))
        if not match:
            continue
        gene = normalize_gene(match.group(1))
        column_map.setdefault(gene, {})[match.group(2)] = col
    return dict((gene, (parts["1"], parts["2"])) for gene, parts in column_map.items() if "1" in parts and "2" in parts)


def extract_wide_confidence_columns(columns):
    column_map = defaultdict(dict)
    for col in columns:
        token = clean_token(col)
        match = re.match(r"(?:HLA[_-])?([A-Za-z0-9]+)[_-](confidence(?:[_-]?score)?|read(?:[_-]?support|s)?|score)$", token, re.IGNORECASE)
        if not match:
            continue
        gene = normalize_gene(match.group(1))
        measure = normalize_column_name(match.group(2))
        if measure in {"confidence", "confidencescore", "score"}:
            column_map[gene]["confidence_score"] = col
        elif measure in {"readsupport", "reads"}:
            column_map[gene]["read_support"] = col
    return dict(column_map)


def long_truth_rows(rows):
    truth = defaultdict(dict)
    if not rows:
        return truth
    key_map = dict((normalize_column_name(name), name) for name in rows[0].keys())
    sample_key = first_existing(key_map, ["sample", "sampleid", "id"])
    gene_key = first_existing(key_map, ["gene", "locus"])
    allele1_key = first_existing(key_map, ["allele1", "truthallele1", "typedallele1", "ggroup1"])
    allele2_key = first_existing(key_map, ["allele2", "truthallele2", "typedallele2", "ggroup2"])
    if not (sample_key and gene_key and allele1_key and allele2_key):
        return truth
    for row in rows:
        sample = clean_token(row[key_map[sample_key]])
        gene = normalize_gene(row[key_map[gene_key]])
        truth[sample][gene] = sort_alleles([row[key_map[allele1_key]], row[key_map[allele2_key]]])
    return truth


def wide_truth_rows(rows):
    truth = defaultdict(dict)
    if not rows:
        return truth
    sample_key = next((col for col in rows[0].keys() if normalize_column_name(col) in {"sample", "sampleid"}), None)
    if not sample_key:
        return truth
    gene_pairs = extract_wide_gene_columns(rows[0].keys())
    for row in rows:
        sample = clean_token(row.get(sample_key, ""))
        if not sample:
            continue
        for gene, pair in gene_pairs.items():
            truth[sample][gene] = sort_alleles([row.get(pair[0], ""), row.get(pair[1], "")])
    return truth


def load_truth(truth_cfg):
    rows = read_table(Path(truth_cfg["path"]), truth_cfg.get("delimiter"))
    truth = long_truth_rows(rows)
    if truth:
        return truth
    truth = wide_truth_rows(rows)
    if truth:
        return truth
    raise ValueError("Could not parse truth table from %s" % truth_cfg["path"])


def parse_gene_list(values):
    if not values:
        return []
    if isinstance(values, str):
        tokens = re.split(r"[,;|\s]+", values)
    else:
        tokens = []
        for value in values:
            tokens.extend(parse_gene_list(value) if isinstance(value, (list, tuple, set)) else [value])
    genes = []
    for token in tokens:
        gene = normalize_gene(token)
        if gene:
            genes.append(gene)
    ordered = []
    for gene in genes:
        if gene not in ordered:
            ordered.append(gene)
    return ordered


def resolve_benchmark_genes(config):
    truth_cfg = config.get("truth", {})
    benchmark_cfg = config.get("benchmark", {})
    genes = parse_gene_list(truth_cfg.get("supported_loci") or benchmark_cfg.get("genes"))
    return sorted(genes, key=gene_sort_key)


def parse_wide_hla_file(path, sample_override=""):
    rows = read_table(path)
    if not rows:
        return []
    header_pairs = extract_wide_gene_columns(rows[0].keys())
    sample_key = next((col for col in rows[0].keys() if normalize_column_name(col) in {"sample", "sampleid"}), None)
    records = []
    for row in rows:
        sample = sample_override or (clean_token(row.get(sample_key, "")) if sample_key else "")
        for gene, pair in header_pairs.items():
            records.append(CallRecord(sample=sample, gene=gene, allele1=row.get(pair[0], ""), allele2=row.get(pair[1], ""), source_file=str(path)))
    return records


def parse_long_hla_file(path, sample_override=""):
    rows = read_table(path)
    if not rows:
        return []
    key_map = dict((normalize_column_name(name), name) for name in rows[0].keys())
    sample_key = first_existing(key_map, ["sample", "sampleid", "id"])
    gene_key = first_existing(key_map, ["gene", "locus"])
    allele1_key = first_existing(key_map, ["allele1", "typedallele1", "truthallele1", "ggroup1"])
    allele2_key = first_existing(key_map, ["allele2", "typedallele2", "truthallele2", "ggroup2"])
    if not (gene_key and allele1_key and allele2_key):
        return []
    records = []
    for row in rows:
        sample = clean_token(row.get(key_map[sample_key], "")) if sample_key else ""
        sample = sample or sample_override
        records.append(CallRecord(
            sample=sample,
            gene=normalize_gene(row[key_map[gene_key]]),
            allele1=row[key_map[allele1_key]],
            allele2=row[key_map[allele2_key]],
            source_file=str(path),
        ))
    return records


def parse_hlahd_file(path, sample_override):
    sample = sample_override or clean_token(path.stem)
    records = []
    for line in read_raw_lines(path):
        if not line.strip() or line.startswith('#'):
            continue
        fields = [clean_token(part) for part in line.split('	')]
        if len(fields) < 3:
            continue
        gene = normalize_gene(fields[0])
        if not gene:
            continue
        records.append(CallRecord(sample=sample, gene=gene, allele1=fields[1], allele2=fields[2], source_file=str(path)))
    return records


def parse_optitype_file(path, sample_override):
    rows = read_table(path, "\t")
    if not rows:
        return []
    row = rows[0]
    sample = sample_override or clean_token(row.get("Sample", "")) or path.stem
    out = []
    for gene, cols in {"A": ("A1", "A2"), "B": ("B1", "B2"), "C": ("C1", "C2")}.items():
        if cols[0] in row and cols[1] in row:
            out.append(CallRecord(sample=sample, gene=gene, allele1=row.get(cols[0], ""), allele2=row.get(cols[1], ""), source_file=str(path)))
    return out


def parse_arcashla_json(path, sample_override):
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    sample = sample_override or clean_token(payload.get("sample", "")) or path.stem
    genotype = payload.get("genotype") or payload.get("genes") or payload
    out = []
    if isinstance(genotype, dict):
        for raw_gene, alleles in genotype.items():
            if not isinstance(alleles, list):
                continue
            allele1 = alleles[0] if len(alleles) > 0 else ""
            allele2 = alleles[1] if len(alleles) > 1 else allele1
            out.append(CallRecord(sample=sample, gene=normalize_gene(raw_gene), allele1=allele1, allele2=allele2, source_file=str(path)))
    return out


def parse_auto_file(path, sample_override):
    if path.suffix.lower() == ".json":
        return parse_arcashla_json(path, sample_override)
    rows = read_table(path)
    if rows:
        header_names = {normalize_column_name(col) for col in rows[0].keys()}
        if {"a1", "a2", "b1", "b2", "c1", "c2"} & header_names:
            return parse_optitype_file(path, sample_override)
        if any(name in header_names for name in {"gene", "locus"}) and any(name in header_names for name in {"allele1", "typedallele1", "truthallele1", "ggroup1"}) and any(name in header_names for name in {"allele2", "typedallele2", "truthallele2", "ggroup2"}):
            return parse_long_hla_file(path, sample_override)
        if extract_wide_gene_columns(rows[0].keys()):
            return parse_wide_hla_file(path, sample_override)
    return []


def parse_result_file(path, parser_name, sample_override):
    parser_name = parser_name.lower()
    if parser_name in {"auto", ""}:
        return parse_auto_file(path, sample_override)
    if parser_name in {"spec_hla_result", "wide_hla_table"}:
        return parse_wide_hla_file(path, sample_override)
    if parser_name == "hlahd_table":
        return parse_hlahd_file(path, sample_override)
    if parser_name in {"long_hla_table", "polysolver_table", "kourami_table", "t1k_table", "seq2hla_table", "locityper_table"}:
        return parse_auto_file(path, sample_override)
    if parser_name == "optitype_tsv":
        return parse_optitype_file(path, sample_override)
    if parser_name == "arcashla_json":
        return parse_arcashla_json(path, sample_override)
    raise ValueError("Unsupported parser: %s" % parser_name)


def parse_time_log(path):
    user_time = None
    sys_time_val = None
    memory_kb = None
    for line in read_raw_lines(path):
        if "User time (seconds)" in line:
            user_time = coerce_float(line)
        elif "System time (seconds)" in line:
            sys_time_val = coerce_float(line)
        elif "Maximum resident set size (kbytes)" in line:
            memory_kb = coerce_float(line)
    runtime_hours = round(((user_time or 0.0) + (sys_time_val or 0.0)) / 3600.0, 4) if (user_time is not None or sys_time_val is not None) else None
    memory_gb = round((memory_kb or 0.0) / 1000000.0, 4) if memory_kb is not None else None
    return runtime_hours, memory_gb


def summarize_resource_logs(paths):
    runtimes = []
    memories = []
    for path in paths:
        runtime_hours, memory_gb = parse_time_log(path)
        if runtime_hours is not None:
            runtimes.append(runtime_hours)
        if memory_gb is not None:
            memories.append(memory_gb)
    return (round(statistics.median(runtimes), 4) if runtimes else None, round(statistics.median(memories), 4) if memories else None)


def normalize_confidence_values(raw_score, read_support, defaults):
    raw_score_val = coerce_float(raw_score)
    read_support_val = coerce_float(read_support)
    target_reads = coerce_float((defaults or {}).get("target_reads")) or DEFAULT_TARGET_READS
    if raw_score_val is not None:
        return clip_unit(raw_score_val), read_support_val, raw_score_val, "native_score", raw_score_val
    if read_support_val is not None:
        return clip_unit(read_support_val / float(target_reads)), read_support_val, None, "read_support", read_support_val
    return None, None, raw_score_val, "missing", None


def build_confidence_entry(sample, gene, raw_score, read_support, defaults, source, source_file):
    confidence_score, read_support_val, raw_confidence_val, raw_score_family, raw_score_value = normalize_confidence_values(raw_score, read_support, defaults)
    return {
        "sample": sample,
        "gene": gene,
        "confidence_score": confidence_score,
        "read_support": read_support_val,
        "raw_confidence": raw_confidence_val,
        "raw_score_family": raw_score_family,
        "raw_score_value": raw_score_value,
        "confidence_source": source if confidence_score is not None or read_support_val is not None or raw_confidence_val is not None else "missing",
        "source_file": str(source_file),
    }


def parse_optitype_objective_confidence(path, sample_override, defaults):
    rows = read_table(path, "\t")
    if not rows:
        return []
    row = rows[0]
    sample = sample_override or clean_token(row.get("Sample", "")) or path.stem
    objective = row.get("Objective", "")
    return [build_confidence_entry(sample, gene, objective, None, defaults, "optitype_result_objective", path) for gene in ["A", "B", "C"]]


def parse_long_confidence_table(path, sample_override, defaults):
    rows = read_table(path)
    if not rows:
        return []
    key_map = dict((normalize_column_name(name), name) for name in rows[0].keys())
    sample_key = first_existing(key_map, ["sample", "sampleid", "id"])
    gene_key = first_existing(key_map, ["gene", "locus"])
    score_key = first_existing(key_map, ["confidencescore", "confidence", "score"])
    reads_key = first_existing(key_map, ["readsupport", "reads"])
    if not (sample_key and gene_key and (score_key or reads_key)):
        return []
    out = []
    for row in rows:
        sample = clean_token(row[key_map[sample_key]]) or sample_override
        gene = normalize_gene(row[key_map[gene_key]])
        out.append(build_confidence_entry(sample, gene, row.get(key_map[score_key], "") if score_key else None, row.get(key_map[reads_key], "") if reads_key else None, defaults, "long_confidence_table", path))
    return out


def parse_wide_confidence_table(path, sample_override, defaults):
    rows = read_table(path)
    if not rows:
        return []
    sample_key = next((col for col in rows[0].keys() if normalize_column_name(col) in {"sample", "sampleid"}), None)
    wide_map = extract_wide_confidence_columns(rows[0].keys())
    out = []
    for row in rows:
        sample = clean_token(row.get(sample_key, "")) if sample_key else sample_override
        for gene, measures in wide_map.items():
            out.append(build_confidence_entry(sample, gene, row.get(measures.get("confidence_score", ""), ""), row.get(measures.get("read_support", ""), ""), defaults, "wide_confidence_table", path))
    return out


def parse_json_confidence_table(path, sample_override, defaults):
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    out = []
    if isinstance(payload, list):
        records = payload
    else:
        records = payload.get("records") or payload.get("rows") or []
        if not records and isinstance(payload.get("genes"), dict):
            sample = clean_token(payload.get("sample", "")) or sample_override or path.stem
            for raw_gene, values in payload["genes"].items():
                if isinstance(values, dict):
                    out.append(build_confidence_entry(sample, normalize_gene(raw_gene), values.get("confidence_score"), values.get("read_support"), defaults, "json_confidence_table", path))
            return out
    for record in records:
        if not isinstance(record, dict):
            continue
        sample = clean_token(record.get("sample", "")) or sample_override or clean_token(payload.get("sample", "")) or path.stem
        gene = normalize_gene(record.get("gene", ""))
        if not gene:
            continue
        out.append(build_confidence_entry(sample, gene, record.get("confidence_score", record.get("confidence")), record.get("read_support", record.get("reads")), defaults, "json_confidence_table", path))
    return out


def parse_t1k_genotype_confidence(path, sample_override, defaults):
    sample = sample_override or clean_token(path.stem.replace('_genotype', ''))
    out = []
    for line in read_raw_lines(path):
        if not line.strip() or line.startswith('#'):
            continue
        fields = [clean_token(part) for part in line.split('	')]
        if len(fields) < 5:
            continue
        gene = normalize_gene(fields[0])
        if not gene:
            continue
        read_support = 0.0
        for value in (fields[4] if len(fields) > 4 else None, fields[7] if len(fields) > 7 else None):
            numeric = coerce_float(value)
            if numeric is not None and numeric >= 0:
                read_support += numeric
        out.append(build_confidence_entry(sample, gene, None, read_support if read_support > 0 else None, defaults, 't1k_genotype_confidence', path))
    return out


def parse_arcashla_genes_confidence(path, sample_override, defaults):
    with path.open('r', encoding='utf-8') as handle:
        payload = json.load(handle)
    sample = sample_override or clean_token(path.stem.replace('.genes', '')) or path.stem
    out = []
    if not isinstance(payload, dict):
        return out
    for raw_gene, values in payload.items():
        gene = normalize_gene(raw_gene)
        if not gene or not isinstance(values, list):
            continue
        read_support = values[1] if len(values) > 1 else None
        raw_score = values[2] if len(values) > 2 else None
        out.append(build_confidence_entry(sample, gene, raw_score, read_support, defaults, 'arcashla_genes_confidence', path))
    return out


def parse_hlahd_read_confidence(path, sample_override, defaults):
    sample = sample_override or clean_token(path.name.split('_', 1)[0]) or path.stem
    gene_token = clean_token(path.stem.split('_')[-1].replace('.read', ''))
    gene = normalize_gene(gene_token)
    if not gene:
        return []
    for line in read_raw_lines(path):
        if not line.strip() or line.startswith('#'):
            continue
        fields = [clean_token(part) for part in line.split('	')]
        if len(fields) < 2:
            continue
        if fields[0].lower().startswith('r1 only'):
            continue
        read_support = coerce_float(fields[1])
        return [build_confidence_entry(sample, gene, None, read_support, defaults, 'hlahd_read_confidence', path)]
    return [build_confidence_entry(sample, gene, None, None, defaults, 'hlahd_read_confidence', path)]


def parse_kourami_result_confidence(path, sample_override, defaults):
    sample = sample_override or clean_token(path.stem.replace('.kourami', '')) or path.stem
    per_gene = {}
    for line in read_raw_lines(path):
        if not line.strip() or line.startswith('#'):
            continue
        fields = [clean_token(part) for part in line.split('	')]
        if len(fields) < 3:
            continue
        allele = fields[0].split(';')[0]
        gene_token = allele.split('*', 1)[0]
        gene = normalize_gene(gene_token)
        if not gene:
            continue
        raw_score = coerce_float(fields[2])
        read_support = coerce_float(fields[1])
        bucket = per_gene.setdefault(gene, {'scores': [], 'supports': []})
        if raw_score is not None:
            bucket['scores'].append(raw_score)
        if read_support is not None:
            bucket['supports'].append(read_support)
    out = []
    for gene, bucket in per_gene.items():
        raw_score = sum(bucket['scores']) / len(bucket['scores']) if bucket['scores'] else None
        read_support = max(bucket['supports']) if bucket['supports'] else None
        out.append(build_confidence_entry(sample, gene, raw_score, read_support, defaults, 'kourami_result_confidence', path))
    return out


def parse_confidence_file(path, parser_name, sample_override, defaults):
    parser_name = clean_token(parser_name).lower()
    if not parser_name:
        return []
    if parser_name == "optitype_result_objective":
        return parse_optitype_objective_confidence(path, sample_override, defaults)
    if parser_name == "json_confidence_table":
        return parse_json_confidence_table(path, sample_override, defaults)
    if parser_name == "long_confidence_table":
        return parse_long_confidence_table(path, sample_override, defaults)
    if parser_name == "wide_confidence_table":
        return parse_wide_confidence_table(path, sample_override, defaults)
    if parser_name == "t1k_genotype_confidence":
        return parse_t1k_genotype_confidence(path, sample_override, defaults)
    if parser_name == "arcashla_genes_confidence":
        return parse_arcashla_genes_confidence(path, sample_override, defaults)
    if parser_name == "hlahd_read_confidence":
        return parse_hlahd_read_confidence(path, sample_override, defaults)
    if parser_name == "kourami_result_confidence":
        return parse_kourami_result_confidence(path, sample_override, defaults)
    raise ValueError("Unsupported confidence parser: %s" % parser_name)


def resolve_config_paths(config, config_dir):
    """Expand relative paths in config against config_dir (the directory of the config file)."""
    truth = config.get("truth", {})
    if truth.get("path") and not Path(truth["path"]).is_absolute():
        truth["path"] = str((config_dir / truth["path"]).resolve())
    for run in config.get("runs", []):
        for key in ("result_glob", "runtime_glob", "confidence_glob"):
            val = run.get(key)
            if val and not Path(val).is_absolute():
                run[key] = str((config_dir / val).resolve())
    return config


def load_run_records(config):
    records = []
    for raw in config.get("runs", []):
        records.append(
            RunRecord(
                tool=raw["tool"],
                modality=clean_token(raw["modality"]).lower(),
                result_glob=raw["result_glob"],
                parser=raw.get("parser", "auto"),
                runtime_glob=raw.get("runtime_glob"),
                sample_from_path_regex=raw.get("sample_from_path_regex"),
                sample_from_path_group=int(raw.get("sample_from_path_group", 1)),
                gene_filter=[normalize_gene(g) for g in raw.get("gene_filter", [])] or None,
                confidence_glob=raw.get("confidence_glob"),
                confidence_parser=raw.get("confidence_parser", ""),
                confidence_sample_from_path_regex=raw.get("confidence_sample_from_path_regex"),
                confidence_sample_from_path_group=int(raw.get("confidence_sample_from_path_group", 1)),
                confidence_defaults=raw.get("confidence_defaults") or {},
            )
        )
    return records


def truth_exclude_tools(config):
    """Tools that generated the truth table and must NOT be evaluated against it.

    Read from `truth.generated_from` (or `truth.exclude_tools`). Enforcing this
    prevents the circularity of scoring a silver-truth source against its own
    truth. Matching is case-insensitive on the tool name.
    """
    truth_cfg = (config or {}).get("truth", {}) if isinstance(config, dict) else {}
    raw = truth_cfg.get("generated_from") or truth_cfg.get("exclude_tools") or []
    if isinstance(raw, str):
        raw = [raw]
    return {clean_token(tool).lower() for tool in raw if clean_token(tool)}


def apply_truth_source_guardrail(runs, config):
    """Drop run records whose tool generated the truth (anti-circularity)."""
    exclude = truth_exclude_tools(config)
    if not exclude:
        return runs
    kept = []
    dropped = []
    for run in runs:
        if clean_token(run.tool).lower() in exclude:
            dropped.append(run.tool)
        else:
            kept.append(run)
    if dropped:
        sys.stderr.write(
            "[guardrail] Excluding truth-source tool(s) from evaluation to avoid "
            "circularity: %s\n" % ", ".join(sorted(set(dropped)))
        )
    return kept


def build_resource_map(run):
    resource_map = defaultdict(list)
    if not run.runtime_glob:
        return {}
    for file_path in glob.glob(run.runtime_glob, recursive=True):
        path = Path(file_path)
        resource_map[derive_sample_from_path(path, run.sample_from_path_regex, run.sample_from_path_group)].append(path)
    return dict((sample, summarize_resource_logs(paths)) for sample, paths in resource_map.items())


def build_confidence_map(run):
    parser_name = clean_token(run.confidence_parser)
    if not parser_name:
        return {}
    search_glob = run.confidence_glob or run.result_glob
    sample_regex = run.confidence_sample_from_path_regex or run.sample_from_path_regex
    sample_group = run.confidence_sample_from_path_group or run.sample_from_path_group
    confidence_map = {}
    for file_path in sorted(glob.glob(search_glob, recursive=True)):
        path = Path(file_path)
        sample_override = derive_sample_from_path(path, sample_regex, sample_group)
        for entry in parse_confidence_file(path, parser_name, sample_override, run.confidence_defaults):
            key = (entry["sample"], entry["gene"])
            current = confidence_map.get(key)
            if current is None:
                confidence_map[key] = entry
                continue
            current_score = current.get("confidence_score")
            next_score = entry.get("confidence_score")
            if next_score is not None and (current_score is None or next_score > current_score):
                confidence_map[key] = entry
    return confidence_map


def resolve_imgt_hla_version(config):
    truth_cfg = config.get("truth", {})
    benchmark_cfg = config.get("benchmark", {})
    return clean_token(truth_cfg.get("imgt_hla_version") or benchmark_cfg.get("imgt_hla_version") or config.get("imgt_hla_version") or "")


def qualify_allele_for_gene(raw, gene):
    value = clean_token(raw)
    if not value or is_missing(value) or '*' in value:
        return value
    if re.match(r'^[0-9]+(?::[0-9A-Za-z]+)+$', value):
        return f"{gene}*{value}"
    return value


_CIWD_CATALOGUE = None


def get_ciwd_catalogue():
    """Lazily load the vendored CIWD 3.0.0 catalogue (assets/ciwd_3.0.0.tsv).

    Standalone module bin/ciwd.py; loaded robustly whether hla_benchmark is executed as a script
    (its dir is on sys.path) or imported via importlib in tests (dir is not). A missing table yields
    an empty catalogue that classifies every allele as `unknown` — the stratification degrades
    gracefully rather than failing the benchmark. See assets/README_CIWD.md."""
    global _CIWD_CATALOGUE
    if _CIWD_CATALOGUE is None:
        bin_dir = str(Path(__file__).resolve().parent)
        if bin_dir not in sys.path:
            sys.path.insert(0, bin_dir)
        import ciwd
        _CIWD_CATALOGUE = ciwd.load_ciwd()
    return _CIWD_CATALOGUE


def annotate_ciwd(rows):
    """Tag each benchmark row with the CIWD commonness of its truth genotype and typed call.

    Additive: `truth_ciwd_stratum` (rarer of the two truth alleles — the genotype's difficulty),
    per-allele truth categories, and `call_implausible` (1 = the call is not-CIWD/absent → candidate
    error or novel allele). Rows are annotated in place and returned."""
    cat = get_ciwd_catalogue()
    for row in rows:
        t1, t2 = row.get("truth_allele1", ""), row.get("truth_allele2", "")
        row["truth_ciwd_1"] = cat.category(t1)
        row["truth_ciwd_2"] = cat.category(t2)
        row["truth_ciwd_stratum"] = cat.genotype_stratum(t1, t2)
        callable_call = row.get("is_callable", "0") == "1"
        implausible = callable_call and (
            cat.is_implausible(row.get("allele1", "")) or cat.is_implausible(row.get("allele2", ""))
        )
        row["call_implausible"] = "1" if implausible else "0"
    return rows


def build_ciwd_stratified_summary(row_sets):
    """Concordance stratified by CIWD commonness of the truth genotype.

    `row_sets` is a list of (method_label, rows); each row needs truth_allele1/2, modality,
    is_callable, is_correct. Aggregates by (method, modality, ciwd_stratum) — additive to the
    existing overall metrics, answering "does accuracy hold on rare/well-documented alleles, not just
    common ones?" """
    cat = get_ciwd_catalogue()
    agg = defaultdict(lambda: {"gene_rows": 0, "callable": 0, "correct": 0})
    for method, rows in row_sets:
        for row in rows:
            stratum = cat.genotype_stratum(row.get("truth_allele1", ""), row.get("truth_allele2", ""))
            bucket = agg[(method, row.get("modality", ""), stratum)]
            bucket["gene_rows"] += 1
            bucket["callable"] += 1 if row.get("is_callable", "0") == "1" else 0
            bucket["correct"] += 1 if row.get("is_correct", "0") == "1" else 0
    order = {name: i for i, name in enumerate(["common", "intermediate", "well_documented", "not_ciwd", "unknown"])}
    out = []
    for (method, modality, stratum), b in agg.items():
        n = b["gene_rows"]
        out.append({
            "method": method,
            "modality": modality,
            "ciwd_stratum": stratum,
            "gene_rows": n,
            "callable_rate": round(b["callable"] / n, 4) if n else "",
            "overall_correct_call_rate": round(b["correct"] / n, 4) if n else "",
        })
    return sorted(out, key=lambda r: (r["method"], modality_sort_key(r["modality"]), order.get(r["ciwd_stratum"], 99)))


def build_ciwd_plausibility_summary(rows):
    """Per-tool QC: how often a callable call is biologically implausible (not-CIWD/novel), and
    whether such calls are enriched for errors — support for using CIWD as a plausibility flag."""
    agg = defaultdict(lambda: {"callable": 0, "implausible": 0, "implausible_incorrect": 0})
    for row in rows:
        if row.get("is_callable", "0") != "1":
            continue
        b = agg[(row.get("tool", ""), row.get("modality", ""))]
        b["callable"] += 1
        if row.get("call_implausible", "0") == "1":
            b["implausible"] += 1
            if row.get("is_correct", "0") != "1":
                b["implausible_incorrect"] += 1
    out = []
    for (tool, modality), b in agg.items():
        n, imp = b["callable"], b["implausible"]
        out.append({
            "tool": tool,
            "modality": modality,
            "callable_calls": n,
            "implausible_calls": imp,
            "implausible_rate": round(imp / n, 4) if n else "",
            "implausible_incorrect": b["implausible_incorrect"],
            "implausible_error_rate": round(b["implausible_incorrect"] / imp, 4) if imp else "",
        })
    return sorted(out, key=lambda r: (modality_sort_key(r["modality"]), r["tool"]))


def build_harmonized_rows(truth, runs, resolution=2, imgt_hla_version=""):
    rows = []
    for run in runs:
        resources = build_resource_map(run)
        confidence_map = build_confidence_map(run)
        for file_path in sorted(glob.glob(run.result_glob, recursive=True)):
            path = Path(file_path)
            sample_override = derive_sample_from_path(path, run.sample_from_path_regex, run.sample_from_path_group)
            for call in parse_result_file(path, run.parser, sample_override):
                sample = call.sample or sample_override
                gene = normalize_gene(call.gene)
                if run.gene_filter and gene not in run.gene_filter:
                    continue
                truth_pair = truth.get(sample, {}).get(gene)
                if not truth_pair:
                    continue
                raw_typed = [call.allele1, call.allele2]
                raw_truth = [qualify_allele_for_gene(truth_pair[0], gene), qualify_allele_for_gene(truth_pair[1], gene)]
                ambiguity = evaluate_ambiguity(raw_typed, raw_truth, resolution)
                typed_pair = ambiguity["typed_configured"]
                truth_norm = ambiguity["truth_configured"]
                callable_status = all(typed_pair)
                correct_status = callable_status and ambiguity["configured_match"]
                runtime_hours, max_ram_gb = resources.get(sample, (None, None))
                confidence = confidence_map.get((sample, gene), {})
                rows.append({
                    "sample": sample,
                    "modality": run.modality,
                    "tool": run.tool,
                    "gene": gene,
                    "truth_allele1_raw": clean_token(truth_pair[0]),
                    "truth_allele2_raw": clean_token(truth_pair[1]),
                    "allele1_raw": clean_token(call.allele1),
                    "allele2_raw": clean_token(call.allele2),
                    "truth_allele1": truth_norm[0],
                    "truth_allele2": truth_norm[1],
                    "allele1": typed_pair[0],
                    "allele2": typed_pair[1],
                    "truth_allele1_3field": ambiguity["truth_3field"][0],
                    "truth_allele2_3field": ambiguity["truth_3field"][1],
                    "allele1_3field": ambiguity["typed_3field"][0],
                    "allele2_3field": ambiguity["typed_3field"][1],
                    "call_status": "callable" if callable_status else "missing",
                    "correct_status": "correct" if correct_status else "incorrect",
                    "is_callable": "1" if callable_status else "0",
                    "is_correct": "1" if correct_status else "0",
                    "is_correct_2field": "1" if ambiguity["match_2field"] else "0",
                    "is_correct_3field": "1" if ambiguity["match_3field"] else "0",
                    "is_correct_g_group": ambiguity_match_value(ambiguity["g_group_comparable"], ambiguity["match_g_group"]),
                    "is_correct_p_group": ambiguity_match_value(ambiguity["p_group_comparable"], ambiguity["match_p_group"]),
                    "is_ambiguity_compatible": "1" if ambiguity["ambiguity_compatible"] else "0",
                    "is_resolution_compatible": "1" if ambiguity["is_resolution_compatible"] else "0",
                    "compatibility_grade": ambiguity["compatibility_grade"],
                    "match_grade": ambiguity["match_grade"],
                    "imgt_hla_version": imgt_hla_version,
                    "runtime_hours": as_string_number(runtime_hours),
                    "max_ram_gb": as_string_number(max_ram_gb),
                    "confidence_score": as_string_number(confidence.get("confidence_score")),
                    "confidence_source": confidence.get("confidence_source", "missing"),
                    "raw_confidence": as_string_number(confidence.get("raw_confidence")),
                    "read_support": as_string_number(confidence.get("read_support")),
                    "raw_score_family": clean_token(confidence.get("raw_score_family", "")),
                    "raw_score_value": as_string_number(confidence.get("raw_score_value")),
                    "calibrated_probability": "",
                    "cv_calibrated_probability": "",
                    "calibration_method": "",
                    "source_file": str(path),
                })
    return rows


def dedupe_rows(rows):
    seen = {}
    for row in rows:
        key = (row["sample"], row["modality"], row["tool"], row["gene"])
        current = seen.get(key)
        if current is None:
            seen[key] = row
            continue
        rank = (int(row["is_correct"]), int(row["is_callable"]), coerce_float(row.get("confidence_score", "")) or -1.0)
        current_rank = (int(current["is_correct"]), int(current["is_callable"]), coerce_float(current.get("confidence_score", "")) or -1.0)
        if rank > current_rank:
            seen[key] = row
    return sorted(seen.values(), key=lambda r: (modality_sort_key(r["modality"]), r["tool"], r["sample"], gene_sort_key(r["gene"])))


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(dict((field, row.get(field, "")) for field in fieldnames))


def calculate_shared_genes(rows):
    tool_gene_sets = defaultdict(set)
    for row in rows:
        tool_gene_sets[(row["tool"], row["modality"])].add(row["gene"])
    if not tool_gene_sets:
        return []
    return sorted(set.intersection(*tool_gene_sets.values()), key=gene_sort_key)


def sample_set(rows, modality, tool):
    return {row["sample"] for row in rows if row["modality"] == modality and row["tool"] == tool}


def filter_rows(rows, genes=None, samples=None, tools=None):
    out = []
    tool_set = set(tools) if tools else None
    for row in rows:
        if genes and row["gene"] not in genes:
            continue
        if samples is not None and row["sample"] not in samples:
            continue
        if tool_set is not None and row["tool"] not in tool_set:
            continue
        out.append(row)
    return out


def ratio(numerator, denominator):
    return round((float(numerator) / denominator), 4) if denominator else 0.0


def wilson_ci(k, n, z=1.96):
    """Wilson score 95% confidence interval for a proportion k/n.

    Returns (ci_lo, ci_hi) rounded to 4 decimal places.
    Returns (0.0, 0.0) when n == 0.
    """
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / denom
    margin = (z * (p * (1.0 - p) / n + z * z / (4.0 * n * n)) ** 0.5) / denom
    return round(max(0.0, centre - margin), 4), round(min(1.0, centre + margin), 4)


def median_or_none(values):
    return round(statistics.median(values), 4) if values else None


def sample_level_metric(rows, field):
    values = defaultdict(list)
    for row in rows:
        metric = coerce_float(row.get(field, ""))
        if metric is not None:
            values[row["sample"]].append(metric)
    return [statistics.median(v) for v in values.values() if v]


def summarize_entries(entries):
    total = len(entries)
    callable_n = sum(int(row["is_callable"]) for row in entries)
    correct_n = sum(int(row["is_correct"]) for row in entries)
    ci_lo, ci_hi = wilson_ci(correct_n, total)
    return {
        "sample_count": len({row["sample"] for row in entries}),
        "gene_rows": total,
        "callable_rate": ratio(callable_n, total),
        "accuracy_among_callable": ratio(correct_n, callable_n),
        "overall_correct_call_rate": ratio(correct_n, total),
        "overall_correct_call_rate_ci_lo": ci_lo,
        "overall_correct_call_rate_ci_hi": ci_hi,
        "median_runtime_hours": median_or_none(sample_level_metric(entries, "runtime_hours")),
        "median_max_ram_gb": median_or_none(sample_level_metric(entries, "max_ram_gb")),
    }


def aggregate_metrics(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["modality"])].append(row)
    summary = {}
    for key, entries in grouped.items():
        record = summarize_entries(entries)
        record.update({"tool": key[0], "modality": key[1]})
        summary[key] = record
    return summary


def build_per_gene_summary(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["modality"], row["gene"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        record = summarize_entries(entries)
        record.update({"tool": key[0], "modality": key[1], "gene": key[2]})
        out.append(record)
    return out


def build_modality_gene_summary(rows):
    grouped = defaultdict(set)
    for row in rows:
        grouped[row["modality"]].add(row["gene"])
    out = []
    for modality, genes in sorted(grouped.items(), key=lambda item: modality_sort_key(item[0])):
        for gene in sorted(genes, key=gene_sort_key):
            out.append({"modality": modality, "gene": gene})
    return out


def summarize_match_field(entries, field_name):
    comparable = [row for row in entries if clean_token(row.get(field_name, "")) != ""]
    if not comparable:
        return 0, ""
    matches = sum(int(row[field_name]) for row in comparable)
    return len(comparable), ratio(matches, len(comparable))


def ambiguity_summary_record(entries):
    summary = summarize_entries(entries)
    g_rows, g_rate = summarize_match_field(entries, "is_correct_g_group")
    p_rows, p_rate = summarize_match_field(entries, "is_correct_p_group")
    record = {
        "sample_count": summary["sample_count"],
        "gene_rows": summary["gene_rows"],
        "exact_2field_rate": ratio(sum(int(row["is_correct_2field"]) for row in entries), len(entries)),
        "exact_3field_rate": ratio(sum(int(row["is_correct_3field"]) for row in entries), len(entries)),
        "g_group_comparable_rows": g_rows,
        "g_group_match_rate": g_rate,
        "p_group_comparable_rows": p_rows,
        "p_group_match_rate": p_rate,
    }
    return record


def build_ambiguity_summary(rows):
    grouped = defaultdict(list)
    gene_grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["modality"])].append(row)
        gene_grouped[(row["tool"], row["modality"], row["gene"])].append(row)
    summary_rows = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0])):
        record = ambiguity_summary_record(entries)
        record.update({"tool": key[0], "modality": key[1], "imgt_hla_version": entries[0].get("imgt_hla_version", "")})
        summary_rows.append(record)
    per_gene_rows = []
    for key, entries in sorted(gene_grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        record = ambiguity_summary_record(entries)
        record.update({"tool": key[0], "modality": key[1], "gene": key[2], "imgt_hla_version": entries[0].get("imgt_hla_version", "")})
        per_gene_rows.append(record)
    return summary_rows, per_gene_rows


def build_method_ambiguity_summary(*method_row_sets):
    grouped = defaultdict(list)
    gene_grouped = defaultdict(list)
    for row_set in method_row_sets:
        for row in row_set:
            grouped[(row["method"], row["modality"])].append(row)
            gene_grouped[(row["method"], row["modality"], row["gene"])].append(row)
    summary_rows = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0])):
        record = ambiguity_summary_record(entries)
        record.update({"method": key[0], "modality": key[1], "imgt_hla_version": entries[0].get("imgt_hla_version", "")})
        summary_rows.append(record)
    per_gene_rows = []
    for key, entries in sorted(gene_grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        record = ambiguity_summary_record(entries)
        record.update({"method": key[0], "modality": key[1], "gene": key[2], "imgt_hla_version": entries[0].get("imgt_hla_version", "")})
        per_gene_rows.append(record)
    return summary_rows, per_gene_rows


def mean_confidence(entries):
    values = [coerce_float(row.get("confidence_score", "")) for row in entries if row.get("is_callable") == "1" and coerce_float(row.get("confidence_score", "")) is not None]
    return round(statistics.mean(values), 4) if values else None


def mean_probability(entries, field_name):
    values = [coerce_float(row.get(field_name, "")) for row in entries if row.get("is_callable") == "1" and coerce_float(row.get(field_name, "")) is not None]
    return round(statistics.mean(values), 4) if values else None


def confidence_coverage_rate(entries):
    callable_entries = [row for row in entries if row.get("is_callable") == "1"]
    if not callable_entries:
        return 0.0
    with_confidence = sum(1 for row in callable_entries if coerce_float(row.get("confidence_score", "")) is not None)
    return ratio(with_confidence, len(callable_entries))


def confidence_guardrail_settings(config):
    benchmark_cfg = (config or {}).get("benchmark", {}) if isinstance(config, dict) else {}
    raw = benchmark_cfg.get("confidence_guardrail", {}) if isinstance(benchmark_cfg, dict) else {}
    settings = dict(DEFAULT_CONFIDENCE_GUARDRAIL)
    if isinstance(raw, dict):
        settings.update(raw)
    settings["enabled"] = bool(settings.get("enabled", True))
    settings["max_expected_calibration_error_for_boost"] = float(settings.get("max_expected_calibration_error_for_boost", DEFAULT_CONFIDENCE_GUARDRAIL["max_expected_calibration_error_for_boost"]))
    settings["max_brier_score_for_boost"] = float(settings.get("max_brier_score_for_boost", DEFAULT_CONFIDENCE_GUARDRAIL["max_brier_score_for_boost"]))
    settings["min_confidence_coverage_for_boost"] = float(settings.get("min_confidence_coverage_for_boost", DEFAULT_CONFIDENCE_GUARDRAIL["min_confidence_coverage_for_boost"]))
    return settings


def stable_sigmoid(value):
    if value >= 0:
        exp_term = math.exp(-value)
        return 1.0 / (1.0 + exp_term)
    exp_term = math.exp(value)
    return exp_term / (1.0 + exp_term)


def score_value_for_calibration(row):
    raw_score = coerce_float(row.get("raw_score_value", ""))
    if raw_score is not None:
        return raw_score
    return coerce_float(row.get("confidence_score", ""))


def fit_platt_calibrator(scores, labels):
    if not scores:
        return None
    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        constant = 1.0 if positives else 0.0
        return {"method": "constant", "constant": constant}
    mean_score = statistics.mean(scores)
    std_score = statistics.pstdev(scores) if len(scores) > 1 else 0.0
    scale = std_score if std_score > 1e-6 else 1.0
    scaled = [(score - mean_score) / scale for score in scores]
    slope = 0.0
    intercept = math.log((positives + 0.5) / (negatives + 0.5))
    learning_rate = 0.1
    for _ in range(250):
        grad_a = 0.0
        grad_b = 0.0
        for score, label in zip(scaled, labels):
            pred = stable_sigmoid(intercept + slope * score)
            error = pred - label
            grad_a += error
            grad_b += error * score
        intercept -= learning_rate * grad_a / len(scores)
        slope -= learning_rate * grad_b / len(scores)
    return {
        "method": "platt",
        "intercept": intercept,
        "slope": slope,
        "mean_score": mean_score,
        "scale": scale,
    }


def fit_probability_calibrator(entries, settings):
    scored = [row for row in entries if score_value_for_calibration(row) is not None]
    if not scored:
        return None
    scores = [score_value_for_calibration(row) for row in scored]
    labels = [int(row.get("is_correct", "0")) for row in scored]
    method = clean_token((settings or {}).get("method", "platt")).lower() or "platt"
    if method in {"platt", "beta", "isotonic"}:
        calibrator = fit_platt_calibrator(scores, labels)
        if calibrator:
            calibrator["requested_method"] = method
            calibrator["sample_size"] = len(scores)
        return calibrator
    calibrator = fit_platt_calibrator(scores, labels)
    if calibrator:
        calibrator["requested_method"] = "platt"
        calibrator["sample_size"] = len(scores)
    return calibrator


def apply_probability_calibrator(calibrator, score):
    if calibrator is None or score is None:
        return None
    if calibrator.get("method") == "constant":
        return round(float(calibrator["constant"]), 4)
    mean_score = float(calibrator.get("mean_score", 0.0))
    scale = float(calibrator.get("scale", 1.0)) or 1.0
    slope = float(calibrator.get("slope", 0.0))
    intercept = float(calibrator.get("intercept", 0.0))
    scaled = (score - mean_score) / scale
    return round(clip_unit(stable_sigmoid(intercept + slope * scaled)), 4)


def attach_probabilistic_calibration(rows, mode, config):
    settings = probabilistic_calibration_settings(config)
    for row in rows:
        row["calibrated_probability"] = ""
        row["cv_calibrated_probability"] = ""
        row["calibration_method"] = ""
    if mode == "legacy_heuristic":
        return rows
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["modality"])].append(row)
    for key, entries in grouped.items():
        calibrator = fit_probability_calibrator(entries, settings)
        method_label = (calibrator or {}).get("requested_method", settings["method"])
        for row in entries:
            row["calibration_method"] = method_label
            row["calibrated_probability"] = as_string_number(apply_probability_calibrator(calibrator, score_value_for_calibration(row)))
        strategy = settings.get("cv_strategy", "loo")
        for idx, row in enumerate(entries):
            cv_probability = None
            if strategy == "loo" and len(entries) >= 3:
                train_entries = entries[:idx] + entries[idx + 1:]
                cv_calibrator = fit_probability_calibrator(train_entries, settings)
                cv_probability = apply_probability_calibrator(cv_calibrator, score_value_for_calibration(row))
            elif len(entries) >= 2:
                cv_probability = apply_probability_calibrator(calibrator, score_value_for_calibration(row))
            row["cv_calibrated_probability"] = as_string_number(cv_probability)
    return rows


def calibration_stats(entries, bin_count=5, score_field="confidence_score"):
    scored = [row for row in entries if coerce_float(row.get(score_field, "")) is not None]
    n = len(scored)
    if not n:
        return {"n_rows": 0, "mean_confidence": None, "observed_accuracy": None, "brier_score": None, "expected_calibration_error": None}
    grouped = defaultdict(list)
    for row in scored:
        grouped[calibration_bin_label(coerce_float(row.get(score_field, "")), bin_count)[0]].append(row)
    brier = sum((coerce_float(row[score_field]) - int(row["is_correct"])) ** 2 for row in scored) / float(n)
    ece = 0.0
    for entries_in_bin in grouped.values():
        mean_conf = statistics.mean(coerce_float(row[score_field]) for row in entries_in_bin)
        observed = statistics.mean(int(row["is_correct"]) for row in entries_in_bin)
        ece += abs(mean_conf - observed) * len(entries_in_bin)
    ece = ece / float(n)
    return {
        "n_rows": n,
        "mean_confidence": round(statistics.mean(coerce_float(row[score_field]) for row in scored), 4),
        "observed_accuracy": round(statistics.mean(int(row["is_correct"]) for row in scored), 4),
        "brier_score": round(brier, 4),
        "expected_calibration_error": round(ece, 4),
    }


def guarded_confidence_value(reliability, calibrated_confidence, coverage_rate, calibration, settings):
    if calibrated_confidence is None:
        return None, None, "no_confidence"
    if not settings.get("enabled", True):
        return calibrated_confidence, 1.0, "disabled"
    if coverage_rate < settings["min_confidence_coverage_for_boost"]:
        return reliability, 0.0, "coverage_too_low"
    brier = calibration.get("brier_score")
    ece = calibration.get("expected_calibration_error")
    if brier is None or ece is None:
        return reliability, 0.0, "missing_calibration"
    max_cal = max(brier, ece)
    steepness = settings.get("guardrail_steepness", 0)
    if steepness > 0:
        # Soft mode: exponential dampening — no binary cutoff; overrides the hard threshold
        shrink_factor = round(max(0.0, 1.0 - max_cal) * math.exp(-steepness * max_cal), 4)
        status = "applied_soft"
    else:
        # Strict mode (default): binary cutoff at configured threshold
        if max_cal > settings["max_brier_score_for_boost"]:
            return reliability, 0.0, "poor_calibration"
        shrink_factor = round(max(0.0, 1.0 - max_cal), 4)
        status = "applied"
    effective_confidence = round(reliability + shrink_factor * (calibrated_confidence - reliability), 4)
    return effective_confidence, shrink_factor, status


def reliability_estimate(entries, mode, config):
    correct_n = sum(int(row["is_correct"]) for row in entries)
    total = len(entries)
    empirical = ratio(correct_n, total)
    if mode != "bayesian_shrinkage":
        return empirical
    prior = bayesian_shrinkage_settings(config)
    posterior = (correct_n + prior["alpha"]) / (total + prior["alpha"] + prior["beta"]) if total else 0.0
    return round(posterior, 4)


def compute_weight_record(entries, tool, modality, gene=None, settings=None, mode=DEFAULT_BENCHMARK_MODE, config=None):
    summary = summarize_entries(entries)
    reliability = reliability_estimate(entries, mode, config)
    coverage_rate = confidence_coverage_rate(entries)
    calibrated_confidence = mean_confidence(entries) if mode == "legacy_heuristic" else mean_probability(entries, "calibrated_probability")
    cv_calibrated_confidence = mean_probability(entries, "cv_calibrated_probability")
    calibration = calibration_stats(entries)
    calibrated_stats = calibration_stats(entries, score_field="calibrated_probability")
    guardrail_calibration = calibration if mode == "legacy_heuristic" else calibrated_stats
    effective_confidence, guardrail_factor, guardrail_status = guarded_confidence_value(reliability, calibrated_confidence, coverage_rate, guardrail_calibration, settings or DEFAULT_CONFIDENCE_GUARDRAIL)
    alpha = (settings or {}).get("weight_alpha", 0.7)
    beta = (settings or {}).get("weight_beta", 0.3)
    final_weight = reliability if effective_confidence is None else round(alpha * reliability + beta * effective_confidence, 4)
    correct_n = sum(int(row["is_correct"]) for row in entries)
    ci_lo, ci_hi = wilson_ci(correct_n, len(entries))
    raw_families = sorted({clean_token(row.get("raw_score_family", "")) for row in entries if clean_token(row.get("raw_score_family", ""))})
    calibration_method = next((clean_token(row.get("calibration_method", "")) for row in entries if clean_token(row.get("calibration_method", ""))), "")
    record = {
        "tool": tool,
        "modality": modality,
        "benchmark_mode": mode,
        "sample_count": summary["sample_count"],
        "gene_rows": summary["gene_rows"],
        "callable_rate": summary["callable_rate"],
        "confidence_coverage_rate": coverage_rate,
        "base_reliability": reliability,
        "reliability_ci_lo": ci_lo,
        "reliability_ci_hi": ci_hi,
        "reliability_ci_width": round(ci_hi - ci_lo, 4),
        "raw_score_family": ",".join(raw_families),
        "calibration_method": calibration_method,
        "calibration_sample_size": calibrated_stats["n_rows"],
        "calibrated_confidence": "" if calibrated_confidence is None else calibrated_confidence,
        "mean_calibrated_probability": "" if calibrated_confidence is None else calibrated_confidence,
        "cv_mean_calibrated_probability": "" if cv_calibrated_confidence is None else cv_calibrated_confidence,
        "effective_confidence": "" if effective_confidence is None else effective_confidence,
        "guardrail_factor": "" if guardrail_factor is None else guardrail_factor,
        "guardrail_status": guardrail_status,
        "brier_score": "" if calibration["brier_score"] is None else calibration["brier_score"],
        "expected_calibration_error": "" if calibration["expected_calibration_error"] is None else calibration["expected_calibration_error"],
        "calibrated_brier_score": "" if calibrated_stats["brier_score"] is None else calibrated_stats["brier_score"],
        "calibrated_expected_calibration_error": "" if calibrated_stats["expected_calibration_error"] is None else calibrated_stats["expected_calibration_error"],
        "final_weight": final_weight,
        "weight_version": WEIGHT_VERSION,
    }
    if gene is not None:
        record["gene"] = gene
    return record


def build_confidence_weights(rows, config=None, weight_alpha=0.7, weight_beta=0.3, mode=DEFAULT_BENCHMARK_MODE):
    settings = confidence_guardrail_settings(config)
    settings["weight_alpha"] = weight_alpha
    settings["weight_beta"] = weight_beta
    grouped = defaultdict(list)
    gene_grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["modality"])].append(row)
        gene_grouped[(row["tool"], row["modality"], row["gene"])].append(row)
    tool_rows = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0])):
        tool_rows.append(compute_weight_record(entries, key[0], key[1], settings=settings, mode=mode, config=config))
    gene_rows = []
    for key, entries in sorted(gene_grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        gene_rows.append(compute_weight_record(entries, key[0], key[1], gene=key[2], settings=settings, mode=mode, config=config))
    return tool_rows, gene_rows


def classify_failure_mode(row):
    """Classify a benchmark row into a failure mode category.

    Categories (mutually exclusive, in priority order):
      abstained   — tool produced no call
      correct     — both alleles correct at 2-field resolution
      partial     — one of the two called alleles matches a truth allele
      wrong_field2 — field 1 of at least one called allele matches truth but field 2 is wrong
      wrong_field1 — allele group (field 1) is entirely wrong
    """
    if not int(row.get("is_callable", "0")):
        return "abstained"
    if int(row.get("is_correct", "0")):
        return "correct"

    truth_a1 = row.get("truth_allele1", "")
    truth_a2 = row.get("truth_allele2", "")
    call_a1 = row.get("allele1", "")
    call_a2 = row.get("allele2", "")

    # Partial: at least one called allele matches a truth allele exactly
    truth_set = {a for a in (truth_a1, truth_a2) if a}
    if (call_a1 and call_a1 in truth_set) or (call_a2 and call_a2 in truth_set):
        return "partial"

    # Field-1 check (everything before the first colon: e.g. "HLA-A*02" from "HLA-A*02:01")
    def field1(allele):
        if not allele or ":" not in allele:
            return allele
        return allele.split(":")[0]

    truth_f1 = {field1(a) for a in (truth_a1, truth_a2) if a}
    if (call_a1 and field1(call_a1) in truth_f1) or (call_a2 and field1(call_a2) in truth_f1):
        return "wrong_field2"
    return "wrong_field1"


def build_failure_mode_summary(rows):
    """Return per-(tool, modality, gene) failure mode counts and fractions."""
    counts = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)
    for row in rows:
        key = (row["tool"], row["modality"], row["gene"])
        counts[key][classify_failure_mode(row)] += 1
        totals[key] += 1
    result = []
    for key in sorted(counts, key=lambda k: (modality_sort_key(k[1]), k[0], gene_sort_key(k[2]))):
        tool, modality, gene = key
        total = totals[key]
        for mode in ("correct", "partial", "wrong_field2", "wrong_field1", "abstained"):
            count = counts[key].get(mode, 0)
            result.append({
                "tool": tool,
                "modality": modality,
                "gene": gene,
                "failure_mode": mode,
                "count": count,
                "fraction": round(count / total, 4) if total > 0 else 0.0,
            })
    return result


def run_cross_validation(rows, n_folds, config, weight_alpha, weight_beta, mode=DEFAULT_BENCHMARK_MODE):
    """Estimate per-(tool, modality) weight uncertainty via k-fold cross-validation.

    Samples are split into n_folds groups; for each fold, weights are computed
    on the remaining folds. Returns mean ± std of final_weight across folds.
    """
    samples = sorted({r["sample"] for r in rows})
    if len(samples) < n_folds:
        n_folds = max(2, len(samples))
    # Round-robin fold assignment (deterministic, no external dependency)
    chunks = [samples[i::n_folds] for i in range(n_folds)]
    per_fold = []
    for held_out in chunks:
        held_out_set = set(held_out)
        train_rows = [dict(r) for r in rows if r["sample"] not in held_out_set]
        if not train_rows:
            continue
        attach_probabilistic_calibration(train_rows, mode, config)
        tool_weights, _ = build_confidence_weights(train_rows, config=config,
                                                    weight_alpha=weight_alpha,
                                                    weight_beta=weight_beta,
                                                    mode=mode)
        per_fold.append({(w["tool"], w["modality"]): w["final_weight"] for w in tool_weights})
    if not per_fold:
        return []
    all_keys = sorted({k for fw in per_fold for k in fw})
    cv_rows = []
    for (tool, modality) in all_keys:
        fold_weights = [fw[(tool, modality)] for fw in per_fold if (tool, modality) in fw]
        n = len(fold_weights)
        mean_w = round(sum(fold_weights) / n, 4)
        std_w = round(statistics.stdev(fold_weights), 4) if n > 1 else 0.0
        cv_rows.append({
            "tool": tool,
            "modality": modality,
            "n_folds": n,
            "weight_mean": mean_w,
            "weight_std": std_w,
            "weight_cv": round(std_w / mean_w, 4) if mean_w > 0 else 0.0,
            "weight_min": round(min(fold_weights), 4),
            "weight_max": round(max(fold_weights), 4),
        })
    return cv_rows


def load_population_manifest(path):
    """Load sample→superpopulation mapping from a TSV with columns sample/superpopulation."""
    pop_map = {}
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = row.get("sample", "").strip()
            superpop = row.get("superpopulation", "").strip() or row.get("population", "").strip()
            if sample and superpop:
                pop_map[sample] = superpop
    return pop_map


def build_population_weights(rows, config=None, weight_alpha=0.7, weight_beta=0.3):
    """Compute per-superpopulation tool weights by grouping rows on the superpopulation field."""
    by_pop = defaultdict(list)
    for row in rows:
        by_pop[row.get("superpopulation", "unknown")].append(row)
    result = []
    for superpop, pop_rows in sorted(by_pop.items()):
        tool_weights, _ = build_confidence_weights(pop_rows, config=config,
                                                    weight_alpha=weight_alpha,
                                                    weight_beta=weight_beta,
                                                    mode=benchmark_mode(config))
        for w in tool_weights:
            w["superpopulation"] = superpop
        result.extend(tool_weights)
    return result


def build_population_diagnostics(rows, mode):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row.get("superpopulation", "unknown"), row["tool"], row["modality"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (item[0][0], modality_sort_key(item[0][2]), item[0][1])):
        callable_n = sum(int(row["is_callable"]) for row in entries)
        correct_n = sum(int(row["is_correct"]) for row in entries)
        coverage = confidence_coverage_rate(entries)
        out.append({
            "superpopulation": key[0],
            "tool": key[1],
            "modality": key[2],
            "benchmark_mode": mode,
            "sample_count": len({row["sample"] for row in entries}),
            "gene_rows": len(entries),
            "callable_rate": ratio(callable_n, len(entries)),
            "overall_correct_call_rate": ratio(correct_n, len(entries)),
            "confidence_coverage_rate": coverage,
        })
    return out


def build_population_calibration(rows, mode):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row.get("superpopulation", "unknown"), row["tool"], row["modality"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (item[0][0], modality_sort_key(item[0][2]), item[0][1])):
        raw_stats = calibration_stats(entries)
        cal_stats = calibration_stats(entries, score_field="calibrated_probability")
        out.append({
            "superpopulation": key[0],
            "tool": key[1],
            "modality": key[2],
            "benchmark_mode": mode,
            "n_rows": raw_stats["n_rows"],
            "mean_confidence": raw_stats["mean_confidence"],
            "brier_score": raw_stats["brier_score"],
            "expected_calibration_error": raw_stats["expected_calibration_error"],
            "mean_calibrated_probability": cal_stats["mean_confidence"],
            "calibrated_brier_score": cal_stats["brier_score"],
            "calibrated_expected_calibration_error": cal_stats["expected_calibration_error"],
        })
    return out


def build_cross_validation_weight_summary(rows, mode):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["modality"])].append(row)
    summary_rows = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0])):
        cv_values = [coerce_float(row.get("cv_calibrated_probability", "")) for row in entries if coerce_float(row.get("cv_calibrated_probability", "")) is not None]
        cal_values = [coerce_float(row.get("calibrated_probability", "")) for row in entries if coerce_float(row.get("calibrated_probability", "")) is not None]
        summary_rows.append({
            "tool": key[0],
            "modality": key[1],
            "benchmark_mode": mode,
            "n_rows": len(entries),
            "cv_n_rows": len(cv_values),
            "calibration_method": next((clean_token(row.get("calibration_method", "")) for row in entries if clean_token(row.get("calibration_method", ""))), ""),
            "mean_calibrated_probability": round(statistics.mean(cal_values), 4) if cal_values else "",
            "cv_mean_calibrated_probability": round(statistics.mean(cv_values), 4) if cv_values else "",
            "calibrated_brier_score": calibration_stats(entries, score_field="calibrated_probability")["brier_score"] if cal_values else "",
            "calibrated_expected_calibration_error": calibration_stats(entries, score_field="calibrated_probability")["expected_calibration_error"] if cal_values else "",
        })
    return summary_rows


def build_runtime_weight_payload(tool_weights, gene_weights, weight_alpha=0.7, weight_beta=0.3, mode=DEFAULT_BENCHMARK_MODE, meta_params=None):
    gene_strength_priors = defaultdict(dict)
    for row in gene_weights:
        gene_strength_priors[row["modality"]][row["gene"]] = row["base_reliability"]
    payload = {
        "weight_version": WEIGHT_VERSION,
        "benchmark_mode": mode,
        "ensemble_method": "WeightedConsensus",
        "agreement_priors_by_gene": {},
        "gene_strength_priors": {modality: dict(values) for modality, values in gene_strength_priors.items()},
        "meta_score_parameters": meta_params or meta_consensus_parameters({"benchmark": {}}),
        "formula": {
            "final_weight": f"{weight_alpha} * overall_correct_call_rate + {weight_beta} * effective_confidence_score",
            "effective_confidence_score": "base_reliability + shrink_factor * (mean_confidence_score - base_reliability)",
            "fallback": "overall_correct_call_rate when confidence is missing or blocked by guardrail",
        },
        "tool_weights": {},
        "gene_weights": {},
    }
    for row in tool_weights:
        payload["tool_weights"].setdefault(row["tool"], {})[row["modality"]] = {
            "final_weight": row["final_weight"],
            "base_reliability": row["base_reliability"],
            "benchmark_mode": row.get("benchmark_mode", mode),
            "raw_score_family": row.get("raw_score_family", ""),
            "calibration_method": row.get("calibration_method", ""),
            "calibration_sample_size": row.get("calibration_sample_size", ""),
            "mean_calibrated_probability": None if row.get("mean_calibrated_probability", "") == "" else row.get("mean_calibrated_probability"),
            "calibrated_confidence": None if row["calibrated_confidence"] == "" else row["calibrated_confidence"],
            "effective_confidence": None if row["effective_confidence"] == "" else row["effective_confidence"],
            "guardrail_factor": None if row["guardrail_factor"] == "" else row["guardrail_factor"],
            "guardrail_status": row["guardrail_status"],
            "brier_score": None if row["brier_score"] == "" else row["brier_score"],
            "expected_calibration_error": None if row["expected_calibration_error"] == "" else row["expected_calibration_error"],
            "calibrated_brier_score": None if row.get("calibrated_brier_score", "") == "" else row.get("calibrated_brier_score"),
            "calibrated_expected_calibration_error": None if row.get("calibrated_expected_calibration_error", "") == "" else row.get("calibrated_expected_calibration_error"),
            "confidence_coverage_rate": row["confidence_coverage_rate"],
            "sample_count": row["sample_count"],
            "gene_rows": row["gene_rows"],
            "callable_rate": row["callable_rate"],
        }
    for row in gene_weights:
        payload["gene_weights"].setdefault(row["tool"], {}).setdefault(row["modality"], {})[row["gene"]] = {
            "final_weight": row["final_weight"],
            "base_reliability": row["base_reliability"],
            "benchmark_mode": row.get("benchmark_mode", mode),
            "raw_score_family": row.get("raw_score_family", ""),
            "calibration_method": row.get("calibration_method", ""),
            "calibration_sample_size": row.get("calibration_sample_size", ""),
            "mean_calibrated_probability": None if row.get("mean_calibrated_probability", "") == "" else row.get("mean_calibrated_probability"),
            "calibrated_confidence": None if row["calibrated_confidence"] == "" else row["calibrated_confidence"],
            "effective_confidence": None if row["effective_confidence"] == "" else row["effective_confidence"],
            "guardrail_factor": None if row["guardrail_factor"] == "" else row["guardrail_factor"],
            "guardrail_status": row["guardrail_status"],
            "brier_score": None if row["brier_score"] == "" else row["brier_score"],
            "expected_calibration_error": None if row["expected_calibration_error"] == "" else row["expected_calibration_error"],
            "calibrated_brier_score": None if row.get("calibrated_brier_score", "") == "" else row.get("calibrated_brier_score"),
            "calibrated_expected_calibration_error": None if row.get("calibrated_expected_calibration_error", "") == "" else row.get("calibrated_expected_calibration_error"),
            "confidence_coverage_rate": row["confidence_coverage_rate"],
            "sample_count": row["sample_count"],
            "gene_rows": row["gene_rows"],
            "callable_rate": row["callable_rate"],
        }
    return payload


def matched_sample_sets(rows):
    tools = sorted({row["tool"] for row in rows})
    modalities = ordered_modalities(rows)
    matched = {}
    for tool in tools:
        for idx, left in enumerate(modalities):
            for right in modalities[idx + 1:]:
                samples = sample_set(rows, left, tool) & sample_set(rows, right, tool)
                if samples:
                    matched[(tool, left, right)] = samples
    return matched


def build_cohort_overview(rows, truth, shared_genes):
    truth_samples = set(truth.keys())
    tools = sorted({row["tool"] for row in rows})
    out = []
    for modality in ordered_modalities(rows):
        modality_rows = [row for row in rows if row["modality"] == modality]
        modality_samples = {row["sample"] for row in modality_rows} & truth_samples
        for tool in tools:
            tool_samples = {row["sample"] for row in modality_rows if row["tool"] == tool}
            out.append({"modality": modality, "tool": tool, "truth_samples": len(truth_samples), "completed_tool_samples": len(tool_samples), "modality_samples": len(modality_samples), "shared_gene_count": len(shared_genes)})
    return out


def build_sample_disagreements(rows):
    return [row for row in rows if row["is_callable"] == "1" and row["is_correct"] == "0"]


def build_missing_patterns(rows):
    grouped = defaultdict(int)
    for row in rows:
        if row["is_callable"] == "0":
            grouped[(row["tool"], row["modality"], row["gene"])] += 1
    out = []
    for key, missing_n in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        total = sum(1 for row in rows if row["tool"] == key[0] and row["modality"] == key[1] and row["gene"] == key[2])
        out.append({"tool": key[0], "modality": key[1], "gene": key[2], "missing_calls": missing_n, "total_rows": total, "missing_rate": ratio(missing_n, total)})
    return out


def build_truth_mismatches(rows):
    out = []
    for row in rows:
        if row["is_correct"] == "0":
            out.append({"sample": row["sample"], "modality": row["modality"], "tool": row["tool"], "gene": row["gene"], "truth_allele1": row["truth_allele1"], "truth_allele2": row["truth_allele2"], "typed_allele1": row["allele1"], "typed_allele2": row["allele2"], "call_status": row["call_status"]})
    return out


def allele_pair(row):
    if row.get("is_callable") != "1":
        return ("", "")
    return (row.get("allele1", ""), row.get("allele2", ""))


def majority_vote_for_group(entries):
    pair_counts = defaultdict(list)
    for row in entries:
        pair = allele_pair(row)
        if all(pair):
            pair_counts[pair].append(row)
    total_callable = sum(len(rows) for rows in pair_counts.values())
    if not pair_counts:
        return None, 0, 0, 0.0, 0.0, False
    ranked = sorted(pair_counts.items(), key=lambda item: (-len(item[1]), "%s|%s" % item[0]))
    top_pair, top_rows = ranked[0]
    second_n = len(ranked[1][1]) if len(ranked) > 1 else 0
    support_fraction = len(top_rows) / float(total_callable) if total_callable else 0.0
    support_margin = (len(top_rows) - second_n) / float(total_callable) if total_callable else 0.0
    is_tie = len(ranked) > 1 and len(ranked[1][1]) == len(top_rows)
    return top_pair, len(top_rows), total_callable, round(support_fraction, 4), round(support_margin, 4), is_tie


def lookup_weight(payload, tool, modality, gene):
    gene_entry = payload.get("gene_weights", {}).get(tool, {}).get(modality, {}).get(gene)
    if gene_entry and gene_entry.get("final_weight") is not None:
        return float(gene_entry["final_weight"])
    tool_entry = payload.get("tool_weights", {}).get(tool, {}).get(modality, {})
    if tool_entry.get("final_weight") is not None:
        return float(tool_entry["final_weight"])
    return 0.0


def lookup_weight_entry(payload, tool, modality, gene):
    gene_entry = payload.get("gene_weights", {}).get(tool, {}).get(modality, {}).get(gene)
    if gene_entry:
        return gene_entry
    return payload.get("tool_weights", {}).get(tool, {}).get(modality, {})


def load_runtime_weight_override(config):
    override_path = clean_token((config or {}).get("benchmark", {}).get("runtime_weight_override", ""))
    if not override_path:
        return None
    path = Path(override_path)
    if not path.exists():
        raise FileNotFoundError(f"runtime_weight_override not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def classify_consensus_status(top_support, margin, thresholds, has_callable, is_tie):
    if not has_callable:
        return "no_call"
    if is_tie:
        return "no_call"
    if top_support < thresholds["min_support"]:
        return "low_confidence"
    if margin < thresholds["min_margin"]:
        return "low_confidence"
    return "called"


def build_majority_vote_rows(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        pair, agreeing_tools, callable_tools, support_fraction, support_margin, is_tie = majority_vote_for_group(entries)
        truth_1 = entries[0]["truth_allele1"]
        truth_2 = entries[0]["truth_allele2"]
        call_status = "called" if pair else "no_call"
        allele1, allele2 = pair if pair else ("", "")
        is_callable = "1" if pair else "0"
        ambiguity = evaluate_ambiguity([allele1, allele2], [truth_1, truth_2], 2)
        is_correct = "1" if pair and ambiguity["match_2field"] else "0"
        discordance_tag = "technical_conflict" if is_tie else ("no_evidence" if callable_tools == 0 else "consensus_call")
        out.append({
            "sample": key[0],
            "superpopulation": entries[0].get("superpopulation", "unknown"),
            "modality": key[1],
            "gene": key[2],
            "method": "MajorityVote",
            "truth_allele1": truth_1,
            "truth_allele2": truth_2,
            "allele1": allele1,
            "allele2": allele2,
            "truth_allele1_3field": ambiguity["truth_3field"][0],
            "truth_allele2_3field": ambiguity["truth_3field"][1],
            "allele1_3field": ambiguity["typed_3field"][0],
            "allele2_3field": ambiguity["typed_3field"][1],
            "call_status": call_status,
            "is_callable": is_callable,
            "is_correct": is_correct,
            "is_correct_2field": "1" if ambiguity["match_2field"] else "0",
            "is_correct_3field": "1" if ambiguity["match_3field"] else "0",
            "is_correct_g_group": ambiguity_match_value(ambiguity["g_group_comparable"], ambiguity["match_g_group"]),
            "is_correct_p_group": ambiguity_match_value(ambiguity["p_group_comparable"], ambiguity["match_p_group"]),
            "is_ambiguity_compatible": "1" if ambiguity["ambiguity_compatible"] else "0",
            "is_resolution_compatible": "1" if ambiguity["is_resolution_compatible"] else "0",
            "compatibility_grade": ambiguity["compatibility_grade"],
            "match_grade": ambiguity["match_grade"],
            "imgt_hla_version": entries[0].get("imgt_hla_version", ""),
            "agreeing_tools": agreeing_tools,
            "contributing_tools": callable_tools,
            "support_fraction": as_string_number(support_fraction),
            "support_margin": as_string_number(support_margin),
            "discordance_tag": discordance_tag,
        })
    return out


def weighted_vote_for_group(entries, weight_payload):
    pair_weights = defaultdict(lambda: {"weight": 0.0, "tools": []})
    total_weight = 0.0
    contributing_tools = 0
    for row in entries:
        pair = allele_pair(row)
        if not all(pair):
            continue
        weight = lookup_weight(weight_payload, row["tool"], row["modality"], row["gene"])
        if weight <= 0:
            continue
        contributing_tools += 1
        total_weight += weight
        pair_weights[pair]["weight"] += weight
        pair_weights[pair]["tools"].append(row["tool"])
    if not pair_weights:
        # No positive-weight tool voted; fall back to equal weighting for any callable tool
        for row in entries:
            pair = allele_pair(row)
            if not all(pair):
                continue
            pair_weights[pair]["weight"] += 1.0
            pair_weights[pair]["tools"].append(row["tool"])
            total_weight += 1.0
            contributing_tools += 1
        if not pair_weights:
            return None, 0.0, 0, 0.0, 0.0, False, 0
    ranked = sorted(pair_weights.items(), key=lambda item: (-item[1]["weight"], "%s|%s" % item[0]))
    top_pair, top_meta = ranked[0]
    second_weight = ranked[1][1]["weight"] if len(ranked) > 1 else 0.0
    top_support = top_meta["weight"] / total_weight if total_weight else 0.0
    margin = (top_meta["weight"] - second_weight) / total_weight if total_weight else 0.0
    is_tie = len(ranked) > 1 and abs(ranked[0][1]["weight"] - ranked[1][1]["weight"]) < 1e-9
    return top_pair, round(top_meta["weight"], 4), contributing_tools, round(top_support, 4), round(margin, 4), is_tie, len(top_meta["tools"])


def ranked_weighted_candidates(entries, weight_payload):
    pair_weights = defaultdict(lambda: {"weight": 0.0, "tools": [], "tool_weights": {}, "probabilities": []})
    total_weight = 0.0
    for row in entries:
        pair = allele_pair(row)
        if not all(pair):
            continue
        weight = lookup_weight(weight_payload, row["tool"], row["modality"], row["gene"])
        if weight <= 0:
            continue
        pair_weights[pair]["weight"] += weight
        pair_weights[pair]["tools"].append(row["tool"])
        pair_weights[pair]["tool_weights"][row["tool"]] = round(weight, 4)
        probability = coerce_float(row.get("calibrated_probability", "")) or coerce_float(row.get("confidence_score", "")) or 0.0
        pair_weights[pair]["probabilities"].append(probability)
        total_weight += weight
    ranked = sorted(pair_weights.items(), key=lambda item: (-item[1]["weight"], "%s|%s" % item[0]))
    return ranked, total_weight


def consensus_thresholds(config):
    benchmark_cfg = config.get("benchmark", {})
    consensus_cfg = benchmark_cfg.get("consensus", {})
    return {
        "min_support": float(consensus_cfg.get("min_support", 0.55)),
        "min_margin": float(consensus_cfg.get("min_margin", 0.15)),
    }


def pair_string_matches_truth(pair_string, truth_1, truth_2):
    if not pair_string:
        return False
    alleles = pair_string.split("+", 1)
    if len(alleles) != 2:
        return False
    ambiguity = evaluate_ambiguity(alleles, [truth_1, truth_2], 2)
    return bool(ambiguity["match_2field"])


def build_weighted_consensus_rows(rows, weight_payload, config):
    thresholds = consensus_thresholds(config)
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        pair, total_weight, contributing_tools, support_fraction, support_margin, is_tie, top_agreeing_tools = weighted_vote_for_group(entries, weight_payload)
        truth_1 = entries[0]["truth_allele1"]
        truth_2 = entries[0]["truth_allele2"]
        call_status = classify_consensus_status(support_fraction, support_margin, thresholds, contributing_tools > 0, is_tie)
        allele1, allele2 = pair if call_status == "called" and pair else ("", "")
        is_callable = "1" if call_status == "called" and pair else "0"
        ambiguity = evaluate_ambiguity([allele1, allele2], [truth_1, truth_2], 2)
        is_correct = "1" if pair and call_status == "called" and ambiguity["match_2field"] else "0"
        if contributing_tools == 0:
            discordance_tag = "no_evidence"
        elif is_tie:
            discordance_tag = "technical_conflict"
        elif call_status == "low_confidence":
            discordance_tag = "low_evidence_conflict"
        else:
            discordance_tag = "consensus_call"
        out.append({
            "sample": key[0],
            "superpopulation": entries[0].get("superpopulation", "unknown"),
            "modality": key[1],
            "gene": key[2],
            "method": "WeightedConsensus",
            "truth_allele1": truth_1,
            "truth_allele2": truth_2,
            "allele1": allele1,
            "allele2": allele2,
            "truth_allele1_3field": ambiguity["truth_3field"][0],
            "truth_allele2_3field": ambiguity["truth_3field"][1],
            "allele1_3field": ambiguity["typed_3field"][0],
            "allele2_3field": ambiguity["typed_3field"][1],
            "call_status": call_status,
            "is_callable": is_callable,
            "is_correct": is_correct,
            "is_correct_2field": "1" if ambiguity["match_2field"] else "0",
            "is_correct_3field": "1" if ambiguity["match_3field"] else "0",
            "is_correct_g_group": ambiguity_match_value(ambiguity["g_group_comparable"], ambiguity["match_g_group"]),
            "is_correct_p_group": ambiguity_match_value(ambiguity["p_group_comparable"], ambiguity["match_p_group"]),
            "is_ambiguity_compatible": "1" if ambiguity["ambiguity_compatible"] else "0",
            "is_resolution_compatible": "1" if ambiguity["is_resolution_compatible"] else "0",
            "compatibility_grade": ambiguity["compatibility_grade"],
            "match_grade": ambiguity["match_grade"],
            "imgt_hla_version": entries[0].get("imgt_hla_version", ""),
            "agreeing_tools": "" if call_status != "called" else top_agreeing_tools,
            "contributing_tools": contributing_tools,
            "support_fraction": as_string_number(support_fraction),
            "support_margin": as_string_number(support_margin),
            "total_weight": as_string_number(total_weight),
            "discordance_tag": discordance_tag,
        })
    return out


def meta_consensus_parameters(config):
    benchmark_cfg = config.get("benchmark", {}) if isinstance(config, dict) else {}
    raw = benchmark_cfg.get("meta_consensus", {}) if isinstance(benchmark_cfg, dict) else {}
    return {
        "a": float(raw.get("a", 0.45)),
        "b": float(raw.get("b", 0.20)),
        "c": float(raw.get("c", 0.15)),
        "d": float(raw.get("d", 0.10)),
        "e": float(raw.get("e", 0.05)),
        "f": float(raw.get("f", 0.05)),
    }


def allele_resolution_score(allele):
    _gene, parts, _suffix = split_allele_components(allele)
    return min(len(parts), 3) / 3.0 if parts else 0.0


def score_meta_candidate(pair, meta, entries, payload, params):
    callable_tools = len({row["tool"] for row in entries if row.get("is_callable") == "1"})
    supporting_tools = len(set(meta["tools"]))
    support_weight_sum = meta["weight"]
    mean_probability = statistics.mean(meta["probabilities"]) if meta["probabilities"] else 0.0
    agreement_bonus = supporting_tools / float(callable_tools) if callable_tools else 0.0
    base_reliabilities = []
    for tool in set(meta["tools"]):
        entry = lookup_weight_entry(payload, tool, entries[0]["modality"], entries[0]["gene"])
        base = coerce_float(entry.get("base_reliability", "")) if isinstance(entry, dict) else None
        if base is not None:
            base_reliabilities.append(base)
    gene_strength_bonus = statistics.mean(base_reliabilities) if base_reliabilities else 0.0
    ambiguity_bonus = 0.0 if any("/" in allele for allele in pair) else statistics.mean([allele_resolution_score(allele) for allele in pair]) * 0.05
    fragmentation_penalty = 1.0 - agreement_bonus
    meta_score = (
        params["a"] * support_weight_sum +
        params["b"] * mean_probability +
        params["c"] * agreement_bonus +
        params["d"] * gene_strength_bonus +
        params["e"] * ambiguity_bonus -
        params["f"] * fragmentation_penalty
    )
    return {
        "pair": pair,
        "support_weight_sum": round(support_weight_sum, 4),
        "mean_calibrated_probability": round(mean_probability, 4),
        "agreement_bonus": round(agreement_bonus, 4),
        "gene_strength_bonus": round(gene_strength_bonus, 4),
        "ambiguity_bonus": round(ambiguity_bonus, 4),
        "fragmentation_penalty": round(fragmentation_penalty, 4),
        "meta_score": round(meta_score, 4),
        "tool_weights": meta["tool_weights"],
        "supporting_tools": supporting_tools,
    }


def build_meta_consensus_rows(rows, weight_payload, config):
    thresholds = consensus_thresholds(config)
    params = meta_consensus_parameters(config)
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    out = []
    traces = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        ranked, total_weight = ranked_weighted_candidates(entries, weight_payload)
        candidate_scores = [score_meta_candidate(pair, meta, entries, weight_payload, params) for pair, meta in ranked]
        candidate_scores.sort(key=lambda item: (-item["meta_score"], pair_to_string(item["pair"])))
        top = candidate_scores[0] if candidate_scores else None
        runner_up = candidate_scores[1] if len(candidate_scores) > 1 else None
        has_callable = bool(candidate_scores)
        top_support_fraction = (top["support_weight_sum"] / total_weight) if (top and total_weight) else 0.0
        runner_up_weight = runner_up["support_weight_sum"] if runner_up else 0.0
        top_margin = ((top["support_weight_sum"] - runner_up_weight) / total_weight) if (top and total_weight) else 0.0
        call_status = classify_consensus_status(round(top_support_fraction, 4), round(top_margin, 4), thresholds, has_callable, False)
        pair = top["pair"] if top and call_status == "called" else ("", "")
        truth_1 = entries[0]["truth_allele1"]
        truth_2 = entries[0]["truth_allele2"]
        ambiguity = evaluate_ambiguity(list(pair), [truth_1, truth_2], 2)
        is_correct = "1" if pair and call_status == "called" and ambiguity["match_2field"] else "0"
        out.append({
            "sample": key[0],
            "superpopulation": entries[0].get("superpopulation", "unknown"),
            "modality": key[1],
            "gene": key[2],
            "method": "MetaConsensus",
            "truth_allele1": truth_1,
            "truth_allele2": truth_2,
            "allele1": pair[0] if pair else "",
            "allele2": pair[1] if pair else "",
            "truth_allele1_3field": ambiguity["truth_3field"][0],
            "truth_allele2_3field": ambiguity["truth_3field"][1],
            "allele1_3field": ambiguity["typed_3field"][0],
            "allele2_3field": ambiguity["typed_3field"][1],
            "call_status": call_status,
            "is_callable": "1" if pair and call_status == "called" else "0",
            "is_correct": is_correct,
            "is_correct_2field": "1" if ambiguity["match_2field"] else "0",
            "is_correct_3field": "1" if ambiguity["match_3field"] else "0",
            "is_correct_g_group": ambiguity_match_value(ambiguity["g_group_comparable"], ambiguity["match_g_group"]),
            "is_correct_p_group": ambiguity_match_value(ambiguity["p_group_comparable"], ambiguity["match_p_group"]),
            "is_ambiguity_compatible": "1" if ambiguity["ambiguity_compatible"] else "0",
            "is_resolution_compatible": "1" if ambiguity["is_resolution_compatible"] else "0",
            "compatibility_grade": ambiguity["compatibility_grade"],
            "match_grade": ambiguity["match_grade"],
            "imgt_hla_version": entries[0].get("imgt_hla_version", ""),
            "agreeing_tools": "" if not top else top["supporting_tools"],
            "contributing_tools": len({row["tool"] for row in entries if row.get("is_callable") == "1"}),
            "support_fraction": as_string_number(top_support_fraction),
            "support_margin": as_string_number(top_margin),
            "total_weight": as_string_number(top["support_weight_sum"] if top else None),
            "discordance_tag": "consensus_call" if pair and call_status == "called" else ("no_evidence" if not has_callable else "low_evidence_conflict"),
        })
        traces.append({
            "sample": key[0],
            "population": entries[0].get("superpopulation", entries[0].get("population", "unknown")),
            "modality": key[1],
            "gene": key[2],
            "benchmark_mode": benchmark_mode(config),
            "method": "MetaConsensus",
            "winning_pair": pair_to_string(top["pair"]) if top else "",
            "runner_up_pair": pair_to_string(runner_up["pair"]) if runner_up else "",
            "winning_support_fraction": as_string_number(top_support_fraction),
            "winning_support_margin": as_string_number(top_margin),
            "contributing_tools": len({row["tool"] for row in entries if row.get("is_callable") == "1"}),
            "top_tool_weights_json": json.dumps(top["tool_weights"] if top else {}, sort_keys=True, separators=(",", ":")),
            "call_status": call_status,
            "is_correct": is_correct,
            "support_weight_sum": as_string_number(top["support_weight_sum"] if top else None),
            "mean_calibrated_probability": as_string_number(top["mean_calibrated_probability"] if top else None),
            "agreement_bonus": as_string_number(top["agreement_bonus"] if top else None),
            "gene_strength_bonus": as_string_number(top["gene_strength_bonus"] if top else None),
            "ambiguity_bonus": as_string_number(top["ambiguity_bonus"] if top else None),
            "fragmentation_penalty": as_string_number(top["fragmentation_penalty"] if top else None),
            "meta_score": as_string_number(top["meta_score"] if top else None),
        })
    return out, traces


def build_bimodal_consensus_rows(rows, weight_payload, config):
    """Compute WES+RNA joint consensus for samples that have both modalities available.

    Groups by (sample, gene) across wes + rnaseq tool calls, producing both
    BimodalMajorityVote and BimodalWeightedConsensus rows. Only samples with at
    least one callable WES call AND one callable RNA call for a given gene are
    included, so the bimodal set is gene-specific.
    """
    thresholds = consensus_thresholds(config)

    # Build per-gene bimodal eligibility: need ≥1 callable call in each modality
    callable_by_sample_gene_mod = defaultdict(set)
    for row in rows:
        if row.get("is_callable") == "1" and row.get("modality") in ("wes", "rnaseq"):
            callable_by_sample_gene_mod[(row["sample"], row["gene"])].add(row["modality"])
    bimodal_keys = {k for k, mods in callable_by_sample_gene_mod.items() if "wes" in mods and "rnaseq" in mods}

    # Group tool calls for bimodal-eligible (sample, gene) keys
    grouped = defaultdict(list)
    for row in rows:
        if row.get("modality") not in ("wes", "rnaseq"):
            continue
        key = (row["sample"], row["gene"])
        if key in bimodal_keys:
            grouped[key].append(row)

    out = []
    for key in sorted(grouped, key=lambda k: (k[0], gene_sort_key(k[1]))):
        entries = grouped[key]
        sample, gene = key
        truth_1 = entries[0]["truth_allele1"]
        truth_2 = entries[0]["truth_allele2"]
        n_wes = sum(1 for r in entries if r.get("modality") == "wes" and r.get("is_callable") == "1")
        n_rna = sum(1 for r in entries if r.get("modality") == "rnaseq" and r.get("is_callable") == "1")

        # ── BimodalMajorityVote ───────────────────────────────────────────────
        pair_mv, agreeing_mv, callable_mv, support_mv, margin_mv, is_tie_mv = majority_vote_for_group(entries)
        call_status_mv = "called" if pair_mv else "no_call"
        allele1_mv, allele2_mv = pair_mv if pair_mv else ("", "")
        ambiguity_mv = evaluate_ambiguity([allele1_mv, allele2_mv], [truth_1, truth_2], 2)
        is_correct_mv = "1" if pair_mv and ambiguity_mv["match_2field"] else "0"
        discord_mv = "technical_conflict" if is_tie_mv else ("no_evidence" if callable_mv == 0 else "consensus_call")
        out.append({
            "sample": sample,
            "modality": "wes+rnaseq",
            "gene": gene,
            "method": "BimodalMajorityVote",
            "truth_allele1": truth_1,
            "truth_allele2": truth_2,
            "allele1": allele1_mv,
            "allele2": allele2_mv,
            "call_status": call_status_mv,
            "is_callable": "1" if pair_mv else "0",
            "is_correct": is_correct_mv,
            "is_correct_2field": "1" if ambiguity_mv["match_2field"] else "0",
            "is_correct_3field": "1" if ambiguity_mv["match_3field"] else "0",
            "is_ambiguity_compatible": "1" if ambiguity_mv["ambiguity_compatible"] else "0",
            "is_resolution_compatible": "1" if ambiguity_mv["is_resolution_compatible"] else "0",
            "compatibility_grade": ambiguity_mv["compatibility_grade"],
            "agreeing_tools": agreeing_mv,
            "contributing_tools": callable_mv,
            "support_fraction": as_string_number(support_mv),
            "support_margin": as_string_number(margin_mv),
            "n_wes_tools": n_wes,
            "n_rna_tools": n_rna,
            "discordance_tag": discord_mv,
        })

        # ── BimodalWeightedConsensus ──────────────────────────────────────────
        pair_wc, total_wt, contrib_wc, support_wc, margin_wc, is_tie_wc, top_agree_wc = weighted_vote_for_group(entries, weight_payload)
        call_status_wc = classify_consensus_status(support_wc, margin_wc, thresholds, contrib_wc > 0, is_tie_wc)
        allele1_wc, allele2_wc = pair_wc if call_status_wc == "called" and pair_wc else ("", "")
        ambiguity_wc = evaluate_ambiguity([allele1_wc, allele2_wc], [truth_1, truth_2], 2)
        is_correct_wc = "1" if pair_wc and call_status_wc == "called" and ambiguity_wc["match_2field"] else "0"
        if contrib_wc == 0:
            discord_wc = "no_evidence"
        elif is_tie_wc:
            discord_wc = "technical_conflict"
        elif call_status_wc == "low_confidence":
            discord_wc = "low_evidence_conflict"
        else:
            discord_wc = "consensus_call"
        out.append({
            "sample": sample,
            "modality": "wes+rnaseq",
            "gene": gene,
            "method": "BimodalWeightedConsensus",
            "truth_allele1": truth_1,
            "truth_allele2": truth_2,
            "allele1": allele1_wc,
            "allele2": allele2_wc,
            "call_status": call_status_wc,
            "is_callable": "1" if call_status_wc == "called" and pair_wc else "0",
            "is_correct": is_correct_wc,
            "is_correct_2field": "1" if ambiguity_wc["match_2field"] else "0",
            "is_correct_3field": "1" if ambiguity_wc["match_3field"] else "0",
            "is_ambiguity_compatible": "1" if ambiguity_wc["ambiguity_compatible"] else "0",
            "is_resolution_compatible": "1" if ambiguity_wc["is_resolution_compatible"] else "0",
            "compatibility_grade": ambiguity_wc["compatibility_grade"],
            "agreeing_tools": "" if call_status_wc != "called" else top_agree_wc,
            "contributing_tools": contrib_wc,
            "support_fraction": as_string_number(support_wc),
            "support_margin": as_string_number(margin_wc),
            "total_weight": as_string_number(total_wt),
            "n_wes_tools": n_wes,
            "n_rna_tools": n_rna,
            "discordance_tag": discord_wc,
        })
    return out


def build_bimodal_accuracy_comparison(main_rows, bimodal_rows, weight_payload, config):
    """Compare WES-only, RNA-only, and bimodal joint consensus on bimodal-eligible samples.

    Restricts per-modality baselines to the same sample set used in bimodal rows so
    the comparison is apples-to-apples.
    """
    bimodal_samples = {row["sample"] for row in bimodal_rows}

    wes_rows = [r for r in main_rows if r["sample"] in bimodal_samples and r.get("modality") == "wes"]
    rna_rows = [r for r in main_rows if r["sample"] in bimodal_samples and r.get("modality") == "rnaseq"]
    wes_mv = build_majority_vote_rows(wes_rows)
    wes_wc = build_weighted_consensus_rows(wes_rows, weight_payload, config)
    rna_mv = build_majority_vote_rows(rna_rows)
    rna_wc = build_weighted_consensus_rows(rna_rows, weight_payload, config)
    bi_mv = [r for r in bimodal_rows if r["method"] == "BimodalMajorityVote"]
    bi_wc = [r for r in bimodal_rows if r["method"] == "BimodalWeightedConsensus"]

    out = []
    for label, subset, modality_label in [
        ("WES_MajorityVote",           wes_mv, "wes"),
        ("WES_WeightedConsensus",      wes_wc, "wes"),
        ("RNA_MajorityVote",           rna_mv, "rnaseq"),
        ("RNA_WeightedConsensus",      rna_wc, "rnaseq"),
        ("Bimodal_MajorityVote",       bi_mv,  "wes+rnaseq"),
        ("Bimodal_WeightedConsensus",  bi_wc,  "wes+rnaseq"),
    ]:
        total = len(subset)
        callable_n = sum(int(r["is_callable"]) for r in subset)
        correct_n = sum(int(r["is_correct"]) for r in subset)
        sample_n = len({r["sample"] for r in subset})
        ci_lo, ci_hi = wilson_ci(correct_n, total)
        out.append({
            "comparison": label,
            "modality": modality_label,
            "sample_count": sample_n,
            "gene_rows": total,
            "callable_rate": ratio(callable_n, total),
            "accuracy_among_callable": ratio(correct_n, callable_n),
            "overall_correct_call_rate": ratio(correct_n, total),
            "overall_correct_call_rate_ci_lo": ci_lo,
            "overall_correct_call_rate_ci_hi": ci_hi,
        })
    return out


def summarize_consensus_method(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["method"], row["modality"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0])):
        total = len(entries)
        callable_n = sum(int(row["is_callable"]) for row in entries)
        correct_n = sum(int(row["is_correct"]) for row in entries)
        ci_lo, ci_hi = wilson_ci(correct_n, total)
        out.append({
            "method": key[0],
            "modality": key[1],
            "sample_count": len({row["sample"] for row in entries}),
            "gene_rows": total,
            "callable_rate": ratio(callable_n, total),
            "accuracy_among_callable": ratio(correct_n, callable_n),
            "overall_correct_call_rate": ratio(correct_n, total),
            "overall_correct_call_rate_ci_lo": ci_lo,
            "overall_correct_call_rate_ci_hi": ci_hi,
        })
    return out


def build_method_comparison(single_tool_summary, *consensus_row_sets):
    out = []
    for record in sorted(single_tool_summary.values(), key=lambda row: (modality_sort_key(row["modality"]), row["tool"])):
        out.append({
            "method": record["tool"],
            "method_type": "single_tool",
            "modality": record["modality"],
            "sample_count": record["sample_count"],
            "gene_rows": record["gene_rows"],
            "callable_rate": record["callable_rate"],
            "accuracy_among_callable": record["accuracy_among_callable"],
            "overall_correct_call_rate": record["overall_correct_call_rate"],
            "overall_correct_call_rate_ci_lo": record.get("overall_correct_call_rate_ci_lo", ""),
            "overall_correct_call_rate_ci_hi": record.get("overall_correct_call_rate_ci_hi", ""),
        })
    for row_set in consensus_row_sets:
        summaries = summarize_consensus_method(row_set)
        for row in summaries:
            method_type = "baseline" if row["method"] == "MajorityVote" else "ensemble"
            out.append(dict(row, method_type=method_type))
    return out


def build_method_per_gene(single_tool_rows, *consensus_row_sets):
    out = []
    out.extend({"method": row["tool"], "method_type": "single_tool", **{k: v for k, v in row.items() if k != "tool"}} for row in single_tool_rows)
    for consensus_rows in consensus_row_sets:
        method_name = consensus_rows[0]["method"] if consensus_rows else "Consensus"
        method_type = "baseline" if method_name == "MajorityVote" else "ensemble"
        grouped = defaultdict(list)
        for row in consensus_rows:
            grouped[(row["method"], row["modality"], row["gene"])].append(row)
        for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
            total = len(entries)
            callable_n = sum(int(row["is_callable"]) for row in entries)
            correct_n = sum(int(row["is_correct"]) for row in entries)
            out.append({
                "method": method_name,
                "method_type": method_type,
                "modality": key[1],
                "gene": key[2],
                "sample_count": len({row["sample"] for row in entries}),
                "gene_rows": total,
                "callable_rate": ratio(callable_n, total),
                "accuracy_among_callable": ratio(correct_n, callable_n),
                "overall_correct_call_rate": ratio(correct_n, total),
            })
    return out


def ensemble_ablation_settings(config):
    benchmark_cfg = (config or {}).get("benchmark", {}) if isinstance(config, dict) else {}
    raw = benchmark_cfg.get("ensemble_ablation", {}) if isinstance(benchmark_cfg, dict) else {}
    if not isinstance(raw, dict):
        return {"enabled": False, "tool_subsets": []}
    enabled = bool(raw.get("enabled"))
    subsets = []
    for subset in raw.get("tool_subsets", []) or []:
        if not isinstance(subset, dict):
            continue
        name = clean_token(subset.get("name", ""))
        tools = [clean_token(tool) for tool in subset.get("tools", []) if clean_token(tool)]
        if not name or not tools:
            continue
        subsets.append({"name": name, "tools": tools})
    return {"enabled": enabled and bool(subsets), "tool_subsets": subsets}


def locus_expert_consensus_settings(config):
    benchmark_cfg = (config or {}).get("benchmark", {}) if isinstance(config, dict) else {}
    raw = benchmark_cfg.get("locus_expert_consensus", {}) if isinstance(benchmark_cfg, dict) else {}
    if not isinstance(raw, dict):
        return {"enabled": False, "parent_method": "weighted_consensus", "panel_sets": []}
    parent_method = clean_token(raw.get("parent_method", "weighted_consensus")).lower() or "weighted_consensus"
    panel_sets = []
    for panel in raw.get("panel_sets", []) or []:
        if not isinstance(panel, dict):
            continue
        name = clean_token(panel.get("name", ""))
        raw_subsets = panel.get("tool_subsets_by_gene", {}) or {}
        tool_subsets_by_gene = {}
        if isinstance(raw_subsets, dict):
            for gene, tools in raw_subsets.items():
                norm_gene = normalize_gene(gene)
                norm_tools = [clean_token(tool) for tool in (tools or []) if clean_token(tool)]
                tool_subsets_by_gene[norm_gene] = norm_tools
        if name and tool_subsets_by_gene:
            panel_sets.append({"name": name, "tool_subsets_by_gene": tool_subsets_by_gene})
    return {"enabled": bool(raw.get("enabled")) and bool(panel_sets), "parent_method": parent_method, "panel_sets": panel_sets}


def champion_challenger_settings(config):
    benchmark_cfg = (config or {}).get("benchmark", {}) if isinstance(config, dict) else {}
    raw = benchmark_cfg.get("champion_challenger", {}) if isinstance(benchmark_cfg, dict) else {}
    if not isinstance(raw, dict):
        return {
            "enabled": False,
            "champion_by_gene": {},
            "fallback_method": "weighted_consensus",
            "override_policy": {
                "min_challenger_support_fraction": 0.65,
                "min_challenger_margin": 0.20,
                "min_supporting_tools": 2,
                "require_non_ambiguity_override": True,
            },
        }
    override_policy = raw.get("override_policy", {}) if isinstance(raw.get("override_policy", {}), dict) else {}
    return {
        "enabled": bool(raw.get("enabled")),
        "champion_by_gene": {normalize_gene(gene): clean_token(tool) for gene, tool in (raw.get("champion_by_gene", {}) or {}).items() if clean_token(tool)},
        "fallback_method": clean_token(raw.get("fallback_method", "weighted_consensus")).lower() or "weighted_consensus",
        "override_policy": {
            "min_challenger_support_fraction": float(override_policy.get("min_challenger_support_fraction", 0.65)),
            "min_challenger_margin": float(override_policy.get("min_challenger_margin", 0.20)),
            "min_supporting_tools": int(override_policy.get("min_supporting_tools", 2)),
            "require_non_ambiguity_override": bool(override_policy.get("require_non_ambiguity_override", True)),
        },
    }


def validate_champion_challenger_settings(settings, rows, benchmark_genes):
    if not settings.get("enabled"):
        return settings
    fallback_method = settings.get("fallback_method", "weighted_consensus")
    if fallback_method != "weighted_consensus":
        raise ValueError("champion_challenger.fallback_method must be weighted_consensus")
    allowed_genes = set(benchmark_genes or [])
    champion_by_gene = settings.get("champion_by_gene", {})
    missing_genes = sorted(allowed_genes.difference(champion_by_gene), key=gene_sort_key)
    if missing_genes:
        raise ValueError("champion_challenger is missing champions for genes: %s" % ",".join(missing_genes))
    available_tools = {clean_token(row.get("tool")) for row in rows if clean_token(row.get("tool"))}
    for gene, champion in champion_by_gene.items():
        if gene not in allowed_genes:
            raise ValueError("champion_challenger contains unsupported gene: %s" % gene)
        if champion not in available_tools:
            raise ValueError("champion_challenger references unknown champion tool %s for gene %s" % (champion, gene))
    return settings


def consensus_gating_settings(config):
    benchmark_cfg = (config or {}).get("benchmark", {}) if isinstance(config, dict) else {}
    raw = benchmark_cfg.get("consensus_gating", {}) if isinstance(benchmark_cfg, dict) else {}
    if not isinstance(raw, dict):
        return {
            "enabled": False,
            "hard_locus_rule": {
                "min_distinct_pairs_for_hard": 4,
                "min_support_margin_for_easy": 0.35,
                "min_support_fraction_for_easy": 0.50,
                "trigger_on_majority_weighted_disagreement": True,
                "trigger_on_duplicated_top_pair_with_alternative": True,
            },
            "easy_policy": "majority_vote",
            "hard_policy": "champion_challenger",
        }
    rule = raw.get("hard_locus_rule", {}) if isinstance(raw.get("hard_locus_rule", {}), dict) else {}
    return {
        "enabled": bool(raw.get("enabled")),
        "hard_locus_rule": {
            "min_distinct_pairs_for_hard": int(rule.get("min_distinct_pairs_for_hard", 4)),
            "min_support_margin_for_easy": float(rule.get("min_support_margin_for_easy", 0.35)),
            "min_support_fraction_for_easy": float(rule.get("min_support_fraction_for_easy", 0.50)),
            "trigger_on_majority_weighted_disagreement": bool(rule.get("trigger_on_majority_weighted_disagreement", True)),
            "trigger_on_duplicated_top_pair_with_alternative": bool(rule.get("trigger_on_duplicated_top_pair_with_alternative", True)),
        },
        "easy_policy": clean_token(raw.get("easy_policy", "majority_vote")).lower() or "majority_vote",
        "hard_policy": clean_token(raw.get("hard_policy", "champion_challenger")).lower() or "champion_challenger",
    }


def validate_consensus_gating_settings(settings, config):
    if not settings.get("enabled"):
        return settings
    if settings.get("easy_policy") != "majority_vote":
        raise ValueError("consensus_gating.easy_policy must be majority_vote")
    if settings.get("hard_policy") != "champion_challenger":
        raise ValueError("consensus_gating.hard_policy must be champion_challenger")
    rule = settings.get("hard_locus_rule", {})
    if rule.get("min_distinct_pairs_for_hard", 0) < 2:
        raise ValueError("consensus_gating.hard_locus_rule.min_distinct_pairs_for_hard must be >= 2")
    for key in ("min_support_margin_for_easy", "min_support_fraction_for_easy"):
        value = float(rule.get(key, 0.0))
        if value < 0.0 or value > 1.0:
            raise ValueError("consensus_gating.hard_locus_rule.%s must be between 0 and 1" % key)
    champion_settings = champion_challenger_settings(config)
    if not champion_settings.get("enabled"):
        raise ValueError("consensus_gating.enabled requires champion_challenger.enabled = true")
    return settings


def is_duplicated_pair(pair):
    return bool(pair and len(pair) == 2 and pair[0] and pair[0] == pair[1])


def compute_gating_features(entries, majority_row, weighted_row):
    callable_rows = [row for row in entries if row.get("is_callable") == "1" and pair_to_string((row.get("allele1", ""), row.get("allele2", "")))]
    distinct_pairs = sorted({pair_to_string((row["allele1"], row["allele2"])) for row in callable_rows})
    weighted_pair = pair_to_string((weighted_row.get("allele1", ""), weighted_row.get("allele2", "")))
    majority_pair = pair_to_string((majority_row.get("allele1", ""), majority_row.get("allele2", "")))
    alternative_nonduplicated_exists = any(
        pair and pair != weighted_pair and not is_duplicated_pair(pair.split("+"))
        for pair in distinct_pairs
    )
    return {
        "distinct_pair_count": len(distinct_pairs),
        "callable_tool_count": len(callable_rows),
        "weighted_winning_pair": weighted_pair,
        "weighted_support_fraction": coerce_float(weighted_row.get("support_fraction", "")) or 0.0,
        "weighted_support_margin": coerce_float(weighted_row.get("support_margin", "")) or 0.0,
        "majority_pair": majority_pair,
        "majority_weighted_disagree": bool(weighted_pair and majority_pair and weighted_pair != majority_pair),
        "weighted_top_pair_is_duplicated": is_duplicated_pair((weighted_row.get("allele1", ""), weighted_row.get("allele2", ""))),
        "alternative_nonduplicated_exists": alternative_nonduplicated_exists,
    }


def classify_hard_locus(features, settings):
    rule = settings["hard_locus_rule"]
    reasons = []
    if features["distinct_pair_count"] >= rule["min_distinct_pairs_for_hard"]:
        reasons.append("distinct_pairs")
    if features["weighted_support_margin"] < rule["min_support_margin_for_easy"]:
        reasons.append("low_margin")
    if features["weighted_support_fraction"] < rule["min_support_fraction_for_easy"]:
        reasons.append("low_support")
    if rule["trigger_on_majority_weighted_disagreement"] and features["majority_weighted_disagree"]:
        reasons.append("majority_weighted_disagree")
    if (
        rule["trigger_on_duplicated_top_pair_with_alternative"]
        and features["weighted_top_pair_is_duplicated"]
        and features["alternative_nonduplicated_exists"]
    ):
        reasons.append("duplicated_top_pair_with_alternative")
    return ("hard" if reasons else "easy"), reasons


def ambiguity_override_allowed(pair, require_non_ambiguity_override):
    if not require_non_ambiguity_override:
        return True
    return not any("/" in allele for allele in pair if allele)


def build_champion_challenger_outputs(rows, runtime_weights, config, mode, benchmark_genes):
    settings = validate_champion_challenger_settings(champion_challenger_settings(config), rows, benchmark_genes)
    if not settings.get("enabled"):
        return [], [], [], [], []
    weighted_baseline = build_weighted_consensus_rows(rows, runtime_weights, config)
    weighted_index = {(row["sample"], row["modality"], row["gene"]): row for row in weighted_baseline}
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    call_rows = []
    trace_rows = []
    policy = settings["override_policy"]
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        sample, modality, gene = key
        champion_tool = settings["champion_by_gene"].get(gene, "")
        champion_entry = next((row for row in entries if row["tool"] == champion_tool), None)
        fallback_row = dict(weighted_index.get(key, empty_consensus_call_row(entries, "ChampionChallenger", [])))
        fallback_row["method"] = "ChampionChallenger"
        fallback_row["champion_tool"] = champion_tool
        fallback_row["override_triggered"] = "0"
        fallback_row["override_reason"] = "fallback_weighted_consensus"
        fallback_row["population"] = entries[0].get("superpopulation", entries[0].get("population", "unknown"))
        ranked, total_weight = ranked_weighted_candidates(entries, runtime_weights)
        challenger = None
        for pair, meta in ranked:
            champion_pair = allele_pair(champion_entry) if champion_entry and champion_entry.get("is_callable") == "1" else ("", "")
            if pair != champion_pair:
                challenger = (pair, meta)
                break
        if not champion_entry or champion_entry.get("is_callable") != "1":
            call_rows.append(fallback_row)
            trace_rows.append({
                "sample": sample,
                "population": entries[0].get("superpopulation", entries[0].get("population", "unknown")),
                "modality": modality,
                "gene": gene,
                "method": "ChampionChallenger",
                "champion_tool": champion_tool,
                "champion_pair": "" if not champion_entry else pair_to_string(allele_pair(champion_entry)),
                "challenger_pair": "" if not challenger else pair_to_string(challenger[0]),
                "override_triggered": "0",
                "override_reason": "champion_missing",
                "call_status": fallback_row["call_status"],
                "is_correct": fallback_row["is_correct"],
                "challenger_support_fraction": "",
                "challenger_support_margin": "",
                "supporting_tools": "",
            })
            continue
        champion_pair = allele_pair(champion_entry)
        champion_ambiguity = evaluate_ambiguity(list(champion_pair), [entries[0]["truth_allele1"], entries[0]["truth_allele2"]], 2)
        chosen_pair = champion_pair
        chosen_status = "called"
        override_triggered = False
        override_reason = "champion_retained"
        challenger_support_fraction = None
        challenger_support_margin = None
        supporting_tools = 0
        top_weights = {}
        if challenger and total_weight:
            pair, meta = challenger
            challenger_support_fraction = round(meta["weight"] / total_weight, 4)
            runner_up_weight = ranked[1][1]["weight"] if len(ranked) > 1 and ranked[0][0] == champion_pair else ranked[0][1]["weight"]
            challenger_support_margin = round((meta["weight"] - runner_up_weight) / total_weight, 4)
            supporting_tools = len(set(meta["tools"]))
            top_weights = meta["tool_weights"]
            if (
                challenger_support_fraction >= policy["min_challenger_support_fraction"]
                and challenger_support_margin >= policy["min_challenger_margin"]
                and supporting_tools >= policy["min_supporting_tools"]
                and ambiguity_override_allowed(pair, policy["require_non_ambiguity_override"])
            ):
                chosen_pair = pair
                override_triggered = True
                override_reason = "challenger_override"
        ambiguity = evaluate_ambiguity(list(chosen_pair), [entries[0]["truth_allele1"], entries[0]["truth_allele2"]], 2)
        call_rows.append({
            "sample": sample,
            "population": entries[0].get("superpopulation", entries[0].get("population", "unknown")),
            "superpopulation": entries[0].get("superpopulation", "unknown"),
            "modality": modality,
            "gene": gene,
            "method": "ChampionChallenger",
            "truth_allele1": entries[0]["truth_allele1"],
            "truth_allele2": entries[0]["truth_allele2"],
            "allele1": chosen_pair[0] if chosen_status == "called" else "",
            "allele2": chosen_pair[1] if chosen_status == "called" else "",
            "truth_allele1_3field": ambiguity["truth_3field"][0],
            "truth_allele2_3field": ambiguity["truth_3field"][1],
            "allele1_3field": ambiguity["typed_3field"][0],
            "allele2_3field": ambiguity["typed_3field"][1],
            "call_status": chosen_status,
            "is_callable": "1" if chosen_status == "called" else "0",
            "is_correct": "1" if ambiguity["match_2field"] and chosen_status == "called" else "0",
            "is_correct_2field": "1" if ambiguity["match_2field"] else "0",
            "is_correct_3field": "1" if ambiguity["match_3field"] else "0",
            "is_correct_g_group": ambiguity_match_value(ambiguity["g_group_comparable"], ambiguity["match_g_group"]),
            "is_correct_p_group": ambiguity_match_value(ambiguity["p_group_comparable"], ambiguity["match_p_group"]),
            "is_ambiguity_compatible": "1" if ambiguity["ambiguity_compatible"] else "0",
            "is_resolution_compatible": "1" if ambiguity["is_resolution_compatible"] else "0",
            "compatibility_grade": ambiguity["compatibility_grade"],
            "match_grade": ambiguity["match_grade"],
            "imgt_hla_version": entries[0].get("imgt_hla_version", ""),
            "agreeing_tools": supporting_tools if override_triggered else 1,
            "contributing_tools": len({row["tool"] for row in entries if row.get("is_callable") == "1"}),
            "support_fraction": as_string_number(challenger_support_fraction if override_triggered else 1.0),
            "support_margin": as_string_number(challenger_support_margin if override_triggered else 1.0),
            "total_weight": as_string_number((challenger[1]["weight"] if override_triggered and challenger else lookup_weight(runtime_weights, champion_tool, modality, gene))),
            "discordance_tag": "challenger_override" if override_triggered else "champion_call",
            "champion_tool": champion_tool,
            "override_triggered": "1" if override_triggered else "0",
            "override_reason": override_reason,
        })
        trace_rows.append({
            "sample": sample,
            "population": entries[0].get("superpopulation", entries[0].get("population", "unknown")),
            "modality": modality,
            "gene": gene,
            "method": "ChampionChallenger",
            "champion_tool": champion_tool,
            "champion_pair": pair_to_string(champion_pair),
            "challenger_pair": "" if not challenger else pair_to_string(challenger[0]),
            "override_triggered": "1" if override_triggered else "0",
            "override_reason": override_reason,
            "call_status": chosen_status,
            "is_correct": "1" if ambiguity["match_2field"] and chosen_status == "called" else "0",
            "challenger_support_fraction": as_string_number(challenger_support_fraction),
            "challenger_support_margin": as_string_number(challenger_support_margin),
            "supporting_tools": supporting_tools,
            "top_tool_weights_json": json.dumps(top_weights, sort_keys=True, separators=(",", ":")),
        })
    comparison_rows = [dict(row, method_type="ensemble") for row in summarize_consensus_method(call_rows)]
    per_gene_rows = build_method_per_gene([], call_rows)
    summary_rows = [{
        "champion_by_gene": settings["champion_by_gene"],
        "fallback_method": settings["fallback_method"],
        "override_policy": settings["override_policy"],
        "overall_metrics": comparison_rows,
        "per_gene_metrics": per_gene_rows,
    }]
    return call_rows, trace_rows, comparison_rows, per_gene_rows, summary_rows


def summarize_override_effects(trace_rows, truth_index):
    override_count = 0
    corrective = 0
    harmful = 0
    neutral = 0
    for row in trace_rows:
        if row.get("override_triggered") != "1":
            continue
        override_count += 1
        truth_1, truth_2 = truth_index.get((row["sample"], row["modality"], row["gene"]), ("", ""))
        champion_correct = pair_string_matches_truth(row.get("champion_pair", ""), truth_1, truth_2)
        challenger_correct = pair_string_matches_truth(row.get("challenger_pair", ""), truth_1, truth_2)
        if not champion_correct and challenger_correct:
            corrective += 1
        elif champion_correct and not challenger_correct:
            harmful += 1
        else:
            neutral += 1
    return {
        "override_count": override_count,
        "corrective_override_count": corrective,
        "harmful_override_count": harmful,
        "neutral_override_count": neutral,
    }


def run_weighted_threshold_sweep(rows, runtime_weights, config):
    settings = weighted_threshold_sweep_settings(config)
    modality = sweep_target_modality(rows)
    target_rows = [row for row in rows if row["modality"] == modality]
    if not settings.get("enabled") or not target_rows or not modality:
        return []
    out = []
    for min_support in settings["min_support_values"]:
        for min_margin in settings["min_margin_values"]:
            temp_config = copy.deepcopy(config)
            temp_config.setdefault("benchmark", {})
            temp_config["benchmark"].setdefault("consensus", {})
            temp_config["benchmark"]["consensus"]["min_support"] = min_support
            temp_config["benchmark"]["consensus"]["min_margin"] = min_margin
            weighted_rows = build_weighted_consensus_rows(target_rows, runtime_weights, temp_config)
            summary_rows = summarize_consensus_method(weighted_rows)
            summary = next((row for row in summary_rows if row["modality"] == modality), None)
            if summary is None:
                summary = {
                    "method": "WeightedConsensus",
                    "modality": modality,
                    "sample_count": len({row["sample"] for row in target_rows}),
                    "gene_rows": len({(row["sample"], row["gene"]) for row in target_rows}),
                    "callable_rate": 0.0,
                    "accuracy_among_callable": 0.0,
                    "overall_correct_call_rate": 0.0,
                    "overall_correct_call_rate_ci_lo": 0.0,
                    "overall_correct_call_rate_ci_hi": 0.0,
                }
            out.append({
                "method": "WeightedConsensus",
                "modality": modality,
                "min_support": min_support,
                "min_margin": min_margin,
                "sample_count": summary["sample_count"],
                "gene_rows": summary["gene_rows"],
                "callable_rate": summary["callable_rate"],
                "accuracy_among_callable": summary["accuracy_among_callable"],
                "overall_correct_call_rate": summary["overall_correct_call_rate"],
                "overall_correct_call_rate_ci_lo": summary["overall_correct_call_rate_ci_lo"],
                "overall_correct_call_rate_ci_hi": summary["overall_correct_call_rate_ci_hi"],
            })
    return out


def run_champion_override_sweep(rows, runtime_weights, config, benchmark_genes, mode):
    settings = champion_override_sweep_settings(config)
    modality = sweep_target_modality(rows)
    target_rows = [row for row in rows if row["modality"] == modality]
    if (
        not settings.get("enabled")
        or not target_rows
        or not modality
        or not settings["champion_by_gene"]
        or not settings["min_challenger_support_fraction_values"]
        or not settings["min_challenger_margin_values"]
        or not settings["min_supporting_tools_values"]
    ):
        return []
    truth_index = {}
    for row in target_rows:
        truth_index[(row["sample"], row["modality"], row["gene"])] = (row["truth_allele1"], row["truth_allele2"])
    out = []
    for min_support in settings["min_challenger_support_fraction_values"]:
        for min_margin in settings["min_challenger_margin_values"]:
            for min_tools in settings["min_supporting_tools_values"]:
                temp_config = copy.deepcopy(config)
                temp_config.setdefault("benchmark", {})
                temp_config["benchmark"]["champion_challenger"] = {
                    "enabled": True,
                    "champion_by_gene": settings["champion_by_gene"],
                    "fallback_method": "weighted_consensus",
                    "override_policy": {
                        "min_challenger_support_fraction": min_support,
                        "min_challenger_margin": min_margin,
                        "min_supporting_tools": min_tools,
                        "require_non_ambiguity_override": settings["require_non_ambiguity_override"],
                    },
                }
                call_rows, trace_rows, comparison_rows, _per_gene_rows, _summary_rows = build_champion_challenger_outputs(
                    target_rows, runtime_weights, temp_config, mode, benchmark_genes
                )
                summary = next((row for row in comparison_rows if row["modality"] == modality), None)
                if summary is None:
                    summary = {
                        "method": "ChampionChallenger",
                        "method_type": "ensemble",
                        "modality": modality,
                        "sample_count": len({row["sample"] for row in target_rows}),
                        "gene_rows": len({(row["sample"], row["gene"]) for row in target_rows}),
                        "callable_rate": 0.0,
                        "accuracy_among_callable": 0.0,
                        "overall_correct_call_rate": 0.0,
                        "overall_correct_call_rate_ci_lo": 0.0,
                        "overall_correct_call_rate_ci_hi": 0.0,
                    }
                override_stats = summarize_override_effects(trace_rows, truth_index)
                out.append({
                    "method": "ChampionChallenger",
                    "modality": modality,
                    "champion_A": settings["champion_by_gene"].get("A", ""),
                    "champion_B": settings["champion_by_gene"].get("B", ""),
                    "champion_C": settings["champion_by_gene"].get("C", ""),
                    "min_challenger_support_fraction": min_support,
                    "min_challenger_margin": min_margin,
                    "min_supporting_tools": min_tools,
                    "require_non_ambiguity_override": settings["require_non_ambiguity_override"],
                    "override_count": override_stats["override_count"],
                    "corrective_override_count": override_stats["corrective_override_count"],
                    "harmful_override_count": override_stats["harmful_override_count"],
                    "neutral_override_count": override_stats["neutral_override_count"],
                    "sample_count": summary["sample_count"],
                    "gene_rows": summary["gene_rows"],
                    "callable_rate": summary["callable_rate"],
                    "accuracy_among_callable": summary["accuracy_among_callable"],
                    "overall_correct_call_rate": summary["overall_correct_call_rate"],
                    "overall_correct_call_rate_ci_lo": summary["overall_correct_call_rate_ci_lo"],
                    "overall_correct_call_rate_ci_hi": summary["overall_correct_call_rate_ci_hi"],
                })
    return out


def validate_locus_expert_consensus_settings(settings, rows, benchmark_genes):
    if not settings.get("enabled"):
        return settings
    parent_method = settings.get("parent_method", "weighted_consensus")
    if parent_method not in {"weighted_consensus", "meta_consensus"}:
        raise ValueError("locus_expert_consensus.parent_method must be weighted_consensus or meta_consensus")
    allowed_genes = set(benchmark_genes or [])
    tools_by_modality = defaultdict(set)
    for row in rows:
        tools_by_modality[row["modality"]].add(row["tool"])
    wgs_tools = tools_by_modality.get("wgs", set())
    for panel in settings.get("panel_sets", []):
        subsets = panel["tool_subsets_by_gene"]
        missing_genes = sorted(allowed_genes.difference(subsets), key=gene_sort_key)
        if missing_genes:
            raise ValueError("locus_expert_consensus panel %s is missing genes: %s" % (panel["name"], ",".join(missing_genes)))
        for gene, tools in subsets.items():
            if gene not in allowed_genes:
                raise ValueError("locus_expert_consensus panel %s contains unsupported gene: %s" % (panel["name"], gene))
            if not tools:
                raise ValueError("locus_expert_consensus panel %s defines an empty tool list for gene %s" % (panel["name"], gene))
            for tool in tools:
                if tool not in wgs_tools:
                    raise ValueError("locus_expert_consensus panel %s references unknown WGS tool %s for gene %s" % (panel["name"], tool, gene))
    return settings


def filter_rows_for_gene_subset(rows, tool_subsets_by_gene):
    out = []
    for row in rows:
        allowed = tool_subsets_by_gene.get(row["gene"])
        if allowed and row["tool"] in allowed:
            out.append(row)
    return out


def empty_consensus_call_row(entries, method_name, selected_tools):
    truth_1 = entries[0]["truth_allele1"]
    truth_2 = entries[0]["truth_allele2"]
    ambiguity = evaluate_ambiguity(["", ""], [truth_1, truth_2], 2)
    return {
        "sample": entries[0]["sample"],
        "population": entries[0].get("superpopulation", entries[0].get("population", "unknown")),
        "superpopulation": entries[0].get("superpopulation", "unknown"),
        "modality": entries[0]["modality"],
        "gene": entries[0]["gene"],
        "method": method_name,
        "truth_allele1": truth_1,
        "truth_allele2": truth_2,
        "allele1": "",
        "allele2": "",
        "truth_allele1_3field": ambiguity["truth_3field"][0],
        "truth_allele2_3field": ambiguity["truth_3field"][1],
        "allele1_3field": "",
        "allele2_3field": "",
        "call_status": "no_call",
        "is_callable": "0",
        "is_correct": "0",
        "is_correct_2field": "0",
        "is_correct_3field": "0",
        "is_correct_g_group": ambiguity_match_value(ambiguity["g_group_comparable"], False),
        "is_correct_p_group": ambiguity_match_value(ambiguity["p_group_comparable"], False),
        "is_ambiguity_compatible": "0",
        "is_resolution_compatible": "0",
        "compatibility_grade": "missing",
        "match_grade": ambiguity["match_grade"],
        "imgt_hla_version": entries[0].get("imgt_hla_version", ""),
        "agreeing_tools": "",
        "contributing_tools": 0,
        "support_fraction": "",
        "support_margin": "",
        "total_weight": "",
        "discordance_tag": "no_evidence",
        "selected_tools": ",".join(selected_tools),
    }


def build_locus_expert_trace_rows(base_rows, routed_rows, consensus_rows, runtime_weights, config, panel_name, tool_subsets_by_gene, parent_method, mode):
    grouped = defaultdict(list)
    routed_grouped = defaultdict(list)
    for row in base_rows:
        grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    for row in routed_rows:
        routed_grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    consensus_index = {(row["sample"], row["modality"], row["gene"]): row for row in consensus_rows}
    result_rows = []
    params = meta_consensus_parameters(config)
    thresholds = consensus_thresholds(config)
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        sample, modality, gene = key
        selected_tools = tool_subsets_by_gene.get(gene, [])
        routed_entries = routed_grouped.get(key, [])
        consensus_row = consensus_index.get(key, {})
        top_pair = ("", "")
        runner_up_pair = ("", "")
        top_support_fraction = None
        top_support_margin = None
        contributing_tools = 0
        top_weights = {}
        meta_fields = {
            "support_weight_sum": "",
            "mean_calibrated_probability": "",
            "agreement_bonus": "",
            "gene_strength_bonus": "",
            "ambiguity_bonus": "",
            "fragmentation_penalty": "",
            "meta_score": "",
        }
        if parent_method == "weighted_consensus":
            ranked, total_support = ranked_weighted_candidates(routed_entries, runtime_weights)
            top_pair = ranked[0][0] if ranked else ("", "")
            runner_up_pair = ranked[1][0] if len(ranked) > 1 else ("", "")
            top_support_fraction = round(ranked[0][1]["weight"] / total_support, 4) if ranked and total_support else None
            runner_up_weight = ranked[1][1]["weight"] if len(ranked) > 1 else 0.0
            top_support_margin = round((ranked[0][1]["weight"] - runner_up_weight) / total_support, 4) if ranked and total_support else None
            contributing_tools = sum(len(meta["tools"]) for _pair, meta in ranked)
            top_weights = ranked[0][1]["tool_weights"] if ranked else {}
        else:
            ranked, total_weight = ranked_weighted_candidates(routed_entries, runtime_weights)
            candidate_scores = [score_meta_candidate(pair, meta, routed_entries, runtime_weights, params) for pair, meta in ranked]
            candidate_scores.sort(key=lambda item: (-item["meta_score"], pair_to_string(item["pair"])))
            top = candidate_scores[0] if candidate_scores else None
            runner_up = candidate_scores[1] if len(candidate_scores) > 1 else None
            top_pair = top["pair"] if top else ("", "")
            runner_up_pair = runner_up["pair"] if runner_up else ("", "")
            top_support_fraction = round((top["support_weight_sum"] / total_weight), 4) if (top and total_weight) else None
            runner_up_weight = runner_up["support_weight_sum"] if runner_up else 0.0
            top_support_margin = round((top["support_weight_sum"] - runner_up_weight) / total_weight, 4) if (top and total_weight) else None
            contributing_tools = len({row["tool"] for row in routed_entries if row.get("is_callable") == "1"})
            top_weights = top["tool_weights"] if top else {}
            if top:
                meta_fields = {
                    "support_weight_sum": as_string_number(top["support_weight_sum"]),
                    "mean_calibrated_probability": as_string_number(top["mean_calibrated_probability"]),
                    "agreement_bonus": as_string_number(top["agreement_bonus"]),
                    "gene_strength_bonus": as_string_number(top["gene_strength_bonus"]),
                    "ambiguity_bonus": as_string_number(top["ambiguity_bonus"]),
                    "fragmentation_penalty": as_string_number(top["fragmentation_penalty"]),
                    "meta_score": as_string_number(top["meta_score"]),
                }
        row = {
            "panel_name": panel_name,
            "parent_method": parent_method,
            "sample": sample,
            "population": entries[0].get("superpopulation", entries[0].get("population", "unknown")),
            "modality": modality,
            "gene": gene,
            "selected_tools": ",".join(selected_tools),
            "winning_pair": pair_to_string(top_pair),
            "runner_up_pair": pair_to_string(runner_up_pair),
            "winning_support_fraction": as_string_number(top_support_fraction),
            "winning_support_margin": as_string_number(top_support_margin),
            "contributing_tools": contributing_tools,
            "top_tool_weights_json": json.dumps(top_weights, sort_keys=True, separators=(",", ":")),
            "call_status": consensus_row.get("call_status", classify_consensus_status(top_support_fraction or 0.0, top_support_margin or 0.0, thresholds, bool(routed_entries), False) if routed_entries else "no_call"),
            "is_correct": consensus_row.get("is_correct", "0"),
        }
        row.update(meta_fields)
        result_rows.append(row)
    return result_rows


def build_locus_expert_consensus_rows(rows, runtime_weights, config, panel_name, tool_subsets_by_gene, parent_method, mode):
    routed_rows = filter_rows_for_gene_subset(rows, tool_subsets_by_gene)
    if parent_method == "meta_consensus":
        parent_rows, _parent_traces = build_meta_consensus_rows(routed_rows, runtime_weights, config)
    else:
        parent_rows = build_weighted_consensus_rows(routed_rows, runtime_weights, config)
    parent_index = {(row["sample"], row["modality"], row["gene"]): row for row in parent_rows}
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        gene = key[2]
        selected_tools = tool_subsets_by_gene.get(gene, [])
        row = parent_index.get(key)
        if row is None:
            row = empty_consensus_call_row(entries, "LocusExpertConsensus", selected_tools)
        else:
            row = dict(row)
            row["method"] = "LocusExpertConsensus"
            row["selected_tools"] = ",".join(selected_tools)
        row["population"] = row.get("superpopulation", entries[0].get("superpopulation", entries[0].get("population", "unknown")))
        row["panel_name"] = panel_name
        row["parent_method"] = parent_method
        out.append(row)
    trace_rows = build_locus_expert_trace_rows(rows, routed_rows, out, runtime_weights, config, panel_name, tool_subsets_by_gene, parent_method, mode)
    return out, trace_rows


def build_locus_expert_outputs(rows, runtime_weights, config, mode, benchmark_genes):
    settings = validate_locus_expert_consensus_settings(locus_expert_consensus_settings(config), rows, benchmark_genes)
    if not settings.get("enabled"):
        return [], [], [], [], []
    comparison_rows = []
    per_gene_rows = []
    call_rows = []
    trace_rows = []
    summary_rows = []
    for panel in settings["panel_sets"]:
        panel_calls, panel_traces = build_locus_expert_consensus_rows(
            rows,
            runtime_weights,
            config,
            panel["name"],
            panel["tool_subsets_by_gene"],
            settings["parent_method"],
            mode,
        )
        panel_summary = summarize_consensus_method(panel_calls)
        panel_per_gene = build_method_per_gene([], panel_calls)
        comparison_rows.extend(dict(row, panel_name=panel["name"], parent_method=settings["parent_method"], method_type="ensemble") for row in panel_summary)
        per_gene_rows.extend(dict(row, panel_name=panel["name"], parent_method=settings["parent_method"]) for row in panel_per_gene)
        call_rows.extend(panel_calls)
        trace_rows.extend(panel_traces)
        summary_rows.append({
            "panel_name": panel["name"],
            "parent_method": settings["parent_method"],
            "tool_subsets_by_gene": panel["tool_subsets_by_gene"],
            "overall_metrics": panel_summary,
            "per_gene_metrics": panel_per_gene,
        })
    return call_rows, trace_rows, comparison_rows, per_gene_rows, summary_rows


def build_gated_consensus_outputs(rows, majority_rows, weighted_rows, champion_rows, config):
    settings = validate_consensus_gating_settings(consensus_gating_settings(config), config)
    if not settings.get("enabled"):
        return [], [], [], [], []
    grouped = defaultdict(list)
    majority_index = {(row["sample"], row["modality"], row["gene"]): row for row in majority_rows}
    weighted_index = {(row["sample"], row["modality"], row["gene"]): row for row in weighted_rows}
    champion_index = {(row["sample"], row["modality"], row["gene"]): row for row in champion_rows}
    call_rows = []
    trace_rows = []
    difficulty_counts = {"easy": 0, "hard": 0}
    for row in rows:
        grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        majority_row = majority_index[key]
        weighted_row = weighted_index[key]
        champion_row = champion_index[key]
        features = compute_gating_features(entries, majority_row, weighted_row)
        difficulty_class, reasons = classify_hard_locus(features, settings)
        difficulty_counts[difficulty_class] += 1
        selected_row = champion_row if difficulty_class == "hard" else majority_row
        selected_policy = settings["hard_policy"] if difficulty_class == "hard" else settings["easy_policy"]
        call_rows.append(dict(
            selected_row,
            method="GatedConsensus",
            population=selected_row.get("population", selected_row.get("superpopulation", "unknown")),
            difficulty_class=difficulty_class,
            selected_policy=selected_policy,
        ))
        trace_rows.append({
            "sample": entries[0]["sample"],
            "population": entries[0].get("superpopulation", entries[0].get("population", "unknown")),
            "modality": entries[0]["modality"],
            "gene": entries[0]["gene"],
            "method": "GatedConsensus",
            "difficulty_class": difficulty_class,
            "selected_policy": selected_policy,
            "distinct_pair_count": features["distinct_pair_count"],
            "callable_tool_count": features["callable_tool_count"],
            "weighted_winning_pair": features["weighted_winning_pair"],
            "weighted_runner_up_pair": next(
                (
                    pair
                    for pair in sorted(
                        {
                            pair_to_string((row["allele1"], row["allele2"]))
                            for row in entries
                            if row.get("is_callable") == "1" and pair_to_string((row["allele1"], row["allele2"])) and pair_to_string((row["allele1"], row["allele2"])) != features["weighted_winning_pair"]
                        }
                    )
                ),
                "",
            ),
            "weighted_support_fraction": as_string_number(features["weighted_support_fraction"]),
            "weighted_support_margin": as_string_number(features["weighted_support_margin"]),
            "majority_pair": features["majority_pair"],
            "majority_weighted_disagree": "1" if features["majority_weighted_disagree"] else "0",
            "weighted_top_pair_is_duplicated": "1" if features["weighted_top_pair_is_duplicated"] else "0",
            "alternative_nonduplicated_exists": "1" if features["alternative_nonduplicated_exists"] else "0",
            "hardness_reasons": ",".join(reasons),
            "call_status": selected_row.get("call_status", ""),
            "is_correct": selected_row.get("is_correct", ""),
        })
    comparison_rows = [dict(row, method_type="ensemble") for row in summarize_consensus_method(call_rows)]
    per_gene_rows = build_method_per_gene([], call_rows)
    summary_rows = [{
        "hard_locus_rule": settings["hard_locus_rule"],
        "easy_policy": settings["easy_policy"],
        "hard_policy": settings["hard_policy"],
        "difficulty_counts": difficulty_counts,
        "overall_metrics": comparison_rows,
        "per_gene_metrics": per_gene_rows,
    }]
    return call_rows, trace_rows, comparison_rows, per_gene_rows, summary_rows


def build_ensemble_ablation_outputs(rows, config, mode, weight_alpha, weight_beta):
    settings = ensemble_ablation_settings(config)
    if not settings["enabled"]:
        return [], [], []
    comparison_rows = []
    per_gene_rows = []
    summary_rows = []
    meta_params = meta_consensus_parameters(config)
    for subset in settings["tool_subsets"]:
        subset_rows = filter_rows(rows, tools=subset["tools"])
        if not subset_rows:
            continue
        subset_summary = aggregate_metrics(subset_rows)
        subset_per_gene = build_per_gene_summary(subset_rows)
        subset_tool_weights, subset_gene_weights = build_confidence_weights(
            subset_rows,
            config=config,
            weight_alpha=weight_alpha,
            weight_beta=weight_beta,
            mode=mode,
        )
        subset_runtime_weights = build_runtime_weight_payload(
            subset_tool_weights,
            subset_gene_weights,
            weight_alpha=weight_alpha,
            weight_beta=weight_beta,
            mode=mode,
            meta_params=meta_params,
        )
        subset_majority_rows = build_majority_vote_rows(subset_rows)
        subset_weighted_rows = build_weighted_consensus_rows(subset_rows, subset_runtime_weights, config)
        subset_meta_rows, _subset_meta_trace_rows = build_meta_consensus_rows(subset_rows, subset_runtime_weights, config)
        subset_method_rows = build_method_comparison(subset_summary, subset_majority_rows, subset_weighted_rows)
        subset_method_rows.extend(dict(row, method_type="ensemble") for row in summarize_consensus_method(subset_meta_rows))
        subset_method_per_gene = build_method_per_gene(subset_per_gene, subset_majority_rows, subset_weighted_rows, subset_meta_rows)
        subset_label = ",".join(subset["tools"])
        for row in subset_method_rows:
            comparison_rows.append(dict(row, ablation_name=subset["name"], tool_subset=subset_label))
        for row in subset_method_per_gene:
            per_gene_rows.append(dict(row, ablation_name=subset["name"], tool_subset=subset_label))
        for row in subset_method_rows:
            summary_rows.append({
                "ablation_name": subset["name"],
                "tool_subset": subset["tools"],
                "method": row["method"],
                "method_type": row["method_type"],
                "modality": row["modality"],
                "sample_count": row["sample_count"],
                "gene_rows": row["gene_rows"],
                "callable_rate": row["callable_rate"],
                "accuracy_among_callable": row["accuracy_among_callable"],
                "overall_correct_call_rate": row["overall_correct_call_rate"],
                "overall_correct_call_rate_ci_lo": row.get("overall_correct_call_rate_ci_lo", ""),
                "overall_correct_call_rate_ci_hi": row.get("overall_correct_call_rate_ci_hi", ""),
            })
    return comparison_rows, per_gene_rows, summary_rows


def build_per_gene_gain_table(method_per_gene_rows):
    out = []
    grouped = defaultdict(list)
    for row in method_per_gene_rows:
        grouped[(row["modality"], row["gene"])].append(row)
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][0]), gene_sort_key(item[0][1]))):
        weighted = next((row for row in entries if row["method"] == "WeightedConsensus"), None)
        majority = next((row for row in entries if row["method"] == "MajorityVote"), None)
        single_tool_best = max([row for row in entries if row["method_type"] == "single_tool"], key=lambda row: row["overall_correct_call_rate"], default=None)
        if not weighted:
            continue
        out.append({
            "modality": key[0],
            "gene": key[1],
            "weighted_correct_call_rate": weighted["overall_correct_call_rate"],
            "majority_correct_call_rate": "" if not majority else majority["overall_correct_call_rate"],
            "best_single_tool": "" if not single_tool_best else single_tool_best["method"],
            "best_single_tool_rate": "" if not single_tool_best else single_tool_best["overall_correct_call_rate"],
            "gain_vs_majority": "" if not majority else round(weighted["overall_correct_call_rate"] - majority["overall_correct_call_rate"], 4),
            "gain_vs_best_single": "" if not single_tool_best else round(weighted["overall_correct_call_rate"] - single_tool_best["overall_correct_call_rate"], 4),
        })
    return out


def pair_to_string(pair):
    if not pair or not all(pair):
        return ""
    return "%s+%s" % (pair[0], pair[1])


def build_pairwise_agreement(rows):
    grouped = defaultdict(dict)
    for row in rows:
        if row.get("is_callable") != "1":
            continue
        grouped[(row["modality"], row["gene"], row["sample"])][row["tool"]] = row
    overall = defaultdict(list)
    per_gene = defaultdict(list)
    for key, tool_map in grouped.items():
        modality, gene, _sample = key
        tools = sorted(tool_map)
        for idx, tool_a in enumerate(tools):
            for tool_b in tools[idx + 1:]:
                row_a = tool_map[tool_a]
                row_b = tool_map[tool_b]
                overall[(modality, tool_a, tool_b)].append((row_a, row_b))
                per_gene[(modality, gene, tool_a, tool_b)].append((row_a, row_b))

    def summarize(pairs, include_gene):
        out = []
        for key, entries in sorted(pairs.items(), key=lambda item: (modality_sort_key(item[0][0]), gene_sort_key(item[0][1]) if include_gene else -1, item[0][-2], item[0][-1])):
            exact2 = sum(1 for a, b in entries if pair_to_string((a["allele1"], a["allele2"])) == pair_to_string((b["allele1"], b["allele2"])))
            exact3 = sum(1 for a, b in entries if pair_to_string((a.get("allele1_3field", ""), a.get("allele2_3field", ""))) == pair_to_string((b.get("allele1_3field", ""), b.get("allele2_3field", ""))) and all([a.get("allele1_3field", ""), a.get("allele2_3field", ""), b.get("allele1_3field", ""), b.get("allele2_3field", "")]))
            g_group = sum(1 for a, b in entries if a.get("is_correct_g_group", "") != "" and b.get("is_correct_g_group", "") != "" and pair_to_string((a["allele1"], a["allele2"])) == pair_to_string((b["allele1"], b["allele2"])))
            p_group = sum(1 for a, b in entries if a.get("is_correct_p_group", "") != "" and b.get("is_correct_p_group", "") != "" and pair_to_string((a["allele1"], a["allele2"])) == pair_to_string((b["allele1"], b["allele2"])))
            record = {
                "modality": key[0],
                "tool_a": key[-2],
                "tool_b": key[-1],
                "n_shared_callable_rows": len(entries),
                "exact_2field_agreement_rate": ratio(exact2, len(entries)),
                "exact_3field_agreement_rate": ratio(exact3, len(entries)),
                "g_group_agreement_rate": ratio(g_group, len(entries)),
                "p_group_agreement_rate": ratio(p_group, len(entries)),
            }
            if include_gene:
                record["gene"] = key[1]
            out.append(record)
        return out

    return summarize(overall, include_gene=False), summarize(per_gene, include_gene=True)


def classify_truth_error(row):
    if row.get("is_callable") != "1":
        return "abstained"
    if row.get("compatibility_grade") == "ambiguity_compatible":
        return "ambiguity_compatible"
    if row.get("is_correct") == "1":
        return "exact_match"
    truth_set = {row.get("truth_allele1", ""), row.get("truth_allele2", "")}
    if row.get("allele1", "") in truth_set or row.get("allele2", "") in truth_set:
        return "one_allele_correct"
    def field1(value):
        return value.split(":")[0] if ":" in value else value
    truth_field1 = {field1(v) for v in truth_set if v}
    call_field1 = [field1(v) for v in (row.get("allele1", ""), row.get("allele2", "")) if v]
    if any(value in truth_field1 for value in call_field1):
        return "field1_right_field2_wrong"
    return "wrong_allele_group"


def build_tool_vs_truth_error_taxonomy(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["modality"], row["gene"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        counts = defaultdict(int)
        for row in entries:
            counts[classify_truth_error(row)] += 1
        total = len(entries)
        for category in ["exact_match", "one_allele_correct", "field1_right_field2_wrong", "wrong_allele_group", "abstained", "ambiguity_compatible"]:
            out.append({
                "tool": key[0],
                "modality": key[1],
                "gene": key[2],
                "error_category": category,
                "count": counts.get(category, 0),
                "fraction": ratio(counts.get(category, 0), total),
            })
    return out


def build_tool_disagreement_events(rows, majority_rows, weighted_rows):
    majority_index = {(row["sample"], row["modality"], row["gene"]): row for row in majority_rows}
    weighted_index = {(row["sample"], row["modality"], row["gene"]): row for row in weighted_rows}
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        callable_rows = [row for row in entries if row.get("is_callable") == "1"]
        distinct_pairs = sorted({pair_to_string((row["allele1"], row["allele2"])) for row in callable_rows if pair_to_string((row["allele1"], row["allele2"]))})
        if len(distinct_pairs) < 2:
            continue
        sample, modality, gene = key
        tool_calls = []
        for row in sorted(callable_rows, key=lambda entry: entry["tool"]):
            tool_calls.append({
                "tool": row["tool"],
                "pair": pair_to_string((row["allele1"], row["allele2"])),
                "compatibility_grade": row.get("compatibility_grade", ""),
                "calibrated_probability": coerce_float(row.get("calibrated_probability", "")),
            })
        majority_row = majority_index.get(key, {})
        weighted_row = weighted_index.get(key, {})
        disagreement_class = "majority_weighted_disagree" if pair_to_string((majority_row.get("allele1", ""), majority_row.get("allele2", ""))) != pair_to_string((weighted_row.get("allele1", ""), weighted_row.get("allele2", ""))) else "tool_only_disagreement"
        out.append({
            "sample": sample,
            "population": entries[0].get("superpopulation", entries[0].get("population", "unknown")),
            "modality": modality,
            "gene": gene,
            "truth_allele1": entries[0]["truth_allele1"],
            "truth_allele2": entries[0]["truth_allele2"],
            "callable_tool_count": len(callable_rows),
            "distinct_pair_count": len(distinct_pairs),
            "tool_calls_json": json.dumps(tool_calls, sort_keys=True, separators=(",", ":")),
            "majority_pair": pair_to_string((majority_row.get("allele1", ""), majority_row.get("allele2", ""))),
            "weighted_pair": pair_to_string((weighted_row.get("allele1", ""), weighted_row.get("allele2", ""))),
            "majority_correct": majority_row.get("is_correct", ""),
            "weighted_correct": weighted_row.get("is_correct", ""),
            "disagreement_class": disagreement_class,
        })
    return out


def build_consensus_decision_trace(rows, majority_rows, weighted_rows, weight_payload, mode):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    result_rows = []
    existing_rows = {"MajorityVote": majority_rows, "WeightedConsensus": weighted_rows}
    for method_name, consensus_rows in existing_rows.items():
        consensus_index = {(row["sample"], row["modality"], row["gene"]): row for row in consensus_rows}
        for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
            sample, modality, gene = key
            if method_name == "MajorityVote":
                pair_counts = defaultdict(list)
                for row in entries:
                    pair = allele_pair(row)
                    if all(pair):
                        pair_counts[pair].append(row)
                ranked = sorted(pair_counts.items(), key=lambda item: (-len(item[1]), "%s|%s" % item[0]))
                top_pair = ranked[0][0] if ranked else ("", "")
                runner_up_pair = ranked[1][0] if len(ranked) > 1 else ("", "")
                contributing_tools = sum(len(v) for v in pair_counts.values())
                total_support = float(contributing_tools) if contributing_tools else 0.0
                top_support_fraction = round(len(ranked[0][1]) / total_support, 4) if ranked and total_support else 0.0
                runner_up_n = len(ranked[1][1]) if len(ranked) > 1 else 0
                top_support_margin = round((len(ranked[0][1]) - runner_up_n) / total_support, 4) if ranked and total_support else 0.0
                top_weights = {row["tool"]: 1.0 for row in ranked[0][1]} if ranked else {}
            else:
                ranked, total_support = ranked_weighted_candidates(entries, weight_payload)
                top_pair = ranked[0][0] if ranked else ("", "")
                runner_up_pair = ranked[1][0] if len(ranked) > 1 else ("", "")
                top_support_fraction = round(ranked[0][1]["weight"] / total_support, 4) if ranked and total_support else 0.0
                runner_up_weight = ranked[1][1]["weight"] if len(ranked) > 1 else 0.0
                top_support_margin = round((ranked[0][1]["weight"] - runner_up_weight) / total_support, 4) if ranked and total_support else 0.0
                contributing_tools = sum(len(meta["tools"]) for _pair, meta in ranked)
                top_weights = ranked[0][1]["tool_weights"] if ranked else {}
            consensus_row = consensus_index.get(key, {})
            result_rows.append({
                "sample": sample,
                "population": entries[0].get("superpopulation", entries[0].get("population", "unknown")),
                "modality": modality,
                "gene": gene,
                "benchmark_mode": mode,
                "method": method_name,
                "winning_pair": pair_to_string(top_pair),
                "runner_up_pair": pair_to_string(runner_up_pair),
                "winning_support_fraction": as_string_number(top_support_fraction),
                "winning_support_margin": as_string_number(top_support_margin),
                "contributing_tools": contributing_tools,
                "top_tool_weights_json": json.dumps(top_weights, sort_keys=True, separators=(",", ":")),
                "call_status": consensus_row.get("call_status", ""),
                "is_correct": consensus_row.get("is_correct", ""),
            })
    return result_rows


def build_population_method_comparison(method_rows, mode):
    out = []
    grouped = defaultdict(list)
    for row in method_rows:
        grouped[(row.get("superpopulation", "unknown"), row["method"], row["modality"])].append(row)
    for key, entries in sorted(grouped.items(), key=lambda item: (item[0][0], modality_sort_key(item[0][2]), item[0][1])):
        total = len(entries)
        callable_n = sum(int(row["is_callable"]) for row in entries)
        correct_n = sum(int(row["is_correct"]) for row in entries)
        out.append({
            "superpopulation": key[0],
            "method": key[1],
            "modality": key[2],
            "benchmark_mode": mode,
            "sample_count": len({row["sample"] for row in entries}),
            "gene_rows": total,
            "callable_rate": ratio(callable_n, total),
            "accuracy_among_callable": ratio(correct_n, callable_n),
            "overall_correct_call_rate": ratio(correct_n, total),
        })
    return out


def build_locus_difficulty_summary(rows, method_per_gene_rows):
    out = []
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["modality"], row["gene"])].append(row)
    method_grouped = defaultdict(list)
    for row in method_per_gene_rows:
        method_grouped[(row["modality"], row["gene"])].append(row)
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][0]), gene_sort_key(item[0][1]))):
        tools = sorted({row["tool"] for row in entries})
        per_tool_accuracy = []
        for tool in tools:
            tool_rows = [row for row in entries if row["tool"] == tool]
            per_tool_accuracy.append(ratio(sum(int(row["is_correct"]) for row in tool_rows), len(tool_rows)))
        disagreement_events = 0
        sample_grouped = defaultdict(list)
        for row in entries:
            sample_grouped[row["sample"]].append(row)
        for sample_entries in sample_grouped.values():
            callable_pairs = {pair_to_string((row["allele1"], row["allele2"])) for row in sample_entries if row.get("is_callable") == "1" and pair_to_string((row["allele1"], row["allele2"]))}
            if len(callable_pairs) > 1:
                disagreement_events += 1
        method_rows = method_grouped.get(key, [])
        weighted = next((row for row in method_rows if row["method"] == "WeightedConsensus"), None)
        best_single = max((row for row in method_rows if row["method_type"] == "single_tool"), key=lambda row: row["overall_correct_call_rate"], default=None)
        out.append({
            "modality": key[0],
            "gene": key[1],
            "n_tools": len(tools),
            "best_single_tool_rate": "" if not best_single else best_single["overall_correct_call_rate"],
            "tool_accuracy_spread": round(max(per_tool_accuracy) - min(per_tool_accuracy), 4) if per_tool_accuracy else "",
            "disagreement_frequency": ratio(disagreement_events, len(sample_grouped)),
            "weighted_vs_best_single_gain": "" if not (weighted and best_single) else round(weighted["overall_correct_call_rate"] - best_single["overall_correct_call_rate"], 4),
        })
    return out


def summarize_multiresolution_rows(rows, name_key, include_gene=False):
    grouped = defaultdict(list)
    for row in rows:
        gene = row.get("gene", "") if include_gene else ""
        grouped[(row[name_key], row["modality"], gene)].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), gene_sort_key(item[0][2]), item[0][0])):
        total = len(entries)
        callable_n = sum(int(row["is_callable"]) for row in entries)
        g_rows = [row for row in entries if row.get("is_correct_g_group", "") != ""]
        p_rows = [row for row in entries if row.get("is_correct_p_group", "") != ""]
        record = {
            name_key: key[0],
            "modality": key[1],
            "sample_count": len({row["sample"] for row in entries}),
            "gene_rows": total,
            "callable_rate": ratio(callable_n, total),
            "exact_2field_rate": ratio(sum(int(row["is_correct_2field"]) for row in entries), total),
            "exact_3field_rate": ratio(sum(int(row["is_correct_3field"]) for row in entries), total),
            "g_group_match_rate": ratio(sum(int(row["is_correct_g_group"]) for row in g_rows), len(g_rows)),
            "p_group_match_rate": ratio(sum(int(row["is_correct_p_group"]) for row in p_rows), len(p_rows)),
            "ambiguity_compatible_rate": ratio(sum(int(row.get("is_ambiguity_compatible", "0")) for row in entries), total),
        }
        if include_gene:
            record["gene"] = key[2]
        out.append(record)
    return out


def build_population_resolution_summary(rows, mode):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row.get("superpopulation", "unknown"), row["tool"], row["modality"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (item[0][0], modality_sort_key(item[0][2]), item[0][1])):
        total = len(entries)
        out.append({
            "superpopulation": key[0],
            "tool": key[1],
            "modality": key[2],
            "benchmark_mode": mode,
            "gene_rows": total,
            "exact_2field_rate": ratio(sum(int(row["is_correct_2field"]) for row in entries), total),
            "exact_3field_rate": ratio(sum(int(row["is_correct_3field"]) for row in entries), total),
            "ambiguity_compatible_rate": ratio(sum(int(row.get("is_ambiguity_compatible", "0")) for row in entries), total),
            "resolution_compatible_rate": ratio(sum(int(row.get("is_resolution_compatible", "0")) for row in entries), total),
        })
    return out


def build_population_conflict_summary(rows, mode):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row.get("superpopulation", "unknown"), row["tool"], row["modality"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (item[0][0], modality_sort_key(item[0][2]), item[0][1])):
        counts = defaultdict(int)
        for row in entries:
            counts[classify_truth_error(row)] += 1
        total = len(entries)
        for category in ["exact_match", "ambiguity_compatible", "one_allele_correct", "field1_right_field2_wrong", "wrong_allele_group", "abstained"]:
            out.append({
                "superpopulation": key[0],
                "tool": key[1],
                "modality": key[2],
                "benchmark_mode": mode,
                "conflict_category": category,
                "count": counts.get(category, 0),
                "fraction": ratio(counts.get(category, 0), total),
            })
    return out


def calibration_bin_label(value, bin_count=5):
    idx = min(bin_count - 1, int(value * bin_count))
    lower = idx / float(bin_count)
    upper = (idx + 1) / float(bin_count)
    return idx + 1, round(lower, 2), round(upper, 2)


def build_confidence_bin_summary(rows, bin_count=5):
    grouped = defaultdict(list)
    for row in rows:
        score = coerce_float(row.get("confidence_score", ""))
        if score is None:
            continue
        grouped[(row["tool"], row["modality"], calibration_bin_label(score, bin_count)[0])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], item[0][2])):
        bin_idx, lower, upper = calibration_bin_label(coerce_float(entries[0]["confidence_score"]), bin_count)
        out.append({
            "tool": key[0],
            "modality": key[1],
            "bin_index": bin_idx,
            "bin_lower": lower,
            "bin_upper": upper,
            "n_rows": len(entries),
            "mean_confidence": round(statistics.mean(coerce_float(row["confidence_score"]) for row in entries), 4),
            "observed_accuracy": round(statistics.mean(int(row["is_correct"]) for row in entries), 4),
        })
    return out


def build_confidence_calibration_summary(rows, bin_count=5, mode=DEFAULT_BENCHMARK_MODE):
    grouped = defaultdict(list)
    for row in rows:
        raw_score = coerce_float(row.get("confidence_score", ""))
        calibrated_score = coerce_float(row.get("calibrated_probability", ""))
        if raw_score is None and calibrated_score is None:
            continue
        grouped[(row["tool"], row["modality"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0])):
        stats = calibration_stats(entries, bin_count=bin_count)
        calibrated_stats = calibration_stats(entries, bin_count=bin_count, score_field="calibrated_probability")
        out.append({
            "tool": key[0],
            "modality": key[1],
            "benchmark_mode": mode,
            "calibration_method": next((clean_token(row.get("calibration_method", "")) for row in entries if clean_token(row.get("calibration_method", ""))), ""),
            "n_rows": stats["n_rows"],
            "mean_confidence": stats["mean_confidence"],
            "observed_accuracy": stats["observed_accuracy"],
            "brier_score": stats["brier_score"],
            "expected_calibration_error": stats["expected_calibration_error"],
            "mean_calibrated_probability": calibrated_stats["mean_confidence"],
            "calibrated_brier_score": calibrated_stats["brier_score"],
            "calibrated_expected_calibration_error": calibrated_stats["expected_calibration_error"],
        })
    return out


def confidence_error_record(entries):
    total = len(entries)
    error_count = sum(1 for row in entries if row.get("is_correct") == "0")
    mean_conf = statistics.mean(coerce_float(row["confidence_score"]) for row in entries) if entries else 0.0
    return {
        "n_rows": total,
        "error_count": error_count,
        "error_rate": ratio(error_count, total),
        "observed_accuracy": ratio(total - error_count, total),
        "mean_confidence": round(mean_conf, 4) if entries else "",
    }


def build_confidence_error_summary(rows, bin_count=5):
    grouped = defaultdict(list)
    gene_grouped = defaultdict(list)
    for row in rows:
        score = coerce_float(row.get("confidence_score", ""))
        if score is None:
            continue
        bin_index, lower, upper = calibration_bin_label(score, bin_count)
        grouped[(row["tool"], row["modality"], bin_index, lower, upper)].append(row)
        gene_grouped[(row["tool"], row["modality"], row["gene"], bin_index, lower, upper)].append(row)
    summary_rows = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], item[0][2])):
        record = confidence_error_record(entries)
        record.update({
            "tool": key[0],
            "modality": key[1],
            "bin_index": key[2],
            "bin_lower": key[3],
            "bin_upper": key[4],
        })
        summary_rows.append(record)
    per_gene_rows = []
    for key, entries in sorted(gene_grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]), item[0][3])):
        record = confidence_error_record(entries)
        record.update({
            "tool": key[0],
            "modality": key[1],
            "gene": key[2],
            "bin_index": key[3],
            "bin_lower": key[4],
            "bin_upper": key[5],
        })
        per_gene_rows.append(record)
    return summary_rows, per_gene_rows


def build_abstention_tradeoff(rows, weight_payload, config):
    base = consensus_thresholds(config)
    tradeoffs = []
    for min_support in [0.45, 0.55, 0.65, 0.75, 0.85]:
        local_cfg = {"benchmark": {"consensus": {"min_support": min_support, "min_margin": base["min_margin"]}}}
        consensus_rows = build_weighted_consensus_rows(rows, weight_payload, local_cfg)
        total = len(consensus_rows)
        called = [row for row in consensus_rows if row["call_status"] == "called"]
        correct_called = sum(int(row["is_correct"]) for row in called)
        tradeoffs.append({
            "min_support": min_support,
            "min_margin": base["min_margin"],
            "call_rate": ratio(len(called), total),
            "no_call_rate": ratio(sum(1 for row in consensus_rows if row["call_status"] == "no_call"), total),
            "low_confidence_rate": ratio(sum(1 for row in consensus_rows if row["call_status"] == "low_confidence"), total),
            "accuracy_among_called": ratio(correct_called, len(called)),
            "overall_correct_call_rate": ratio(correct_called, total),
        })
    return tradeoffs


def _modality_mean_confidence(rows):
    """Return mean confidence_score across callable tools in a pre-filtered row list."""
    values = [
        coerce_float(r.get("confidence_score", ""))
        for r in rows
        if r.get("is_callable") == "1" and coerce_float(r.get("confidence_score", "")) is not None
    ]
    return round(sum(values) / len(values), 4) if values else None


def _modality_mean_read_support(rows):
    """Return mean read_support across callable tools in a pre-filtered row list."""
    values = [
        coerce_float(r.get("read_support", ""))
        for r in rows
        if r.get("is_callable") == "1" and coerce_float(r.get("read_support", "")) is not None
    ]
    return round(sum(values) / len(values), 4) if values else None


def arbitrate_dna_rna_discordance(dna_summary, rna_summary):
    """Apply biological priority rules to resolve a WES/WGS vs RNA-seq discordance.

    Parameters
    ----------
    dna_summary : dict
        Keys: pair (tuple|None), callable_tools (int), support_fraction (float),
        support_margin (float), mean_confidence (float|None), mean_read_support (float|None).
    rna_summary : dict
        Same keys as dna_summary but for the rnaseq modality.

    Returns
    -------
    dict with keys:
        arbitrated_pair     : tuple|None  — resolved allele pair, or None if abstained
        arbitration_rule    : str         — identifier of the winning rule
        arbitration_outcome : str         — "dna_wins" | "rna_wins" | "abstain"
    """
    rna_callable = rna_summary.get("callable_tools", 0)
    rna_frac = rna_summary.get("support_fraction") or 0.0
    rna_reads = rna_summary.get("mean_read_support")
    dna_callable = dna_summary.get("callable_tools", 0)
    dna_frac = dna_summary.get("support_fraction") or 0.0
    dna_pair = dna_summary.get("pair")
    rna_pair = rna_summary.get("pair")

    # Rule 1 — RNA has no callable tools
    if rna_callable == 0:
        return {"arbitrated_pair": dna_pair, "arbitration_rule": "rna_no_callable_tools", "arbitration_outcome": "dna_wins"}

    # Rule 2 — RNA read support is below reliability floor (skip if not reported)
    if rna_reads is not None and rna_reads < RNA_MIN_READ_SUPPORT:
        return {"arbitrated_pair": dna_pair, "arbitration_rule": "rna_low_read_support", "arbitration_outcome": "dna_wins"}

    # Rule 3 — DNA evidence is thin but RNA consensus is confident
    if dna_callable <= 1 and rna_frac >= 0.65:
        return {"arbitrated_pair": rna_pair, "arbitration_rule": "dna_insufficient_evidence", "arbitration_outcome": "rna_wins"}

    # Rule 4 — DNA consensus is strong, RNA consensus is weak
    if dna_frac >= 0.70 and rna_frac < 0.55:
        return {"arbitrated_pair": dna_pair, "arbitration_rule": "dna_stronger_consensus", "arbitration_outcome": "dna_wins"}

    # Rule 5 — RNA consensus is strong, DNA consensus is weak
    if rna_frac >= 0.70 and dna_frac < 0.55:
        return {"arbitrated_pair": rna_pair, "arbitration_rule": "rna_stronger_consensus", "arbitration_outcome": "rna_wins"}

    # Rule 6 — No rule fired; genuine ambiguity
    return {"arbitrated_pair": None, "arbitration_rule": "high_confidence_conflict", "arbitration_outcome": "abstain"}


def build_discordance_rows(harmonized_rows, weighted_rows):
    weighted_index = {(row["sample"], row["modality"], row["gene"]): row for row in weighted_rows}
    grouped = defaultdict(list)
    for row in harmonized_rows:
        grouped[(row["sample"], row["gene"])].append(row)
    # Pre-index by (sample, gene, modality) so signal-extraction helpers are O(1)
    hrx = defaultdict(list)
    for row in harmonized_rows:
        hrx[(row["sample"], row["gene"], row["modality"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (item[0][0], gene_sort_key(item[0][1]))):
        modality_pairs = {}
        for modality in sorted({row["modality"] for row in entries}, key=modality_sort_key):
            modality_entries = [row for row in entries if row["modality"] == modality]
            pair, agreeing_tools, callable_tools, support_fraction, support_margin, is_tie = majority_vote_for_group(modality_entries)
            # For discordance detection, treat ties as absent — a contested call should not
            # mask genuine cross-modality discordance or missing-RNA signals.
            modality_pairs[modality] = {"pair": None if is_tie else pair, "callable_tools": callable_tools, "support_fraction": support_fraction, "support_margin": support_margin}
            weighted_row = weighted_index.get((key[0], modality, key[1]))
            if weighted_row and weighted_row["discordance_tag"] in {"low_evidence_conflict", "technical_conflict", "no_evidence"}:
                out.append({
                    "sample": key[0],
                    "gene": key[1],
                    "scope": modality,
                    "tag": weighted_row["discordance_tag"],
                    "detail": weighted_row["call_status"],
                })
        dna_pairs = {modality_pairs[m]["pair"] for m in modality_pairs if m in {"wes", "wgs"} and modality_pairs[m]["pair"]}
        rna_pair = modality_pairs.get("rnaseq", {}).get("pair")
        if len(dna_pairs) > 1:
            out.append({"sample": key[0], "gene": key[1], "scope": "dna", "tag": "technical_conflict", "detail": "wes_wgs_disagree"})
        if rna_pair and dna_pairs and rna_pair not in dna_pairs:
            # Select the strongest DNA modality (most callable tools; prefer wes on tie)
            dna_mods = [m for m in modality_pairs if m in {"wes", "wgs"}]
            best_dna_mod = max(dna_mods, key=lambda m: (modality_pairs[m]["callable_tools"], m == "wes")) if dna_mods else None
            if best_dna_mod:
                dna_sum = dict(modality_pairs[best_dna_mod])
                dna_sum["mean_confidence"] = _modality_mean_confidence(hrx[(key[0], key[1], best_dna_mod)])
                dna_sum["mean_read_support"] = None  # DNA tools do not emit read_support
            else:
                dna_sum = {"pair": None, "callable_tools": 0, "support_fraction": 0.0, "support_margin": 0.0, "mean_confidence": None, "mean_read_support": None}
            rna_sum = dict(modality_pairs.get("rnaseq", {"pair": rna_pair, "callable_tools": 0, "support_fraction": 0.0, "support_margin": 0.0}))
            rna_sum["mean_confidence"] = _modality_mean_confidence(hrx[(key[0], key[1], "rnaseq")])
            rna_sum["mean_read_support"] = _modality_mean_read_support(hrx[(key[0], key[1], "rnaseq")])
            arb = arbitrate_dna_rna_discordance(dna_sum, rna_sum)
            arb_pair = arb["arbitrated_pair"]
            out.append({
                "sample": key[0],
                "gene": key[1],
                "scope": "cross_modality",
                "tag": "dna_rna_discordance",
                "detail": "rna_pair_differs_from_dna",
                "arbitrated_pair": "%s+%s" % arb_pair if arb_pair else "",
                "arbitration_rule": arb["arbitration_rule"],
                "arbitration_outcome": arb["arbitration_outcome"],
                "dna_callable_tools": dna_sum.get("callable_tools", ""),
                "dna_support_fraction": as_string_number(dna_sum.get("support_fraction")),
                "rna_callable_tools": rna_sum.get("callable_tools", ""),
                "rna_support_fraction": as_string_number(rna_sum.get("support_fraction")),
                "rna_mean_read_support": as_string_number(rna_sum.get("mean_read_support")),
            })
        if not rna_pair and dna_pairs:
            out.append({"sample": key[0], "gene": key[1], "scope": "rnaseq", "tag": "possible_expression_bias", "detail": "rna_missing_dna_present"})
    return out


def summarize_discordance(rows):
    grouped = defaultdict(int)
    for row in rows:
        grouped[(row["scope"], row["tag"])] += 1
    out = []
    for key, count in sorted(grouped.items()):
        out.append({"scope": key[0], "tag": key[1], "n_events": count})
    return out


def build_sample_discordance_flags(discordance_rows):
    """Aggregate per-sample cross-modality discordance events and emit an integrity flag.

    Only samples with at least one dna_rna_discordance event are included.
    Samples with only possible_expression_bias events are excluded — missing
    RNA calls are normal quality variability, not active cross-modality conflict.
    """
    discordant_by_sample = defaultdict(list)
    bias_by_sample = defaultdict(list)
    for row in discordance_rows:
        if row["tag"] == "dna_rna_discordance":
            discordant_by_sample[row["sample"]].append(row["gene"])
        elif row["tag"] == "possible_expression_bias":
            bias_by_sample[row["sample"]].append(row["gene"])

    out = []
    for sample in sorted(discordant_by_sample):
        discordant_genes = sorted(discordant_by_sample[sample], key=gene_sort_key)
        bias_genes = sorted(bias_by_sample.get(sample, []), key=gene_sort_key)
        n_discordant = len(discordant_genes)
        if n_discordant >= 3:
            integrity_flag = "critical"
            flag_reason = (
                f"{n_discordant}-locus WES/RNA conflict ({','.join(discordant_genes)}); "
                "simultaneous multi-locus discordance indicates likely sample integrity problem "
                "(swap, contamination, or extraction failure)"
            )
        elif n_discordant == 2:
            integrity_flag = "warning"
            flag_reason = (
                f"2-locus WES/RNA conflict ({','.join(discordant_genes)}); "
                "elevated suspicion of sample integrity issue"
            )
        else:
            integrity_flag = "nominal"
            flag_reason = (
                f"1-locus WES/RNA conflict ({','.join(discordant_genes)}); "
                "within background tool error rate"
            )
        out.append({
            "sample": sample,
            "n_discordant_loci": n_discordant,
            "discordant_loci": ",".join(discordant_genes),
            "n_expression_bias_loci": len(bias_genes),
            "expression_bias_loci": ",".join(bias_genes),
            "integrity_flag": integrity_flag,
            "flag_reason": flag_reason,
        })
    return out


def build_arbitration_accuracy_table(discordance_rows, harmonized_rows):
    """Evaluate arbitration accuracy vs ground truth for all dna_rna_discordance rows.

    Parameters
    ----------
    discordance_rows : list[dict]
        Output of build_discordance_rows().
    harmonized_rows : list[dict]
        Full harmonized benchmark rows (contains truth_allele1, truth_allele2).

    Returns
    -------
    list[dict] — one row per dna_rna_discordance event with correctness annotation.
    """
    truth_index = {}
    for row in harmonized_rows:
        k = (row["sample"], row["gene"])
        if k not in truth_index:
            truth_index[k] = (row["truth_allele1"], row["truth_allele2"])

    out = []
    for row in discordance_rows:
        if row.get("tag") != "dna_rna_discordance":
            continue
        sample, gene = row["sample"], row["gene"]
        truth_pair = truth_index.get((sample, gene))
        arb_pair_str = row.get("arbitrated_pair", "")
        outcome = row.get("arbitration_outcome", "")

        if not arb_pair_str or outcome == "abstain":
            is_correct = ""
        else:
            arb_alleles = arb_pair_str.split("+", 1)
            truth_list = list(truth_pair) if truth_pair else ["", ""]
            ambiguity = evaluate_ambiguity(arb_alleles, truth_list, 2)
            is_correct = "1" if ambiguity["match_2field"] else "0"

        out.append({
            "sample": sample,
            "gene": gene,
            "arbitration_rule": row.get("arbitration_rule", ""),
            "arbitration_outcome": outcome,
            "arbitrated_pair": arb_pair_str,
            "truth_allele1": truth_pair[0] if truth_pair else "",
            "truth_allele2": truth_pair[1] if truth_pair else "",
            "is_correct": is_correct,
        })
    return out


def summarize_arbitration_rules(arbitration_accuracy_rows):
    """Aggregate per-rule counts: n_resolved, n_correct, n_abstain, accuracy."""
    stats = defaultdict(lambda: {"n_resolved": 0, "n_correct": 0, "n_abstain": 0})
    for row in arbitration_accuracy_rows:
        rule = row["arbitration_rule"]
        if row["arbitration_outcome"] == "abstain":
            stats[rule]["n_abstain"] += 1
        else:
            stats[rule]["n_resolved"] += 1
            if row["is_correct"] == "1":
                stats[rule]["n_correct"] += 1
    out = []
    for rule, s in sorted(stats.items()):
        n = s["n_resolved"]
        out.append({
            "arbitration_rule": rule,
            "n_resolved": n,
            "n_correct": s["n_correct"],
            "n_abstain": s["n_abstain"],
            "accuracy_among_resolved": as_string_number(ratio(s["n_correct"], n)) if n > 0 else "",
        })
    return out


def build_metadata(config, shared_genes, harmonized_rows, tool_weights, weight_alpha=0.7, weight_beta=0.3, mode=DEFAULT_BENCHMARK_MODE):
    parser_coverage = {}
    for run in config.get("runs", []):
        parser_coverage.setdefault(run.get("tool", ""), {})[clean_token(run.get("modality", "")).lower()] = {
            "confidence_parser": run.get("confidence_parser", ""),
            "confidence_glob": run.get("confidence_glob", ""),
        }
    benchmark_cfg = config.get("benchmark", {})
    imgt_hla_version = resolve_imgt_hla_version(config)
    return {
        "truth_source": config.get("truth", {}).get("source", config.get("truth", {}).get("path", "")),
        "truth_path": config.get("truth", {}).get("path", ""),
        "truth_panel_name": config.get("truth", {}).get("panel_name", "1000G public HLA truth"),
        "truth_panel_date": clean_token(config.get("truth", {}).get("acquisition_date", "")),
        "truth_resolution_limit": config.get("truth", {}).get("resolution_limit", "2-field"),
        "resolution_asymmetry_warning": config.get("truth", {}).get("resolution_limit", "2-field") == "2-field",
        "truth_hierarchy": config.get("truth", {}).get("hierarchy", ["orthogonal_clinical_typing", "targeted_hla_ngs", "public_reference_truth"]),
        "resolution": benchmark_cfg.get("resolution", 2),
        "imgt_hla_version": imgt_hla_version,
        "success_definition": "correct call only at configured primary resolution",
        "shared_genes_main": list(shared_genes),
        "modalities": sorted({clean_token(run.get("modality", "")).lower() for run in config.get("runs", [])}, key=modality_sort_key),
        "tools": [run.get("tool", "") for run in config.get("runs", [])],
        "ambiguity_evaluation": {
            "primary_resolution": benchmark_cfg.get("resolution", 2),
            "secondary_resolutions": benchmark_cfg.get("secondary_resolutions", [2, 3]),
            "group_levels": ["G", "P"],
            "notes": "G-group and P-group rates are only reported for rows where both truth and call include comparable group-suffixed alleles.",
        },
        "configured_benchmark_genes": resolve_benchmark_genes(config),
        "confidence_weighting": {
            "benchmark_mode": mode,
            "weight_version": WEIGHT_VERSION,
            "default_target_reads": DEFAULT_TARGET_READS,
            "formula": f"{weight_alpha} * overall_correct_call_rate + {weight_beta} * effective_confidence_score",
            "effective_confidence_score": "base_reliability + shrink_factor * (mean_confidence_score - base_reliability)",
            "fallback": "overall_correct_call_rate when confidence is missing or blocked by guardrail",
            "guardrail": confidence_guardrail_settings(config),
            "probabilistic_calibration": probabilistic_calibration_settings(config),
            "bayesian_shrinkage": bayesian_shrinkage_settings(config),
            "parser_coverage": parser_coverage,
            "tool_weight_count": len(tool_weights),
            "harmonized_rows_with_confidence": sum(1 for row in harmonized_rows if coerce_float(row.get("confidence_score", "")) is not None),
            "harmonized_rows_with_calibrated_probability": sum(1 for row in harmonized_rows if coerce_float(row.get("calibrated_probability", "")) is not None),
        },
    }


def save_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def svg_header(width, height):
    return [
        '<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s">' % (width, height, width, height),
        '<style>',
        'text { font-family: Arial, Helvetica, sans-serif; fill: #1f2933; }',
        '.title { font-size: 20px; font-weight: bold; }',
        '.label { font-size: 12px; }',
        '.small { font-size: 10px; }',
        '.axis { stroke: #36454F; stroke-width: 1; }',
        '.grid { stroke: #d9e2ec; stroke-width: 1; }',
        '</style>',
    ]


def write_svg(path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(lines + ["</svg>"]))


def draw_method_accuracy_chart(path, title, rows):
    width, height = 1200, 560
    margin_left, margin_top, margin_bottom = 110, 80, 120
    chart_width = width - margin_left - 40
    chart_height = height - margin_top - margin_bottom
    methods = [row["method"] for row in rows]
    modalities = ordered_modalities(rows)
    lines = svg_header(width, height)
    lines.append('<text x="%s" y="32" class="title">%s</text>' % (margin_left, title))
    for tick in range(6):
        y = margin_top + chart_height - (chart_height * tick / 5.0)
        lines.append('<line x1="%s" y1="%s" x2="%s" y2="%s" class="grid" />' % (margin_left, y, margin_left + chart_width, y))
    group_width = chart_width / float(max(len(methods), 1))
    bar_width = min(24, group_width / float(max(len(modalities), 1) + 1))
    for idx, row in enumerate(rows):
        center_x = margin_left + group_width * idx + group_width / 2.0
        lines.append('<text x="%s" y="%s" text-anchor="middle" class="small">%s</text>' % (center_x, margin_top + chart_height + 26, row["method"]))
        value = float(row["overall_correct_call_rate"])
        modality = row["modality"]
        x = center_x - bar_width / 2.0
        bar_height = chart_height * value
        y = margin_top + chart_height - bar_height
        lines.append('<rect x="%s" y="%s" width="%s" height="%s" fill="%s" rx="4" />' % (x, y, bar_width, bar_height, modality_style(modality)["fill"]))
    write_svg(path, lines)


def heat_color(value):
    value = max(-0.5, min(0.5, value))
    if value >= 0:
        start = (229, 231, 235)
        end = (15, 118, 110)
        frac = value / 0.5 if value else 0.0
    else:
        start = (229, 231, 235)
        end = (220, 38, 38)
        frac = abs(value) / 0.5
    r = int(start[0] + (end[0] - start[0]) * frac)
    g = int(start[1] + (end[1] - start[1]) * frac)
    b = int(start[2] + (end[2] - start[2]) * frac)
    return "#%02X%02X%02X" % (r, g, b)


def draw_per_gene_gain_heatmap(path, title, rows):
    genes = sorted({row["gene"] for row in rows}, key=gene_sort_key)
    modalities = sorted({row["modality"] for row in rows}, key=modality_sort_key)
    width = 220 + max(1, len(genes)) * 90
    height = 120 + max(1, len(modalities)) * 48
    margin_left, margin_top = 180, 80
    lines = svg_header(width, height)
    lines.append('<text x="%s" y="32" class="title">%s</text>' % (margin_left, title))
    for idx, gene in enumerate(genes):
        lines.append('<text x="%s" y="%s" text-anchor="middle" class="label">%s</text>' % (margin_left + idx * 90 + 40, margin_top - 12, gene))
    for row_idx, modality in enumerate(modalities):
        lines.append('<text x="20" y="%s" class="label">%s</text>' % (margin_top + row_idx * 48 + 24, modality_style(modality)["label"]))
        for col_idx, gene in enumerate(genes):
            record = next((row for row in rows if row["modality"] == modality and row["gene"] == gene), None)
            value = float(record["gain_vs_majority"]) if record and record["gain_vs_majority"] != "" else 0.0
            x = margin_left + col_idx * 90
            y = margin_top + row_idx * 48
            lines.append('<rect x="%s" y="%s" width="80" height="34" fill="%s" rx="4" />' % (x, y, heat_color(value)))
            lines.append('<text x="%s" y="%s" text-anchor="middle" class="small">%.2f</text>' % (x + 40, y + 21, value))
    write_svg(path, lines)


def draw_calibration_chart(path, title, rows):
    width, height = 1200, 560
    margin_left, margin_top, margin_bottom = 100, 80, 80
    chart_width = width - margin_left - 40
    chart_height = height - margin_top - margin_bottom
    lines = svg_header(width, height)
    lines.append('<text x="%s" y="32" class="title">%s</text>' % (margin_left, title))
    lines.append('<line x1="%s" y1="%s" x2="%s" y2="%s" class="axis" />' % (margin_left, margin_top + chart_height, margin_left + chart_width, margin_top + chart_height))
    lines.append('<line x1="%s" y1="%s" x2="%s" y2="%s" class="axis" />' % (margin_left, margin_top, margin_left, margin_top + chart_height))
    lines.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="#9ca3af" stroke-dasharray="4 4" />' % (margin_left, margin_top + chart_height, margin_left + chart_width, margin_top))
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["modality"])].append(row)
    offsets = ["#1F6FEB", "#F97316", "#0F766E", "#9333EA", "#DC2626"]
    for idx, key in enumerate(sorted(grouped.keys(), key=lambda item: (modality_sort_key(item[1]), item[0]))):
        points = []
        for row in sorted(grouped[key], key=lambda item: item["bin_index"]):
            x = margin_left + chart_width * float(row["mean_confidence"])
            y = margin_top + chart_height - chart_height * float(row["observed_accuracy"])
            points.append((x, y))
        if len(points) > 1:
            lines.append('<polyline fill="none" stroke="%s" stroke-width="2" points="%s" />' % (offsets[idx % len(offsets)], " ".join("%s,%s" % p for p in points)))
        for x, y in points:
            lines.append('<circle cx="%s" cy="%s" r="4" fill="%s" />' % (x, y, offsets[idx % len(offsets)]))
        lines.append('<text x="%s" y="%s" class="small">%s %s</text>' % (width - 240, 90 + idx * 16, key[0], key[1]))
    write_svg(path, lines)


def draw_abstention_tradeoff(path, title, rows):
    width, height = 1100, 560
    margin_left, margin_top, margin_bottom = 100, 80, 90
    chart_width = width - margin_left - 40
    chart_height = height - margin_top - margin_bottom
    lines = svg_header(width, height)
    lines.append('<text x="%s" y="32" class="title">%s</text>' % (margin_left, title))
    lines.append('<line x1="%s" y1="%s" x2="%s" y2="%s" class="axis" />' % (margin_left, margin_top + chart_height, margin_left + chart_width, margin_top + chart_height))
    lines.append('<line x1="%s" y1="%s" x2="%s" y2="%s" class="axis" />' % (margin_left, margin_top, margin_left, margin_top + chart_height))
    pts = []
    for row in rows:
        x = margin_left + chart_width * float(row["call_rate"])
        y = margin_top + chart_height - chart_height * float(row["accuracy_among_called"])
        pts.append((x, y, row["min_support"]))
    if len(pts) > 1:
        lines.append('<polyline fill="none" stroke="#2563eb" stroke-width="2" points="%s" />' % " ".join("%s,%s" % (x, y) for x, y, _ in pts))
    for x, y, threshold in pts:
        lines.append('<circle cx="%s" cy="%s" r="4" fill="#2563eb" />' % (x, y))
        lines.append('<text x="%s" y="%s" class="small">%.2f</text>' % (x + 6, y - 6, threshold))
    write_svg(path, lines)


def draw_discordance_taxonomy(path, title, rows):
    width, height = 1100, 560
    margin_left, margin_top, margin_bottom = 100, 80, 100
    chart_width = width - margin_left - 40
    chart_height = height - margin_top - margin_bottom
    lines = svg_header(width, height)
    lines.append('<text x="%s" y="32" class="title">%s</text>' % (margin_left, title))
    bars = rows or [{"scope": "none", "tag": "none", "n_events": 0}]
    max_value = max(int(row["n_events"]) for row in bars) or 1
    bar_width = chart_width / float(max(len(bars), 1)) * 0.7
    for idx, row in enumerate(bars):
        x = margin_left + idx * (chart_width / float(max(len(bars), 1))) + 10
        value = int(row["n_events"])
        bar_h = chart_height * (value / float(max_value))
        y = margin_top + chart_height - bar_h
        lines.append('<rect x="%s" y="%s" width="%s" height="%s" fill="#7c3aed" rx="4" />' % (x, y, bar_width, bar_h))
        lines.append('<text x="%s" y="%s" text-anchor="middle" class="small">%s</text>' % (x + bar_width / 2.0, y - 6 if value else y + 14, value))
        lines.append('<text x="%s" y="%s" text-anchor="middle" class="small">%s</text>' % (x + bar_width / 2.0, margin_top + chart_height + 20, row["tag"]))
    write_svg(path, lines)


def draw_grouped_rate_chart(path, title, records, value_key, ylabel):
    width, height = 1080, 520
    margin_left, margin_top, margin_bottom = 110, 80, 90
    chart_width = width - margin_left - 40
    chart_height = height - margin_top - margin_bottom
    tools = sorted({row["tool"] for row in records})
    modalities = ordered_modalities(records)
    lines = svg_header(width, height)
    lines.append('<text x="%s" y="32" class="title">%s</text>' % (margin_left, title))
    group_width = chart_width / float(max(len(tools), 1))
    bar_width = min(26, group_width / float(max(len(modalities), 1) + 1))
    for idx, tool in enumerate(tools):
        center_x = margin_left + group_width * idx + group_width / 2.0
        lines.append('<text x="%s" y="%s" text-anchor="middle" class="label">%s</text>' % (center_x, margin_top + chart_height + 28, tool))
        start_x = center_x - ((len(modalities) * bar_width + (len(modalities) - 1) * 6) / 2.0)
        for mod_idx, modality in enumerate(modalities):
            record = next((row for row in records if row["tool"] == tool and row["modality"] == modality), None)
            value = float(record[value_key]) if record and record[value_key] != "" else 0.0
            bar_height = chart_height * value
            x = start_x + mod_idx * (bar_width + 6)
            y = margin_top + chart_height - bar_height
            lines.append('<rect x="%s" y="%s" width="%s" height="%s" fill="%s" rx="4" />' % (x, y, bar_width, bar_height, modality_style(modality)["fill"]))
    lines.append('<text x="%s" y="%s" class="label">%s</text>' % (margin_left - 70, margin_top - 20, ylabel))
    write_svg(path, lines)


def write_caption_stub(path, metadata):
    text = """# Figure Captions

## Benchmark Figure Set

Truth source: %s
Primary allele resolution: two-field
IMGT/HLA version: %s
Weight version: %s

### Figure 2. Method comparison
Compares single-tool baselines, strict majority vote, and weighted consensus across modalities.

### Figure 3. Per-gene gains
Shows the gain of weighted consensus over majority vote for each gene and modality.

### Figure 4. Confidence calibration
Plots observed correctness against calibrated confidence bins for each tool and modality.

### Figure 5. Abstention tradeoff
Shows how increasing the minimum support threshold trades call rate for accuracy among emitted calls.

### Figure 6. Discordance taxonomy
Summarizes low-evidence, technical, and DNA/RNA discordance events detected in the benchmark cohort.

### Figure 7. Confidence-calibrated weights
Visualizes the learned per-tool/per-modality benchmark weights used by the weighted consensus model.
""" % (metadata["truth_source"], metadata.get("imgt_hla_version", "unspecified"), WEIGHT_VERSION)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def generate_figures(output_dir, method_comparison_rows, per_gene_gain_rows, calibration_bin_rows, abstention_rows, discordance_summary_rows, tool_weights):
    figures_dir = output_dir / "figures"
    draw_method_accuracy_chart(figures_dir / "figure_2_accuracy_comparison.svg", "Figure 2. Accuracy Comparison", method_comparison_rows)
    draw_per_gene_gain_heatmap(figures_dir / "figure_3_per_gene_gains.svg", "Figure 3. Per-Gene Gains", per_gene_gain_rows)
    draw_calibration_chart(figures_dir / "figure_4_confidence_calibration.svg", "Figure 4. Confidence Calibration", calibration_bin_rows)
    draw_abstention_tradeoff(figures_dir / "figure_5_abstention_tradeoff.svg", "Figure 5. Abstention Tradeoff", abstention_rows)
    draw_discordance_taxonomy(figures_dir / "figure_6_discordance_taxonomy.svg", "Figure 6. Discordance Taxonomy", discordance_summary_rows)
    draw_grouped_rate_chart(figures_dir / "figure_7_confidence_weights.svg", "Figure 7. Confidence-Calibrated Tool Weights", tool_weights, "final_weight", "Final weight")


def main():
    args = parse_args()
    config_path = Path(args.config)
    config = load_yaml(config_path)
    config = resolve_config_paths(config, config_path.parent)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    truth = load_truth(config["truth"])
    runs = apply_truth_source_guardrail(load_run_records(config), config)
    resolution = int(config.get("benchmark", {}).get("resolution", 2))
    imgt_hla_version = resolve_imgt_hla_version(config)
    mode = benchmark_mode(config, args.benchmark_mode)
    harmonized_rows = dedupe_rows(build_harmonized_rows(truth, runs, resolution=resolution, imgt_hla_version=imgt_hla_version))
    if not harmonized_rows:
        raise SystemExit("No harmonized benchmark rows were produced. Check the config globs and parser settings.")
    benchmark_genes = resolve_benchmark_genes(config)
    filtered_rows = filter_rows(harmonized_rows, genes=benchmark_genes or None)
    if not filtered_rows:
        raise SystemExit("No benchmark rows remained after applying the configured benchmark loci.")
    shared_genes = calculate_shared_genes(filtered_rows)
    main_rows = filtered_rows
    annotate_ciwd(main_rows)
    attach_probabilistic_calibration(main_rows, mode, config)
    # Optionally annotate rows with superpopulation for ancestry-stratified analysis
    pop_map = {}
    if args.population_manifest:
        pop_map = load_population_manifest(Path(args.population_manifest))
    for row in main_rows:
        row["superpopulation"] = pop_map.get(row["sample"], "unknown")
    summary = aggregate_metrics(main_rows)
    per_gene = build_per_gene_summary(main_rows)
    modality_gene = build_modality_gene_summary(main_rows)
    cohort = build_cohort_overview(main_rows, truth, benchmark_genes or shared_genes)
    ambiguity_summary_rows, ambiguity_summary_gene_rows = build_ambiguity_summary(main_rows)
    tool_weights, gene_weights = build_confidence_weights(main_rows, config=config, weight_alpha=args.weight_alpha, weight_beta=args.weight_beta, mode=mode)
    failure_mode_rows = build_failure_mode_summary(main_rows)
    cv_weight_rows = run_cross_validation(main_rows, args.cv_folds, config, args.weight_alpha, args.weight_beta, mode=mode) if args.cv_folds > 1 else []
    pop_weight_rows = build_population_weights(main_rows, config=config, weight_alpha=args.weight_alpha, weight_beta=args.weight_beta) if args.population_manifest else []
    pop_diagnostic_rows = build_population_diagnostics(main_rows, mode) if args.population_manifest else []
    pop_calibration_rows = build_population_calibration(main_rows, mode) if args.population_manifest else []
    runtime_weights = build_runtime_weight_payload(tool_weights, gene_weights, weight_alpha=args.weight_alpha, weight_beta=args.weight_beta, mode=mode, meta_params=meta_consensus_parameters(config))
    runtime_override = load_runtime_weight_override(config)
    if runtime_override:
        runtime_weights = runtime_override
    majority_rows = build_majority_vote_rows(main_rows)
    weighted_rows = build_weighted_consensus_rows(main_rows, runtime_weights, config)
    meta_rows, meta_trace_rows = build_meta_consensus_rows(main_rows, runtime_weights, config)
    champion_challenger_rows, champion_challenger_trace_rows, champion_challenger_method_rows, champion_challenger_method_per_gene_rows, champion_challenger_summary_rows = build_champion_challenger_outputs(
        main_rows,
        runtime_weights,
        config,
        mode,
        benchmark_genes or shared_genes,
    )
    weighted_threshold_sweep_rows = run_weighted_threshold_sweep(main_rows, runtime_weights, config)
    champion_override_sweep_rows = run_champion_override_sweep(main_rows, runtime_weights, config, benchmark_genes or shared_genes, mode)
    gated_rows, gated_trace_rows, gated_method_rows, gated_method_per_gene_rows, gated_summary_rows = build_gated_consensus_outputs(
        main_rows,
        majority_rows,
        weighted_rows,
        champion_challenger_rows,
        config,
    )
    locus_expert_rows, locus_expert_trace_rows, locus_expert_method_rows, locus_expert_method_per_gene_rows, locus_expert_summary_rows = build_locus_expert_outputs(
        main_rows,
        runtime_weights,
        config,
        mode,
        benchmark_genes or shared_genes,
    )
    bimodal_rows = build_bimodal_consensus_rows(main_rows, runtime_weights, config)
    bimodal_comparison_rows = build_bimodal_accuracy_comparison(main_rows, bimodal_rows, runtime_weights, config)
    method_ambiguity_rows, method_ambiguity_gene_rows = build_method_ambiguity_summary(majority_rows, weighted_rows, meta_rows)
    method_comparison_rows = build_method_comparison(summary, majority_rows, weighted_rows)
    meta_method_comparison_rows = summarize_consensus_method(meta_rows)
    method_per_gene_rows = build_method_per_gene(per_gene, majority_rows, weighted_rows, meta_rows)
    ablation_method_rows, ablation_method_per_gene_rows, ablation_summary_rows = build_ensemble_ablation_outputs(
        main_rows,
        config,
        mode,
        args.weight_alpha,
        args.weight_beta,
    )
    per_gene_gain_rows = build_per_gene_gain_table(method_per_gene_rows)
    multi_resolution_tool_rows = summarize_multiresolution_rows(main_rows, "tool", include_gene=False)
    multi_resolution_tool_gene_rows = summarize_multiresolution_rows(main_rows, "tool", include_gene=True)
    multi_resolution_method_rows = summarize_multiresolution_rows(majority_rows + weighted_rows + meta_rows, "method", include_gene=False)
    calibration_bin_rows = build_confidence_bin_summary(main_rows)
    calibration_summary_rows = build_confidence_calibration_summary(main_rows, mode=mode)
    cv_calibration_rows = build_cross_validation_weight_summary(main_rows, mode)
    confidence_error_rows, confidence_error_gene_rows = build_confidence_error_summary(main_rows)
    abstention_rows = build_abstention_tradeoff(main_rows, runtime_weights, config)
    discordance_rows = build_discordance_rows(main_rows, weighted_rows)
    discordance_summary_rows = summarize_discordance(discordance_rows)
    arbitration_accuracy_rows = build_arbitration_accuracy_table(discordance_rows, main_rows)
    arbitration_rule_summary_rows = summarize_arbitration_rules(arbitration_accuracy_rows)
    pairwise_agreement_rows, pairwise_agreement_gene_rows = build_pairwise_agreement(main_rows)
    error_taxonomy_rows = build_tool_vs_truth_error_taxonomy(main_rows)
    disagreement_event_rows = build_tool_disagreement_events(main_rows, majority_rows, weighted_rows)
    consensus_trace_rows = build_consensus_decision_trace(main_rows, majority_rows, weighted_rows, runtime_weights, mode)
    single_tool_method_rows = []
    for row in main_rows:
        single_tool_method_rows.append(dict(row, method=row["tool"]))
    population_method_rows = build_population_method_comparison(single_tool_method_rows + majority_rows + weighted_rows + meta_rows, mode)
    locus_difficulty_rows = build_locus_difficulty_summary(main_rows, method_per_gene_rows)
    population_resolution_rows = build_population_resolution_summary(main_rows, mode)
    population_conflict_rows = build_population_conflict_summary(main_rows, mode)
    agreement_priors_by_gene = defaultdict(dict)
    gene_group = defaultdict(list)
    for row in pairwise_agreement_gene_rows:
        gene_group[(row["modality"], row["gene"])].append(float(row["exact_2field_agreement_rate"]))
    for (modality, gene), values in gene_group.items():
        agreement_priors_by_gene[modality][gene] = round(statistics.mean(values), 4)
    runtime_weights["agreement_priors_by_gene"] = {modality: dict(values) for modality, values in agreement_priors_by_gene.items()}
    metadata = build_metadata(config, shared_genes, harmonized_rows, tool_weights, weight_alpha=args.weight_alpha, weight_beta=args.weight_beta, mode=mode)
    sweep_modality = sweep_target_modality(main_rows)
    sweep_prefix = sweep_file_prefix(sweep_modality)
    observed_sweep_tools = sorted({row["tool"] for row in main_rows if row["modality"] == sweep_modality})
    expected_sweep_tools = sorted({run.tool for run in runs if run.modality == sweep_modality})
    missing_sweep_tools = [tool for tool in expected_sweep_tools if tool not in observed_sweep_tools]
    weighted_best_by_accuracy = max(weighted_threshold_sweep_rows, key=lambda row: (float(row["overall_correct_call_rate"]), float(row["callable_rate"])), default={})
    weighted_best_by_callable = max(weighted_threshold_sweep_rows, key=lambda row: (float(row["callable_rate"]), float(row["overall_correct_call_rate"])), default={})
    champion_best_by_accuracy = max(champion_override_sweep_rows, key=lambda row: (float(row["overall_correct_call_rate"]), float(row["callable_rate"])), default={})
    champion_best_by_override = max(champion_override_sweep_rows, key=lambda row: (int(row["corrective_override_count"]) - int(row["harmful_override_count"]), -int(row["harmful_override_count"]), float(row["overall_correct_call_rate"])), default={})
    weighted_baseline = next((row for row in method_comparison_rows if row.get("method") == "WeightedConsensus" and row.get("modality") == sweep_modality), {})
    champion_baseline = next((row for row in champion_challenger_method_rows if row.get("method") == "ChampionChallenger" and row.get("modality") == sweep_modality), {})
    weighted_baseline_accuracy = coerce_float(weighted_baseline.get("overall_correct_call_rate"))
    champion_baseline_accuracy = coerce_float(champion_baseline.get("overall_correct_call_rate"))
    weighted_best_accuracy = coerce_float(weighted_best_by_accuracy.get("overall_correct_call_rate"))
    champion_best_accuracy = coerce_float(champion_best_by_accuracy.get("overall_correct_call_rate"))
    modality_label = modality_style(sweep_modality)["label"] if sweep_modality else "benchmark"
    forced_spechla_accuracy = spechla_forced_accuracy(config, sweep_modality)
    interpretation_parts = []
    if weighted_baseline_accuracy is not None and weighted_best_accuracy is not None:
        if weighted_best_accuracy > weighted_baseline_accuracy:
            interpretation_parts.append("WeightedConsensus is threshold-limited on %s: support/margin retuning improves overall correct-call rate from %.4f to %.4f." % (modality_label, weighted_baseline_accuracy, weighted_best_accuracy))
        else:
            interpretation_parts.append("WeightedConsensus appears near its current %s threshold ceiling: the sweep does not improve overall correct-call rate beyond %.4f." % (modality_label, weighted_baseline_accuracy))
    if champion_baseline_accuracy is not None and champion_best_accuracy is not None:
        if champion_best_accuracy > champion_baseline_accuracy:
            interpretation_parts.append("ChampionChallenger has a useful %s override regime: the best override setting improves overall correct-call rate from %.4f to %.4f." % (modality_label, champion_baseline_accuracy, champion_best_accuracy))
        else:
            interpretation_parts.append("ChampionChallenger does not show a better %s override regime than its baseline setting." % modality_label)
    if missing_sweep_tools:
        interpretation_parts.append("%s remains coverage-limited for tools missing benchmark rows: %s." % (modality_label, ",".join(missing_sweep_tools)))
    elif (
        weighted_baseline_accuracy is not None and weighted_best_accuracy is not None and weighted_best_accuracy > weighted_baseline_accuracy
    ) or (
        champion_baseline_accuracy is not None and champion_best_accuracy is not None and champion_best_accuracy > champion_baseline_accuracy
    ):
        interpretation_parts.append("Residual %s limitations are better explained by decision thresholds than by missing tool coverage on this cohort." % modality_label)
    else:
        interpretation_parts.append("Residual %s limitations look more structural than coverage-driven on this cohort." % modality_label)
    if forced_spechla_accuracy:
        interpretation_parts.append("SpecHLA was forcibly included in headline RNA accuracy scoring for this run.")
    sweep_summary = {
        "modality": sweep_modality,
        "benchmark_config": str(config_path),
        "%s_run_root" % sweep_prefix: str(output_dir),
        "forced_spechla_accuracy": forced_spechla_accuracy,
        "weighted_sweep_best_by_accuracy": weighted_best_by_accuracy,
        "weighted_sweep_best_by_callable_rate": weighted_best_by_callable,
        "champion_sweep_best_by_accuracy": champion_best_by_accuracy,
        "champion_sweep_best_by_override_quality": champion_best_by_override,
        "coverage_notes": {
            "expected_%s_tools" % sweep_prefix: expected_sweep_tools,
            "observed_%s_tools" % sweep_prefix: observed_sweep_tools,
            "missing_%s_tools" % sweep_prefix: missing_sweep_tools,
            "spechla_forced_into_accuracy": forced_spechla_accuracy,
        },
        "interpretation": " ".join(interpretation_parts),
    }

    write_tsv(output_dir / "tables" / "harmonized_benchmark_rows.tsv", main_rows, ["sample", "superpopulation", "modality", "tool", "gene", "truth_allele1_raw", "truth_allele2_raw", "allele1_raw", "allele2_raw", "truth_allele1", "truth_allele2", "allele1", "allele2", "truth_allele1_3field", "truth_allele2_3field", "allele1_3field", "allele2_3field", "call_status", "correct_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_correct_g_group", "is_correct_p_group", "is_ambiguity_compatible", "is_resolution_compatible", "compatibility_grade", "match_grade", "truth_ciwd_1", "truth_ciwd_2", "truth_ciwd_stratum", "call_implausible", "imgt_hla_version", "runtime_hours", "max_ram_gb", "confidence_score", "confidence_source", "raw_confidence", "read_support", "raw_score_family", "raw_score_value", "calibrated_probability", "cv_calibrated_probability", "calibration_method", "source_file"])
    write_tsv(output_dir / "tables" / "summary_full_cohort.tsv", sorted(summary.values(), key=lambda row: (modality_sort_key(row["modality"]), row["tool"])), ["tool", "modality", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi", "median_runtime_hours", "median_max_ram_gb"])
    write_tsv(output_dir / "tables" / "summary_per_gene.tsv", per_gene, ["tool", "modality", "gene", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"])
    write_tsv(output_dir / "tables" / "summary_full_cohort_multiresolution.tsv", multi_resolution_tool_rows, ["tool", "modality", "sample_count", "gene_rows", "callable_rate", "exact_2field_rate", "exact_3field_rate", "g_group_match_rate", "p_group_match_rate", "ambiguity_compatible_rate"])
    write_tsv(output_dir / "tables" / "summary_per_gene_multiresolution.tsv", multi_resolution_tool_gene_rows, ["tool", "modality", "gene", "sample_count", "gene_rows", "callable_rate", "exact_2field_rate", "exact_3field_rate", "g_group_match_rate", "p_group_match_rate", "ambiguity_compatible_rate"])
    write_tsv(output_dir / "tables" / "summary_modality_gene_coverage.tsv", modality_gene, ["modality", "gene"])
    write_tsv(output_dir / "tables" / "cohort_overview.tsv", cohort, ["modality", "tool", "truth_samples", "completed_tool_samples", "modality_samples", "shared_gene_count"])
    write_tsv(output_dir / "tables" / "ambiguity_summary.tsv", ambiguity_summary_rows, ["tool", "modality", "imgt_hla_version", "sample_count", "gene_rows", "exact_2field_rate", "exact_3field_rate", "g_group_comparable_rows", "g_group_match_rate", "p_group_comparable_rows", "p_group_match_rate"])
    write_tsv(output_dir / "tables" / "ambiguity_summary_by_gene.tsv", ambiguity_summary_gene_rows, ["tool", "modality", "gene", "imgt_hla_version", "sample_count", "gene_rows", "exact_2field_rate", "exact_3field_rate", "g_group_comparable_rows", "g_group_match_rate", "p_group_comparable_rows", "p_group_match_rate"])
    write_tsv(output_dir / "tables" / "reference_metadata.tsv", [{"truth_source": metadata.get("truth_source", ""), "truth_path": metadata.get("truth_path", ""), "imgt_hla_version": metadata.get("imgt_hla_version", ""), "primary_resolution": metadata.get("resolution", ""), "secondary_resolutions": ",".join(str(x) for x in metadata.get("ambiguity_evaluation", {}).get("secondary_resolutions", [])), "benchmark_mode": mode}], ["truth_source", "truth_path", "imgt_hla_version", "primary_resolution", "secondary_resolutions", "benchmark_mode"])
    write_tsv(output_dir / "tables" / "sample_level_disagreements.tsv", build_sample_disagreements(harmonized_rows), ["sample", "modality", "tool", "gene", "truth_allele1", "truth_allele2", "allele1", "allele2", "call_status", "correct_status", "source_file"])
    write_tsv(output_dir / "tables" / "missing_call_patterns.tsv", build_missing_patterns(harmonized_rows), ["tool", "modality", "gene", "missing_calls", "total_rows", "missing_rate"])
    write_tsv(output_dir / "tables" / "truth_mismatches.tsv", build_truth_mismatches(harmonized_rows), ["sample", "modality", "tool", "gene", "truth_allele1", "truth_allele2", "typed_allele1", "typed_allele2", "call_status"])
    write_tsv(output_dir / "tables" / "tool_pairwise_agreement.tsv", pairwise_agreement_rows, ["modality", "tool_a", "tool_b", "n_shared_callable_rows", "exact_2field_agreement_rate", "exact_3field_agreement_rate", "g_group_agreement_rate", "p_group_agreement_rate"])
    write_tsv(output_dir / "tables" / "tool_pairwise_agreement_by_gene.tsv", pairwise_agreement_gene_rows, ["modality", "gene", "tool_a", "tool_b", "n_shared_callable_rows", "exact_2field_agreement_rate", "exact_3field_agreement_rate", "g_group_agreement_rate", "p_group_agreement_rate"])
    write_tsv(output_dir / "tables" / "tool_vs_truth_error_taxonomy.tsv", error_taxonomy_rows, ["tool", "modality", "gene", "error_category", "count", "fraction"])
    write_tsv(output_dir / "tables" / "tool_disagreement_events.tsv", disagreement_event_rows, ["sample", "population", "modality", "gene", "truth_allele1", "truth_allele2", "callable_tool_count", "distinct_pair_count", "tool_calls_json", "majority_pair", "weighted_pair", "majority_correct", "weighted_correct", "disagreement_class"])
    write_tsv(output_dir / "tables" / "consensus_decision_trace.tsv", consensus_trace_rows, ["sample", "population", "modality", "gene", "benchmark_mode", "method", "winning_pair", "runner_up_pair", "winning_support_fraction", "winning_support_margin", "contributing_tools", "top_tool_weights_json", "call_status", "is_correct"])
    write_tsv(output_dir / "tables" / "tool_confidence_weights.tsv", tool_weights, ["tool", "modality", "benchmark_mode", "sample_count", "gene_rows", "callable_rate", "confidence_coverage_rate", "base_reliability", "reliability_ci_lo", "reliability_ci_hi", "reliability_ci_width", "raw_score_family", "calibration_method", "calibration_sample_size", "calibrated_confidence", "mean_calibrated_probability", "cv_mean_calibrated_probability", "effective_confidence", "guardrail_factor", "guardrail_status", "brier_score", "expected_calibration_error", "calibrated_brier_score", "calibrated_expected_calibration_error", "final_weight", "weight_version"])
    write_tsv(output_dir / "tables" / "tool_confidence_weights_by_gene.tsv", gene_weights, ["tool", "modality", "gene", "benchmark_mode", "sample_count", "gene_rows", "callable_rate", "confidence_coverage_rate", "base_reliability", "reliability_ci_lo", "reliability_ci_hi", "reliability_ci_width", "raw_score_family", "calibration_method", "calibration_sample_size", "calibrated_confidence", "mean_calibrated_probability", "cv_mean_calibrated_probability", "effective_confidence", "guardrail_factor", "guardrail_status", "brier_score", "expected_calibration_error", "calibrated_brier_score", "calibrated_expected_calibration_error", "final_weight", "weight_version"])
    write_tsv(output_dir / "tables" / "failure_mode_summary.tsv", failure_mode_rows, ["tool", "modality", "gene", "failure_mode", "count", "fraction"])
    if cv_weight_rows:
        write_tsv(output_dir / "tables" / "tool_confidence_weights_cv.tsv", cv_weight_rows, ["tool", "modality", "n_folds", "weight_mean", "weight_std", "weight_cv", "weight_min", "weight_max"])
    write_tsv(output_dir / "tables" / "cross_validation_weight_summary.tsv", cv_calibration_rows, ["tool", "modality", "benchmark_mode", "n_rows", "cv_n_rows", "calibration_method", "mean_calibrated_probability", "cv_mean_calibrated_probability", "calibrated_brier_score", "calibrated_expected_calibration_error"])
    if pop_weight_rows:
        write_tsv(output_dir / "tables" / "tool_confidence_weights_by_population.tsv", pop_weight_rows, ["superpopulation", "tool", "modality", "benchmark_mode", "sample_count", "gene_rows", "callable_rate", "confidence_coverage_rate", "base_reliability", "reliability_ci_lo", "reliability_ci_hi", "reliability_ci_width", "raw_score_family", "calibration_method", "calibration_sample_size", "calibrated_confidence", "mean_calibrated_probability", "cv_mean_calibrated_probability", "effective_confidence", "guardrail_factor", "guardrail_status", "brier_score", "expected_calibration_error", "calibrated_brier_score", "calibrated_expected_calibration_error", "final_weight", "weight_version"])
    if pop_diagnostic_rows:
        write_tsv(output_dir / "tables" / "tool_population_diagnostics.tsv", pop_diagnostic_rows, ["superpopulation", "tool", "modality", "benchmark_mode", "sample_count", "gene_rows", "callable_rate", "overall_correct_call_rate", "confidence_coverage_rate"])
    if pop_calibration_rows:
        write_tsv(output_dir / "tables" / "tool_population_calibration.tsv", pop_calibration_rows, ["superpopulation", "tool", "modality", "benchmark_mode", "n_rows", "mean_confidence", "brier_score", "expected_calibration_error", "mean_calibrated_probability", "calibrated_brier_score", "calibrated_expected_calibration_error"])
    write_tsv(output_dir / "tables" / "population_resolution_summary.tsv", population_resolution_rows, ["superpopulation", "tool", "modality", "benchmark_mode", "gene_rows", "exact_2field_rate", "exact_3field_rate", "ambiguity_compatible_rate", "resolution_compatible_rate"])
    write_tsv(output_dir / "tables" / "population_conflict_summary.tsv", population_conflict_rows, ["superpopulation", "tool", "modality", "benchmark_mode", "conflict_category", "count", "fraction"])
    write_tsv(output_dir / "tables" / "majority_vote_baseline.tsv", majority_rows, ["sample", "modality", "gene", "method", "truth_allele1", "truth_allele2", "allele1", "allele2", "truth_allele1_3field", "truth_allele2_3field", "allele1_3field", "allele2_3field", "call_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_correct_g_group", "is_correct_p_group", "is_ambiguity_compatible", "is_resolution_compatible", "compatibility_grade", "match_grade", "imgt_hla_version", "agreeing_tools", "contributing_tools", "support_fraction", "support_margin", "discordance_tag"])
    write_tsv(output_dir / "tables" / "weighted_consensus_calls.tsv", weighted_rows, ["sample", "modality", "gene", "method", "truth_allele1", "truth_allele2", "allele1", "allele2", "truth_allele1_3field", "truth_allele2_3field", "allele1_3field", "allele2_3field", "call_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_correct_g_group", "is_correct_p_group", "is_ambiguity_compatible", "is_resolution_compatible", "compatibility_grade", "match_grade", "imgt_hla_version", "agreeing_tools", "contributing_tools", "support_fraction", "support_margin", "total_weight", "discordance_tag"])
    write_tsv(output_dir / "tables" / "meta_consensus_calls.tsv", meta_rows, ["sample", "modality", "gene", "method", "truth_allele1", "truth_allele2", "allele1", "allele2", "truth_allele1_3field", "truth_allele2_3field", "allele1_3field", "allele2_3field", "call_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_correct_g_group", "is_correct_p_group", "is_ambiguity_compatible", "is_resolution_compatible", "compatibility_grade", "match_grade", "imgt_hla_version", "agreeing_tools", "contributing_tools", "support_fraction", "support_margin", "total_weight", "discordance_tag"])
    write_tsv(output_dir / "tables" / "meta_consensus_decision_trace.tsv", meta_trace_rows, ["sample", "population", "modality", "gene", "benchmark_mode", "method", "winning_pair", "runner_up_pair", "winning_support_fraction", "winning_support_margin", "contributing_tools", "top_tool_weights_json", "call_status", "is_correct", "support_weight_sum", "mean_calibrated_probability", "agreement_bonus", "gene_strength_bonus", "ambiguity_bonus", "fragmentation_penalty", "meta_score"])
    write_tsv(output_dir / "tables" / "champion_challenger_calls.tsv", champion_challenger_rows, ["sample", "population", "modality", "gene", "method", "truth_allele1", "truth_allele2", "allele1", "allele2", "truth_allele1_3field", "truth_allele2_3field", "allele1_3field", "allele2_3field", "call_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_correct_g_group", "is_correct_p_group", "is_ambiguity_compatible", "is_resolution_compatible", "compatibility_grade", "match_grade", "imgt_hla_version", "agreeing_tools", "contributing_tools", "support_fraction", "support_margin", "total_weight", "discordance_tag", "champion_tool", "override_triggered", "override_reason"])
    write_tsv(output_dir / "tables" / "champion_challenger_trace.tsv", champion_challenger_trace_rows, ["sample", "population", "modality", "gene", "method", "champion_tool", "champion_pair", "challenger_pair", "override_triggered", "override_reason", "call_status", "is_correct", "challenger_support_fraction", "challenger_support_margin", "supporting_tools", "top_tool_weights_json"])
    write_tsv(output_dir / "tables" / "gated_consensus_calls.tsv", gated_rows, ["sample", "population", "modality", "gene", "method", "truth_allele1", "truth_allele2", "allele1", "allele2", "call_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_ambiguity_compatible", "is_resolution_compatible", "compatibility_grade", "agreeing_tools", "contributing_tools", "support_fraction", "support_margin", "total_weight", "difficulty_class", "selected_policy"])
    write_tsv(output_dir / "tables" / "gated_consensus_trace.tsv", gated_trace_rows, ["sample", "population", "modality", "gene", "method", "difficulty_class", "selected_policy", "distinct_pair_count", "callable_tool_count", "weighted_winning_pair", "weighted_runner_up_pair", "weighted_support_fraction", "weighted_support_margin", "majority_pair", "majority_weighted_disagree", "weighted_top_pair_is_duplicated", "alternative_nonduplicated_exists", "hardness_reasons", "call_status", "is_correct"])
    write_tsv(output_dir / "tables" / "locus_expert_consensus_calls.tsv", locus_expert_rows, ["panel_name", "parent_method", "sample", "population", "modality", "gene", "method", "truth_allele1", "truth_allele2", "allele1", "allele2", "call_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_ambiguity_compatible", "is_resolution_compatible", "compatibility_grade", "agreeing_tools", "contributing_tools", "support_fraction", "support_margin", "total_weight", "selected_tools"])
    write_tsv(output_dir / "tables" / "locus_expert_decision_trace.tsv", locus_expert_trace_rows, ["panel_name", "parent_method", "sample", "population", "modality", "gene", "selected_tools", "winning_pair", "runner_up_pair", "winning_support_fraction", "winning_support_margin", "contributing_tools", "top_tool_weights_json", "call_status", "is_correct", "support_weight_sum", "mean_calibrated_probability", "agreement_bonus", "gene_strength_bonus", "ambiguity_bonus", "fragmentation_penalty", "meta_score"])
    write_tsv(output_dir / "tables" / "bimodal_wes_rna_consensus.tsv", bimodal_rows, ["sample", "modality", "gene", "method", "truth_allele1", "truth_allele2", "allele1", "allele2", "call_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_ambiguity_compatible", "is_resolution_compatible", "compatibility_grade", "agreeing_tools", "contributing_tools", "support_fraction", "support_margin", "total_weight", "n_wes_tools", "n_rna_tools", "discordance_tag"])
    write_tsv(output_dir / "tables" / "bimodal_accuracy_comparison.tsv", bimodal_comparison_rows, ["comparison", "modality", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"])
    write_tsv(output_dir / "tables" / "method_comparison.tsv", method_comparison_rows, ["method", "method_type", "modality", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"])
    # CIWD 3.0.0 auxiliary layer: concordance stratified by allele commonness, and a plausibility QC
    # flag. CIWD is a classification catalogue, not a truth cohort — these are additive views over the
    # existing ground-truth benchmark (see assets/README_CIWD.md, docs/NEW_HLA_BENCHMARK_DATASETS.md).
    ciwd_stratified_rows = build_ciwd_stratified_summary(
        [(tool, [r for r in main_rows if r["tool"] == tool]) for tool in sorted({r["tool"] for r in main_rows})]
        + [("MajorityVote", majority_rows), ("WeightedConsensus", weighted_rows),
           ("ChampionChallenger", champion_challenger_rows)]
    )
    write_tsv(output_dir / "tables" / "summary_ciwd_stratified.tsv", ciwd_stratified_rows, ["method", "modality", "ciwd_stratum", "gene_rows", "callable_rate", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "summary_ciwd_plausibility.tsv", build_ciwd_plausibility_summary(main_rows), ["tool", "modality", "callable_calls", "implausible_calls", "implausible_rate", "implausible_incorrect", "implausible_error_rate"])
    write_tsv(output_dir / "tables" / "meta_method_comparison.tsv", [dict(row, method_type="ensemble") for row in meta_method_comparison_rows], ["method", "method_type", "modality", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"])
    write_tsv(output_dir / "tables" / "champion_challenger_method_comparison.tsv", champion_challenger_method_rows, ["method", "method_type", "modality", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"])
    write_tsv(output_dir / "tables" / "gated_method_comparison.tsv", gated_method_rows, ["method", "method_type", "modality", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"])
    write_tsv(output_dir / "tables" / "locus_expert_method_comparison.tsv", locus_expert_method_rows, ["panel_name", "parent_method", "method", "method_type", "modality", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"])
    write_tsv(output_dir / "tables" / "method_per_gene.tsv", method_per_gene_rows, ["method", "method_type", "modality", "gene", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "champion_challenger_method_per_gene.tsv", champion_challenger_method_per_gene_rows, ["method", "method_type", "modality", "gene", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "gated_method_per_gene.tsv", gated_method_per_gene_rows, ["method", "method_type", "modality", "gene", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "locus_expert_method_per_gene.tsv", locus_expert_method_per_gene_rows, ["panel_name", "parent_method", "method", "method_type", "modality", "gene", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "ablation_method_comparison.tsv", ablation_method_rows, ["ablation_name", "tool_subset", "method", "method_type", "modality", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"])
    write_tsv(output_dir / "tables" / "ablation_method_per_gene.tsv", ablation_method_per_gene_rows, ["ablation_name", "tool_subset", "method", "method_type", "modality", "gene", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "method_comparison_multiresolution.tsv", multi_resolution_method_rows, ["method", "modality", "sample_count", "gene_rows", "callable_rate", "exact_2field_rate", "exact_3field_rate", "g_group_match_rate", "p_group_match_rate", "ambiguity_compatible_rate"])
    write_tsv(output_dir / "tables" / "method_ambiguity_summary.tsv", method_ambiguity_rows, ["method", "modality", "imgt_hla_version", "sample_count", "gene_rows", "exact_2field_rate", "exact_3field_rate", "g_group_comparable_rows", "g_group_match_rate", "p_group_comparable_rows", "p_group_match_rate"])
    write_tsv(output_dir / "tables" / "method_ambiguity_summary_by_gene.tsv", method_ambiguity_gene_rows, ["method", "modality", "gene", "imgt_hla_version", "sample_count", "gene_rows", "exact_2field_rate", "exact_3field_rate", "g_group_comparable_rows", "g_group_match_rate", "p_group_comparable_rows", "p_group_match_rate"])
    write_tsv(output_dir / "tables" / "per_gene_gain.tsv", per_gene_gain_rows, ["modality", "gene", "weighted_correct_call_rate", "majority_correct_call_rate", "best_single_tool", "best_single_tool_rate", "gain_vs_majority", "gain_vs_best_single"])
    write_tsv(output_dir / "tables" / "population_method_comparison.tsv", population_method_rows, ["superpopulation", "method", "modality", "benchmark_mode", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "locus_difficulty_summary.tsv", locus_difficulty_rows, ["modality", "gene", "n_tools", "best_single_tool_rate", "tool_accuracy_spread", "disagreement_frequency", "weighted_vs_best_single_gain"])
    write_tsv(output_dir / "tables" / "confidence_bin_summary.tsv", calibration_bin_rows, ["tool", "modality", "bin_index", "bin_lower", "bin_upper", "n_rows", "mean_confidence", "observed_accuracy"])
    write_tsv(output_dir / "tables" / "confidence_calibration_summary.tsv", calibration_summary_rows, ["tool", "modality", "benchmark_mode", "calibration_method", "n_rows", "mean_confidence", "observed_accuracy", "brier_score", "expected_calibration_error", "mean_calibrated_probability", "calibrated_brier_score", "calibrated_expected_calibration_error"])
    write_tsv(output_dir / "tables" / "confidence_error_summary.tsv", confidence_error_rows, ["tool", "modality", "bin_index", "bin_lower", "bin_upper", "n_rows", "error_count", "error_rate", "observed_accuracy", "mean_confidence"])
    write_tsv(output_dir / "tables" / "confidence_error_summary_by_gene.tsv", confidence_error_gene_rows, ["tool", "modality", "gene", "bin_index", "bin_lower", "bin_upper", "n_rows", "error_count", "error_rate", "observed_accuracy", "mean_confidence"])
    write_tsv(output_dir / "tables" / "abstention_tradeoff.tsv", abstention_rows, ["min_support", "min_margin", "call_rate", "no_call_rate", "low_confidence_rate", "accuracy_among_called", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "discordance_tags.tsv", discordance_rows, [
        "sample", "gene", "scope", "tag", "detail",
        "arbitrated_pair", "arbitration_rule", "arbitration_outcome",
        "dna_callable_tools", "dna_support_fraction",
        "rna_callable_tools", "rna_support_fraction", "rna_mean_read_support",
    ])
    write_tsv(output_dir / "tables" / "discordance_summary.tsv", discordance_summary_rows, ["scope", "tag", "n_events"])
    write_tsv(output_dir / "tables" / "sample_discordance_flags.tsv",
              build_sample_discordance_flags(discordance_rows),
              ["sample", "n_discordant_loci", "discordant_loci",
               "n_expression_bias_loci", "expression_bias_loci",
               "integrity_flag", "flag_reason"])
    write_tsv(output_dir / "tables" / "arbitration_accuracy.tsv", arbitration_accuracy_rows, [
        "sample", "gene", "arbitration_rule", "arbitration_outcome",
        "arbitrated_pair", "truth_allele1", "truth_allele2", "is_correct",
    ])
    write_tsv(output_dir / "tables" / "arbitration_rule_summary.tsv", arbitration_rule_summary_rows, [
        "arbitration_rule", "n_resolved", "n_correct", "n_abstain", "accuracy_among_resolved",
    ])
    save_json(output_dir / "tables" / "benchmark_metadata.json", metadata)
    save_json(output_dir / "tables" / "consensus_runtime_weights.json", runtime_weights)
    save_json(output_dir / "tables" / "ablation_summary.json", ablation_summary_rows)
    save_json(output_dir / "tables" / "champion_challenger_summary.json", champion_challenger_summary_rows)
    save_json(output_dir / "tables" / "gated_summary.json", gated_summary_rows)
    save_json(output_dir / "tables" / "locus_expert_summary.json", locus_expert_summary_rows)
    if weighted_threshold_sweep_rows or champion_override_sweep_rows:
        write_tsv(output_dir / "sweeps" / ("%s_weighted_threshold_sweep.tsv" % sweep_prefix), weighted_threshold_sweep_rows, ["method", "modality", "min_support", "min_margin", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"])
        write_tsv(output_dir / "sweeps" / ("%s_champion_override_sweep.tsv" % sweep_prefix), champion_override_sweep_rows, ["method", "modality", "champion_A", "champion_B", "champion_C", "min_challenger_support_fraction", "min_challenger_margin", "min_supporting_tools", "require_non_ambiguity_override", "override_count", "corrective_override_count", "harmful_override_count", "neutral_override_count", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"])
        save_json(output_dir / "sweeps" / ("%s_sweep_summary.json" % sweep_prefix), sweep_summary)
    generate_figures(output_dir, method_comparison_rows, per_gene_gain_rows, calibration_bin_rows, abstention_rows, discordance_summary_rows, tool_weights)
    write_caption_stub(output_dir / "figures" / "captions.md", metadata)
    return 0


if __name__ == "__main__":
    sys.exit(main())
