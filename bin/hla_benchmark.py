#!/usr/bin/env python3
"""Benchmark HLA typing runs, cohort-aware confidence weights, and figure-ready outputs."""

import argparse
import csv
import glob
import json
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
DEFAULT_CONFIDENCE_GUARDRAIL = {
    "enabled": True,
    "max_expected_calibration_error_for_boost": 0.35,
    "max_brier_score_for_boost": 0.35,
    "min_confidence_coverage_for_boost": 0.5,
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
    else:
        match_grade = "mismatch"

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


def load_runtime_weight_override(config):
    benchmark_cfg = config.get("benchmark", {})
    path_value = clean_token(benchmark_cfg.get("runtime_weight_override", ""))
    if not path_value:
        return None
    with Path(path_value).open("r", encoding="utf-8") as handle:
        return json.load(handle)


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
    if parser_name in {"long_hla_table", "polysolver_table", "kourami_table", "t1k_table", "seq2hla_table"}:
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
        return clip_unit(raw_score_val), read_support_val, raw_score_val
    if read_support_val is not None:
        return clip_unit(read_support_val / float(target_reads)), read_support_val, None
    return None, None, raw_score_val


def build_confidence_entry(sample, gene, raw_score, read_support, defaults, source, source_file):
    confidence_score, read_support_val, raw_confidence_val = normalize_confidence_values(raw_score, read_support, defaults)
    return {
        "sample": sample,
        "gene": gene,
        "confidence_score": confidence_score,
        "read_support": read_support_val,
        "raw_confidence": raw_confidence_val,
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
                    "match_grade": ambiguity["match_grade"],
                    "imgt_hla_version": imgt_hla_version,
                    "runtime_hours": as_string_number(runtime_hours),
                    "max_ram_gb": as_string_number(max_ram_gb),
                    "confidence_score": as_string_number(confidence.get("confidence_score")),
                    "confidence_source": confidence.get("confidence_source", "missing"),
                    "raw_confidence": as_string_number(confidence.get("raw_confidence")),
                    "read_support": as_string_number(confidence.get("read_support")),
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


def filter_rows(rows, genes=None, samples=None):
    out = []
    for row in rows:
        if genes and row["gene"] not in genes:
            continue
        if samples is not None and row["sample"] not in samples:
            continue
        out.append(row)
    return out


def ratio(numerator, denominator):
    return round((float(numerator) / denominator), 4) if denominator else 0.0


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
    return {
        "sample_count": len({row["sample"] for row in entries}),
        "gene_rows": total,
        "callable_rate": ratio(callable_n, total),
        "accuracy_among_callable": ratio(correct_n, callable_n),
        "overall_correct_call_rate": ratio(correct_n, total),
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


def build_method_ambiguity_summary(majority_rows, weighted_rows):
    grouped = defaultdict(list)
    gene_grouped = defaultdict(list)
    for row in majority_rows + weighted_rows:
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


def calibration_stats(entries, bin_count=5):
    scored = [row for row in entries if coerce_float(row.get("confidence_score", "")) is not None]
    n = len(scored)
    if not n:
        return {"n_rows": 0, "mean_confidence": None, "observed_accuracy": None, "brier_score": None, "expected_calibration_error": None}
    grouped = defaultdict(list)
    for row in scored:
        grouped[calibration_bin_label(coerce_float(row.get("confidence_score", "")), bin_count)[0]].append(row)
    brier = sum((coerce_float(row["confidence_score"]) - int(row["is_correct"])) ** 2 for row in scored) / float(n)
    ece = 0.0
    for entries_in_bin in grouped.values():
        mean_conf = statistics.mean(coerce_float(row["confidence_score"]) for row in entries_in_bin)
        observed = statistics.mean(int(row["is_correct"]) for row in entries_in_bin)
        ece += abs(mean_conf - observed) * len(entries_in_bin)
    ece = ece / float(n)
    return {
        "n_rows": n,
        "mean_confidence": round(statistics.mean(coerce_float(row["confidence_score"]) for row in scored), 4),
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
    if brier > settings["max_brier_score_for_boost"] or ece > settings["max_expected_calibration_error_for_boost"]:
        return reliability, 0.0, "poor_calibration"
    shrink_factor = round(max(0.0, 1.0 - max(brier, ece)), 4)
    effective_confidence = round(reliability + shrink_factor * (calibrated_confidence - reliability), 4)
    return effective_confidence, shrink_factor, "applied"


def compute_weight_record(entries, tool, modality, gene=None, settings=None):
    summary = summarize_entries(entries)
    reliability = summary["overall_correct_call_rate"]
    coverage_rate = confidence_coverage_rate(entries)
    calibrated_confidence = mean_confidence(entries)
    calibration = calibration_stats(entries)
    effective_confidence, guardrail_factor, guardrail_status = guarded_confidence_value(reliability, calibrated_confidence, coverage_rate, calibration, settings or DEFAULT_CONFIDENCE_GUARDRAIL)
    final_weight = reliability if effective_confidence is None else round(0.7 * reliability + 0.3 * effective_confidence, 4)
    record = {
        "tool": tool,
        "modality": modality,
        "sample_count": summary["sample_count"],
        "gene_rows": summary["gene_rows"],
        "callable_rate": summary["callable_rate"],
        "confidence_coverage_rate": coverage_rate,
        "base_reliability": reliability,
        "calibrated_confidence": "" if calibrated_confidence is None else calibrated_confidence,
        "effective_confidence": "" if effective_confidence is None else effective_confidence,
        "guardrail_factor": "" if guardrail_factor is None else guardrail_factor,
        "guardrail_status": guardrail_status,
        "brier_score": "" if calibration["brier_score"] is None else calibration["brier_score"],
        "expected_calibration_error": "" if calibration["expected_calibration_error"] is None else calibration["expected_calibration_error"],
        "final_weight": final_weight,
        "weight_version": WEIGHT_VERSION,
    }
    if gene is not None:
        record["gene"] = gene
    return record


def build_confidence_weights(rows, config=None):
    settings = confidence_guardrail_settings(config)
    grouped = defaultdict(list)
    gene_grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["modality"])].append(row)
        gene_grouped[(row["tool"], row["modality"], row["gene"])].append(row)
    tool_rows = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0])):
        tool_rows.append(compute_weight_record(entries, key[0], key[1], settings=settings))
    gene_rows = []
    for key, entries in sorted(gene_grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        gene_rows.append(compute_weight_record(entries, key[0], key[1], gene=key[2], settings=settings))
    return tool_rows, gene_rows


def build_runtime_weight_payload(tool_weights, gene_weights):
    payload = {
        "weight_version": WEIGHT_VERSION,
        "formula": {
            "final_weight": "0.7 * overall_correct_call_rate + 0.3 * effective_confidence_score",
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
            "calibrated_confidence": None if row["calibrated_confidence"] == "" else row["calibrated_confidence"],
            "effective_confidence": None if row["effective_confidence"] == "" else row["effective_confidence"],
            "guardrail_factor": None if row["guardrail_factor"] == "" else row["guardrail_factor"],
            "guardrail_status": row["guardrail_status"],
            "brier_score": None if row["brier_score"] == "" else row["brier_score"],
            "expected_calibration_error": None if row["expected_calibration_error"] == "" else row["expected_calibration_error"],
            "confidence_coverage_rate": row["confidence_coverage_rate"],
            "sample_count": row["sample_count"],
            "gene_rows": row["gene_rows"],
            "callable_rate": row["callable_rate"],
        }
    for row in gene_weights:
        payload["gene_weights"].setdefault(row["tool"], {}).setdefault(row["modality"], {})[row["gene"]] = {
            "final_weight": row["final_weight"],
            "base_reliability": row["base_reliability"],
            "calibrated_confidence": None if row["calibrated_confidence"] == "" else row["calibrated_confidence"],
            "effective_confidence": None if row["effective_confidence"] == "" else row["effective_confidence"],
            "guardrail_factor": None if row["guardrail_factor"] == "" else row["guardrail_factor"],
            "guardrail_status": row["guardrail_status"],
            "brier_score": None if row["brier_score"] == "" else row["brier_score"],
            "expected_calibration_error": None if row["expected_calibration_error"] == "" else row["expected_calibration_error"],
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
        return None, 0, 0, 0.0, 0.0
    ranked = sorted(pair_counts.items(), key=lambda item: (-len(item[1]), "%s|%s" % item[0]))
    top_pair, top_rows = ranked[0]
    second_n = len(ranked[1][1]) if len(ranked) > 1 else 0
    support_fraction = len(top_rows) / float(total_callable) if total_callable else 0.0
    support_margin = (len(top_rows) - second_n) / float(total_callable) if total_callable else 0.0
    is_tie = len(ranked) > 1 and len(ranked[1][1]) == len(top_rows)
    return (None if is_tie else top_pair), len(top_rows), total_callable, round(support_fraction, 4), round(support_margin, 4)


def lookup_weight(payload, tool, modality, gene):
    gene_entry = payload.get("gene_weights", {}).get(tool, {}).get(modality, {}).get(gene)
    if gene_entry and gene_entry.get("final_weight") is not None:
        return float(gene_entry["final_weight"])
    tool_entry = payload.get("tool_weights", {}).get(tool, {}).get(modality, {})
    if tool_entry.get("final_weight") is not None:
        return float(tool_entry["final_weight"])
    return 0.0


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
        pair, agreeing_tools, callable_tools, support_fraction, support_margin = majority_vote_for_group(entries)
        truth_1 = entries[0]["truth_allele1"]
        truth_2 = entries[0]["truth_allele2"]
        is_tie = pair is None and callable_tools > 0
        call_status = "called" if pair else "no_call"
        allele1, allele2 = pair if pair else ("", "")
        is_callable = "1" if pair else "0"
        ambiguity = evaluate_ambiguity([allele1, allele2], [truth_1, truth_2], 2)
        is_correct = "1" if pair and ambiguity["match_2field"] else "0"
        discordance_tag = "technical_conflict" if is_tie else ("no_evidence" if callable_tools == 0 else "consensus_call")
        out.append({
            "sample": key[0],
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
        return None, 0.0, 0, 0.0, 0.0, False
    ranked = sorted(pair_weights.items(), key=lambda item: (-item[1]["weight"], "%s|%s" % item[0]))
    top_pair, top_meta = ranked[0]
    second_weight = ranked[1][1]["weight"] if len(ranked) > 1 else 0.0
    top_support = top_meta["weight"] / total_weight if total_weight else 0.0
    margin = (top_meta["weight"] - second_weight) / total_weight if total_weight else 0.0
    is_tie = len(ranked) > 1 and abs(ranked[0][1]["weight"] - ranked[1][1]["weight"]) < 1e-9
    return top_pair, round(top_meta["weight"], 4), contributing_tools, round(top_support, 4), round(margin, 4), is_tie


def consensus_thresholds(config):
    benchmark_cfg = config.get("benchmark", {})
    consensus_cfg = benchmark_cfg.get("consensus", {})
    return {
        "min_support": float(consensus_cfg.get("min_support", 0.55)),
        "min_margin": float(consensus_cfg.get("min_margin", 0.15)),
    }


def build_weighted_consensus_rows(rows, weight_payload, config):
    thresholds = consensus_thresholds(config)
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["sample"], row["modality"], row["gene"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0], gene_sort_key(item[0][2]))):
        pair, total_weight, contributing_tools, support_fraction, support_margin, is_tie = weighted_vote_for_group(entries, weight_payload)
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
            "match_grade": ambiguity["match_grade"],
            "imgt_hla_version": entries[0].get("imgt_hla_version", ""),
            "agreeing_tools": "" if call_status != "called" else contributing_tools,
            "contributing_tools": contributing_tools,
            "support_fraction": as_string_number(support_fraction),
            "support_margin": as_string_number(support_margin),
            "total_weight": as_string_number(total_weight),
            "discordance_tag": discordance_tag,
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
        out.append({
            "method": key[0],
            "modality": key[1],
            "sample_count": len({row["sample"] for row in entries}),
            "gene_rows": total,
            "callable_rate": ratio(callable_n, total),
            "accuracy_among_callable": ratio(correct_n, callable_n),
            "overall_correct_call_rate": ratio(correct_n, total),
        })
    return out


def build_method_comparison(single_tool_summary, majority_rows, weighted_rows):
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
        })
    out.extend(dict(row, method_type="baseline") for row in summarize_consensus_method(majority_rows))
    out.extend(dict(row, method_type="ensemble") for row in summarize_consensus_method(weighted_rows))
    return out


def build_method_per_gene(single_tool_rows, majority_rows, weighted_rows):
    out = []
    out.extend({"method": row["tool"], "method_type": "single_tool", **{k: v for k, v in row.items() if k != "tool"}} for row in single_tool_rows)
    for consensus_rows, method_name, method_type in [(majority_rows, "MajorityVote", "baseline"), (weighted_rows, "WeightedConsensus", "ensemble")]:
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


def build_confidence_calibration_summary(rows, bin_count=5):
    grouped = defaultdict(list)
    for row in rows:
        score = coerce_float(row.get("confidence_score", ""))
        if score is None:
            continue
        grouped[(row["tool"], row["modality"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (modality_sort_key(item[0][1]), item[0][0])):
        stats = calibration_stats(entries, bin_count=bin_count)
        out.append({
            "tool": key[0],
            "modality": key[1],
            "n_rows": stats["n_rows"],
            "mean_confidence": stats["mean_confidence"],
            "observed_accuracy": stats["observed_accuracy"],
            "brier_score": stats["brier_score"],
            "expected_calibration_error": stats["expected_calibration_error"],
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


def build_discordance_rows(harmonized_rows, weighted_rows):
    weighted_index = {(row["sample"], row["modality"], row["gene"]): row for row in weighted_rows}
    grouped = defaultdict(list)
    for row in harmonized_rows:
        grouped[(row["sample"], row["gene"])].append(row)
    out = []
    for key, entries in sorted(grouped.items(), key=lambda item: (item[0][0], gene_sort_key(item[0][1]))):
        modality_pairs = {}
        for modality in sorted({row["modality"] for row in entries}, key=modality_sort_key):
            modality_entries = [row for row in entries if row["modality"] == modality]
            pair, agreeing_tools, callable_tools, support_fraction, support_margin = majority_vote_for_group(modality_entries)
            modality_pairs[modality] = {"pair": pair, "callable_tools": callable_tools, "support_fraction": support_fraction, "support_margin": support_margin}
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
            out.append({"sample": key[0], "gene": key[1], "scope": "cross_modality", "tag": "dna_rna_discordance", "detail": "rna_pair_differs_from_dna"})
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


def build_metadata(config, shared_genes, harmonized_rows, tool_weights):
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
        "truth_hierarchy": config.get("truth", {}).get("hierarchy", ["orthogonal_clinical_typing", "targeted_hla_ngs", "public_reference_truth"]),
        "benchmark_splits": benchmark_cfg.get("splits", {"training": "calibration", "validation": "threshold_tuning", "holdout": "final_reporting"}),
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
            "weight_version": WEIGHT_VERSION,
            "default_target_reads": DEFAULT_TARGET_READS,
            "formula": "0.7 * overall_correct_call_rate + 0.3 * effective_confidence_score",
            "effective_confidence_score": "base_reliability + shrink_factor * (mean_confidence_score - base_reliability)",
            "fallback": "overall_correct_call_rate when confidence is missing or blocked by guardrail",
            "guardrail": confidence_guardrail_settings(config),
            "parser_coverage": parser_coverage,
            "runtime_weight_override": benchmark_cfg.get("runtime_weight_override", ""),
            "tool_weight_count": len(tool_weights),
            "harmonized_rows_with_confidence": sum(1 for row in harmonized_rows if coerce_float(row.get("confidence_score", "")) is not None),
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
    config = load_yaml(Path(args.config))
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    truth = load_truth(config["truth"])
    runs = load_run_records(config)
    resolution = int(config.get("benchmark", {}).get("resolution", 2))
    imgt_hla_version = resolve_imgt_hla_version(config)
    harmonized_rows = dedupe_rows(build_harmonized_rows(truth, runs, resolution=resolution, imgt_hla_version=imgt_hla_version))
    if not harmonized_rows:
        raise SystemExit("No harmonized benchmark rows were produced. Check the config globs and parser settings.")
    benchmark_genes = resolve_benchmark_genes(config)
    filtered_rows = filter_rows(harmonized_rows, genes=benchmark_genes or None)
    if not filtered_rows:
        raise SystemExit("No benchmark rows remained after applying the configured benchmark loci.")
    shared_genes = calculate_shared_genes(filtered_rows)
    main_rows = filtered_rows
    summary = aggregate_metrics(main_rows)
    per_gene = build_per_gene_summary(main_rows)
    modality_gene = build_modality_gene_summary(main_rows)
    cohort = build_cohort_overview(main_rows, truth, benchmark_genes or shared_genes)
    ambiguity_summary_rows, ambiguity_summary_gene_rows = build_ambiguity_summary(main_rows)
    tool_weights, gene_weights = build_confidence_weights(main_rows, config=config)
    runtime_weights = load_runtime_weight_override(config) or build_runtime_weight_payload(tool_weights, gene_weights)
    majority_rows = build_majority_vote_rows(main_rows)
    weighted_rows = build_weighted_consensus_rows(main_rows, runtime_weights, config)
    method_ambiguity_rows, method_ambiguity_gene_rows = build_method_ambiguity_summary(majority_rows, weighted_rows)
    method_comparison_rows = build_method_comparison(summary, majority_rows, weighted_rows)
    method_per_gene_rows = build_method_per_gene(per_gene, majority_rows, weighted_rows)
    per_gene_gain_rows = build_per_gene_gain_table(method_per_gene_rows)
    calibration_bin_rows = build_confidence_bin_summary(main_rows)
    calibration_summary_rows = build_confidence_calibration_summary(main_rows)
    confidence_error_rows, confidence_error_gene_rows = build_confidence_error_summary(main_rows)
    abstention_rows = build_abstention_tradeoff(main_rows, runtime_weights, config)
    discordance_rows = build_discordance_rows(main_rows, weighted_rows)
    discordance_summary_rows = summarize_discordance(discordance_rows)
    metadata = build_metadata(config, shared_genes, harmonized_rows, tool_weights)

    write_tsv(output_dir / "tables" / "harmonized_benchmark_rows.tsv", main_rows, ["sample", "modality", "tool", "gene", "truth_allele1_raw", "truth_allele2_raw", "allele1_raw", "allele2_raw", "truth_allele1", "truth_allele2", "allele1", "allele2", "truth_allele1_3field", "truth_allele2_3field", "allele1_3field", "allele2_3field", "call_status", "correct_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_correct_g_group", "is_correct_p_group", "match_grade", "imgt_hla_version", "runtime_hours", "max_ram_gb", "confidence_score", "confidence_source", "raw_confidence", "read_support", "source_file"])
    write_tsv(output_dir / "tables" / "summary_full_cohort.tsv", sorted(summary.values(), key=lambda row: (modality_sort_key(row["modality"]), row["tool"])), ["tool", "modality", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate", "median_runtime_hours", "median_max_ram_gb"])
    write_tsv(output_dir / "tables" / "summary_per_gene.tsv", per_gene, ["tool", "modality", "gene", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "summary_modality_gene_coverage.tsv", modality_gene, ["modality", "gene"])
    write_tsv(output_dir / "tables" / "cohort_overview.tsv", cohort, ["modality", "tool", "truth_samples", "completed_tool_samples", "modality_samples", "shared_gene_count"])
    write_tsv(output_dir / "tables" / "ambiguity_summary.tsv", ambiguity_summary_rows, ["tool", "modality", "imgt_hla_version", "sample_count", "gene_rows", "exact_2field_rate", "exact_3field_rate", "g_group_comparable_rows", "g_group_match_rate", "p_group_comparable_rows", "p_group_match_rate"])
    write_tsv(output_dir / "tables" / "ambiguity_summary_by_gene.tsv", ambiguity_summary_gene_rows, ["tool", "modality", "gene", "imgt_hla_version", "sample_count", "gene_rows", "exact_2field_rate", "exact_3field_rate", "g_group_comparable_rows", "g_group_match_rate", "p_group_comparable_rows", "p_group_match_rate"])
    write_tsv(output_dir / "tables" / "reference_metadata.tsv", [{"truth_source": metadata.get("truth_source", ""), "truth_path": metadata.get("truth_path", ""), "imgt_hla_version": metadata.get("imgt_hla_version", ""), "primary_resolution": metadata.get("resolution", ""), "secondary_resolutions": ",".join(str(x) for x in metadata.get("ambiguity_evaluation", {}).get("secondary_resolutions", []))}], ["truth_source", "truth_path", "imgt_hla_version", "primary_resolution", "secondary_resolutions"])
    write_tsv(output_dir / "tables" / "sample_level_disagreements.tsv", build_sample_disagreements(harmonized_rows), ["sample", "modality", "tool", "gene", "truth_allele1", "truth_allele2", "allele1", "allele2", "call_status", "correct_status", "source_file"])
    write_tsv(output_dir / "tables" / "missing_call_patterns.tsv", build_missing_patterns(harmonized_rows), ["tool", "modality", "gene", "missing_calls", "total_rows", "missing_rate"])
    write_tsv(output_dir / "tables" / "truth_mismatches.tsv", build_truth_mismatches(harmonized_rows), ["sample", "modality", "tool", "gene", "truth_allele1", "truth_allele2", "typed_allele1", "typed_allele2", "call_status"])
    write_tsv(output_dir / "tables" / "tool_confidence_weights.tsv", tool_weights, ["tool", "modality", "sample_count", "gene_rows", "callable_rate", "confidence_coverage_rate", "base_reliability", "calibrated_confidence", "effective_confidence", "guardrail_factor", "guardrail_status", "brier_score", "expected_calibration_error", "final_weight", "weight_version"])
    write_tsv(output_dir / "tables" / "tool_confidence_weights_by_gene.tsv", gene_weights, ["tool", "modality", "gene", "sample_count", "gene_rows", "callable_rate", "confidence_coverage_rate", "base_reliability", "calibrated_confidence", "effective_confidence", "guardrail_factor", "guardrail_status", "brier_score", "expected_calibration_error", "final_weight", "weight_version"])
    write_tsv(output_dir / "tables" / "majority_vote_baseline.tsv", majority_rows, ["sample", "modality", "gene", "method", "truth_allele1", "truth_allele2", "allele1", "allele2", "truth_allele1_3field", "truth_allele2_3field", "allele1_3field", "allele2_3field", "call_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_correct_g_group", "is_correct_p_group", "match_grade", "imgt_hla_version", "agreeing_tools", "contributing_tools", "support_fraction", "support_margin", "discordance_tag"])
    write_tsv(output_dir / "tables" / "weighted_consensus_calls.tsv", weighted_rows, ["sample", "modality", "gene", "method", "truth_allele1", "truth_allele2", "allele1", "allele2", "truth_allele1_3field", "truth_allele2_3field", "allele1_3field", "allele2_3field", "call_status", "is_callable", "is_correct", "is_correct_2field", "is_correct_3field", "is_correct_g_group", "is_correct_p_group", "match_grade", "imgt_hla_version", "agreeing_tools", "contributing_tools", "support_fraction", "support_margin", "total_weight", "discordance_tag"])
    write_tsv(output_dir / "tables" / "method_comparison.tsv", method_comparison_rows, ["method", "method_type", "modality", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "method_per_gene.tsv", method_per_gene_rows, ["method", "method_type", "modality", "gene", "sample_count", "gene_rows", "callable_rate", "accuracy_among_callable", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "method_ambiguity_summary.tsv", method_ambiguity_rows, ["method", "modality", "imgt_hla_version", "sample_count", "gene_rows", "exact_2field_rate", "exact_3field_rate", "g_group_comparable_rows", "g_group_match_rate", "p_group_comparable_rows", "p_group_match_rate"])
    write_tsv(output_dir / "tables" / "method_ambiguity_summary_by_gene.tsv", method_ambiguity_gene_rows, ["method", "modality", "gene", "imgt_hla_version", "sample_count", "gene_rows", "exact_2field_rate", "exact_3field_rate", "g_group_comparable_rows", "g_group_match_rate", "p_group_comparable_rows", "p_group_match_rate"])
    write_tsv(output_dir / "tables" / "per_gene_gain.tsv", per_gene_gain_rows, ["modality", "gene", "weighted_correct_call_rate", "majority_correct_call_rate", "best_single_tool", "best_single_tool_rate", "gain_vs_majority", "gain_vs_best_single"])
    write_tsv(output_dir / "tables" / "confidence_bin_summary.tsv", calibration_bin_rows, ["tool", "modality", "bin_index", "bin_lower", "bin_upper", "n_rows", "mean_confidence", "observed_accuracy"])
    write_tsv(output_dir / "tables" / "confidence_calibration_summary.tsv", calibration_summary_rows, ["tool", "modality", "n_rows", "mean_confidence", "observed_accuracy", "brier_score", "expected_calibration_error"])
    write_tsv(output_dir / "tables" / "confidence_error_summary.tsv", confidence_error_rows, ["tool", "modality", "bin_index", "bin_lower", "bin_upper", "n_rows", "error_count", "error_rate", "observed_accuracy", "mean_confidence"])
    write_tsv(output_dir / "tables" / "confidence_error_summary_by_gene.tsv", confidence_error_gene_rows, ["tool", "modality", "gene", "bin_index", "bin_lower", "bin_upper", "n_rows", "error_count", "error_rate", "observed_accuracy", "mean_confidence"])
    write_tsv(output_dir / "tables" / "abstention_tradeoff.tsv", abstention_rows, ["min_support", "min_margin", "call_rate", "no_call_rate", "low_confidence_rate", "accuracy_among_called", "overall_correct_call_rate"])
    write_tsv(output_dir / "tables" / "discordance_tags.tsv", discordance_rows, ["sample", "gene", "scope", "tag", "detail"])
    write_tsv(output_dir / "tables" / "discordance_summary.tsv", discordance_summary_rows, ["scope", "tag", "n_events"])
    save_json(output_dir / "tables" / "benchmark_metadata.json", metadata)
    save_json(output_dir / "tables" / "consensus_runtime_weights.json", runtime_weights)
    generate_figures(output_dir, method_comparison_rows, per_gene_gain_rows, calibration_bin_rows, abstention_rows, discordance_summary_rows, tool_weights)
    write_caption_stub(output_dir / "figures" / "captions.md", metadata)
    return 0


if __name__ == "__main__":
    sys.exit(main())
