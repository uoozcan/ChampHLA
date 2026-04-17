#!/usr/bin/env python3
"""Build canonical 1000 Genomes HLA benchmark manifests."""

import argparse
import csv
import random
from collections import defaultdict
from pathlib import Path

import hla_benchmark as hb

REQUIRED_MODALITIES = ["wgs", "wes", "rnaseq"]


def parse_args():
    parser = argparse.ArgumentParser(description="Build 1000 Genomes benchmark manifests.")
    parser.add_argument("--truth", required=True, help="Truth TSV from the 1000 Genomes HLA reference set.")
    parser.add_argument("--sequencing", required=True, help="Long TSV listing available WGS/WES/RNA-seq data by sample.")
    parser.add_argument("--output-dir", required=True, help="Output directory for truth/sequencing/cohort manifests.")
    parser.add_argument("--truth-source", default="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140725_hla_genotypes/", help="Source URI for the truth set.")
    parser.add_argument("--acquisition-date", default="", help="Acquisition date recorded in the manifest metadata.")
    parser.add_argument("--supported-loci", default="A,B,C,DRB1,DQB1", help="Comma-separated loci for formal benchmarking.")
    parser.add_argument("--population-manifest", default="", help="Optional TSV mapping sample to population.")
    parser.add_argument("--split-seed", type=int, default=1000, help="Seed for deterministic split assignment.")
    return parser.parse_args()


def parse_population_map(path_value):
    if not path_value:
        return {}
    rows = hb.read_table(Path(path_value))
    if not rows:
        return {}
    key_map = {hb.normalize_column_name(name): name for name in rows[0].keys()}
    sample_key = hb.first_existing(key_map, ["sample", "sampleid", "id"])
    population_key = hb.first_existing(key_map, ["population", "pop", "superpopulation", "superpop"])
    if not (sample_key and population_key):
        return {}
    out = {}
    for row in rows:
        sample = hb.clean_token(row[key_map[sample_key]])
        if sample:
            out[sample] = hb.clean_token(row[key_map[population_key]])
    return out


def load_truth_manifest_rows(truth_path, supported_loci, truth_source, acquisition_date, population_map):
    truth = hb.load_truth({"path": truth_path})
    rows = []
    for sample in sorted(truth):
        genes = sorted([gene for gene in truth[sample] if gene in supported_loci], key=hb.gene_sort_key)
        rows.append({
            "sample": sample,
            "population": population_map.get(sample, ""),
            "truth_source": truth_source,
            "acquisition_date": acquisition_date,
            "truth_supported_loci": ",".join(genes),
            "truth_gene_count": len(genes),
        })
    return rows


def load_sequencing_manifest_rows(path_value, population_map):
    rows = hb.read_table(Path(path_value))
    if not rows:
        return []
    key_map = {hb.normalize_column_name(name): name for name in rows[0].keys()}
    sample_key = hb.first_existing(key_map, ["sample", "sampleid", "id"])
    modality_key = hb.first_existing(key_map, ["modality", "datatype", "assay"])
    locator_key = hb.first_existing(key_map, ["datalocator", "path", "url", "source", "publicurl", "resolvedpath"])
    available_key = hb.first_existing(key_map, ["available", "availability", "present"])
    population_key = hb.first_existing(key_map, ["population", "pop", "superpopulation", "superpop"])
    if not (sample_key and modality_key and locator_key):
        raise ValueError("Sequencing manifest source must include sample, modality, and data locator columns")
    normalized = []
    for row in rows:
        sample = hb.clean_token(row[key_map[sample_key]])
        modality = hb.clean_token(row[key_map[modality_key]]).lower()
        locator = hb.clean_token(row[key_map[locator_key]])
        if not sample or modality not in REQUIRED_MODALITIES:
            continue
        available_token = hb.clean_token(row.get(key_map[available_key], "1") if available_key else "1").lower()
        available = "1" if available_token not in {"0", "false", "no", "missing"} and locator else "0"
        population = hb.clean_token(row.get(key_map[population_key], "")) if population_key else population_map.get(sample, "")
        normalized.append({
            "sample": sample,
            "population": population or population_map.get(sample, ""),
            "modality": modality,
            "data_locator": locator,
            "available": available,
        })
    return sorted(normalized, key=lambda row: (row["sample"], hb.modality_sort_key(row["modality"])))


def assign_split_group(samples, seed):
    samples = list(samples)
    random.Random(seed).shuffle(samples)
    n = len(samples)
    if n <= 1:
        return {sample: "holdout" for sample in samples}
    if n == 2:
        return {samples[0]: "training", samples[1]: "holdout"}
    if n == 3:
        return {samples[0]: "training", samples[1]: "validation", samples[2]: "holdout"}
    n_train = max(1, int(round(n * 0.6)))
    n_val = max(1, int(round(n * 0.2)))
    if n_train + n_val >= n:
        n_train = max(1, n - 2)
        n_val = 1
    assignments = {}
    for idx, sample in enumerate(samples):
        if idx < n_train:
            assignments[sample] = "training"
        elif idx < n_train + n_val:
            assignments[sample] = "validation"
        else:
            assignments[sample] = "holdout"
    return assignments


def build_cohort_manifest(truth_manifest_rows, sequencing_manifest_rows, supported_loci, split_seed, required_modalities=None):
    if required_modalities is None:
        required_modalities = REQUIRED_MODALITIES
    truth_index = {row["sample"]: row for row in truth_manifest_rows}
    sequencing_by_sample = defaultdict(dict)
    for row in sequencing_manifest_rows:
        sequencing_by_sample[row["sample"]][row["modality"]] = row
    included_samples = []
    rows = []
    for sample in sorted(truth_index):
        truth_row = truth_index[sample]
        loci = [gene for gene in hb.parse_gene_list(truth_row.get("truth_supported_loci", "")) if gene in supported_loci]
        modalities = sequencing_by_sample.get(sample, {})
        available = {modality: modalities.get(modality, {}).get("available", "0") == "1" for modality in REQUIRED_MODALITIES}
        missing_modalities = [modality for modality in required_modalities if not available.get(modality)]
        include = not missing_modalities and bool(loci)
        reason = "" if include else ("missing_modality:" + ",".join(missing_modalities) if missing_modalities else "no_supported_truth_loci")
        rows.append({
            "sample": sample,
            "population": truth_row.get("population", "") or next((modalities[m].get("population", "") for m in REQUIRED_MODALITIES if m in modalities), ""),
            "include": "1" if include else "0",
            "split": "",
            "truth_supported_loci": ",".join(loci),
            "wgs_available": "1" if available.get("wgs") else "0",
            "wes_available": "1" if available.get("wes") else "0",
            "rnaseq_available": "1" if available.get("rnaseq") else "0",
            "excluded_reason": reason,
        })
        if include:
            included_samples.append(sample)
    by_population = defaultdict(list)
    for row in rows:
        if row["include"] == "1":
            by_population[row["population"] or "unknown"].append(row["sample"])
    split_map = {}
    if by_population and max(len(samples) for samples in by_population.values()) >= 2:
        for offset, population in enumerate(sorted(by_population)):
            split_map.update(assign_split_group(sorted(by_population[population]), split_seed + offset))
    else:
        split_map.update(assign_split_group(sorted([row["sample"] for row in rows if row["include"] == "1"]), split_seed))
    for row in rows:
        if row["include"] == "1":
            row["split"] = split_map.get(row["sample"], "holdout")
    return rows


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    supported_loci = hb.parse_gene_list(args.supported_loci)
    population_map = parse_population_map(args.population_manifest)
    truth_manifest_rows = load_truth_manifest_rows(args.truth, supported_loci, args.truth_source, args.acquisition_date, population_map)
    sequencing_manifest_rows = load_sequencing_manifest_rows(args.sequencing, population_map)
    cohort_manifest_rows = build_cohort_manifest(truth_manifest_rows, sequencing_manifest_rows, supported_loci, args.split_seed)
    write_tsv(output_dir / "truth_manifest.tsv", truth_manifest_rows, ["sample", "population", "truth_source", "acquisition_date", "truth_supported_loci", "truth_gene_count"])
    write_tsv(output_dir / "sequencing_manifest.tsv", sequencing_manifest_rows, ["sample", "population", "modality", "data_locator", "available"])
    write_tsv(output_dir / "cohort_manifest.tsv", cohort_manifest_rows, ["sample", "population", "include", "split", "truth_supported_loci", "wgs_available", "wes_available", "rnaseq_available", "excluded_reason"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
