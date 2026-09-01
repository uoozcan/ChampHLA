from __future__ import annotations

import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path

from .io import read_tsv, sha256, write_json, write_tsv


def _read_csv(path: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle))


def _split(value: str) -> set[str]:
    return {token.strip().upper() for token in value.replace(";", ",").split(",") if token.strip()}


def _metadata_tokens(row: dict[str, str]) -> set[str]:
    tokens = set()
    for field in ("sample_id", "alternative_id", "biosample_id", "paternal_id", "maternal_id"):
        tokens.update(_split(row.get(field, "")))
    return tokens


def build_hprc_release2_candidates(
    assembly_csv: str,
    illumina_csv: str,
    sample_metadata_csv: str,
    development_tsv: str,
    historical_exclusions: str,
    output_tsv: str,
    summary_json: str,
    source_commit: str = "",
) -> dict:
    """Build a truth-free HPRC R2 candidate roster from official metadata.

    Eligibility requires both assembly haplotypes and paired Illumina WGS. Exact
    identifiers, aliases, explicit parents, and shared family IDs with a
    development subject are excluded before any truth is available.
    """
    assemblies = _read_csv(assembly_csv)
    illumina = _read_csv(illumina_csv)
    metadata = _read_csv(sample_metadata_csv)
    development = read_tsv(development_tsv)

    dev_ids = set()
    for row in development:
        for field in ("subject", "sample", "sample_id", "alias", "aliases", "coriell_id",
                      "biosample", "biosample_accession", "individual_id"):
            dev_ids.update(_split(row.get(field, "")))
    with Path(historical_exclusions).open(encoding="utf-8", errors="replace") as handle:
        dev_ids.update(line.strip().upper() for line in handle
                       if line.strip() and not line.lstrip().startswith("#"))

    meta_by_sample = {row["sample_id"].strip().upper(): row for row in metadata
                      if row.get("sample_id", "").strip()}
    dev_family_ids = {
        row.get("family_id", "").strip().upper()
        for row in metadata
        if row.get("family_id", "").strip()
        and (_metadata_tokens(row) & dev_ids)
    }

    assembly_by_sample: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in assemblies:
        sample = row.get("sample_id", "").strip().upper()
        haplotype = row.get("haplotype", "").strip()
        if sample and haplotype in {"1", "2"}:
            assembly_by_sample[sample][haplotype] = row

    illumina_by_sample: dict[str, dict[str, str]] = {}
    for row in illumina:
        sample = row.get("sample_id", "").strip().upper()
        valid = (row.get("filetype", "").strip().lower() == "cram"
                 and row.get("library_strategy", "").strip().upper() == "WGS"
                 and row.get("library_layout", "").strip().upper() == "PAIRED")
        if not sample or not valid:
            continue
        try:
            coverage = float(row.get("coverage", "0") or 0)
            old_coverage = float(illumina_by_sample.get(sample, {}).get("coverage", "0") or 0)
        except ValueError:
            coverage, old_coverage = 0.0, 0.0
        if sample not in illumina_by_sample or coverage > old_coverage:
            illumina_by_sample[sample] = row

    rows = []
    all_samples = sorted(set(assembly_by_sample) | set(illumina_by_sample))
    for sample in all_samples:
        meta = meta_by_sample.get(sample, {})
        haps = assembly_by_sample.get(sample, {})
        seq = illumina_by_sample.get(sample)
        aliases = _metadata_tokens(meta) | {sample}
        exact_overlap = bool(aliases & dev_ids)
        family_id = meta.get("family_id", "").strip().upper()
        explicit_relatives = (_split(meta.get("paternal_id", ""))
                              | _split(meta.get("maternal_id", ""))
                              | _split(meta.get("siblings", "")))
        relative_overlap = bool((family_id and family_id in dev_family_ids)
                                or (explicit_relatives & dev_ids))
        complete_assemblies = set(haps) == {"1", "2"}
        paired_wgs = seq is not None
        eligible = complete_assemblies and paired_wgs and not exact_overlap and not relative_overlap
        reasons = []
        if not complete_assemblies:
            reasons.append("missing_phased_assembly")
        if not paired_wgs:
            reasons.append("missing_paired_illumina_wgs")
        if exact_overlap:
            reasons.append("development_identifier_or_alias_overlap")
        if relative_overlap:
            reasons.append("known_development_relative")
        cram_path = seq.get("path", "") if seq else ""
        rows.append({
            "cohort": "HPRC_R2", "subject": sample, "sample_id": sample,
            "aliases": ",".join(sorted(aliases - {sample})),
            "coriell_id": meta.get("alternative_id", "") or sample,
            "biosample_accession": meta.get("biosample_id", ""),
            "individual_id": sample, "family_id": family_id,
            "paternal_id": meta.get("paternal_id", ""),
            "maternal_id": meta.get("maternal_id", ""),
            "population": meta.get("population_abbreviation", ""),
            "development_relative": int(relative_overlap),
            "input_type": "cram" if seq else "", "input_url": cram_path,
            "index_url": f"{cram_path}.crai" if cram_path else "",
            "genome_build": "GRCh38DH_header_verification_required",
            "paired_end": int(paired_wgs), "read_length": seq.get("read_length", "") if seq else "",
            "coverage": seq.get("coverage", "") if seq else "",
            "assembly_hap1": haps.get("1", {}).get("assembly", ""),
            "assembly_hap2": haps.get("2", {}).get("assembly", ""),
            "truth_source": "two independent blinded assembly extraction/alignment procedures",
            "eligible_nonoverlap": int(eligible), "exclusion_reason": ";".join(reasons),
        })

    fields = [
        "cohort", "subject", "sample_id", "aliases", "coriell_id", "biosample_accession",
        "individual_id", "family_id", "paternal_id", "maternal_id", "population",
        "development_relative", "input_type", "input_url", "index_url", "genome_build",
        "paired_end", "read_length", "coverage", "assembly_hap1", "assembly_hap2",
        "truth_source", "eligible_nonoverlap", "exclusion_reason",
    ]
    write_tsv(output_tsv, rows, fields)
    eligible_count = sum(row["eligible_nonoverlap"] for row in rows)
    summary = {
        "schema_version": "hprc-release2-candidates-1",
        "truth_blind": True,
        "source_commit": source_commit,
        "source_sha256": {
            "assemblies": sha256(assembly_csv),
            "illumina": sha256(illumina_csv),
            "sample_metadata": sha256(sample_metadata_csv),
            "development": sha256(development_tsv),
            "historical_exclusions": sha256(historical_exclusions),
        },
        "official_subjects_seen": len(rows),
        "subjects_with_two_haplotypes": sum(set(assembly_by_sample.get(s, {})) == {"1", "2"}
                                            for s in all_samples),
        "subjects_with_paired_wgs": len(illumina_by_sample),
        "eligible_nonoverlap": eligible_count,
        "target_120_available": eligible_count >= 120,
        "maximum_180_available": eligible_count >= 180,
        "prediction_or_truth_used": False,
    }
    write_json(summary_json, summary)
    return summary


def select_hprc_confirmation_roster(candidate_tsv: str, target: int, seed: str,
                                    output_tsv: str, summary_json: str) -> dict:
    """Select a deterministic population-proportional roster without truth."""
    if target < 1:
        raise ValueError("target must be positive")
    source = read_tsv(candidate_tsv)
    eligible = [row for row in source if row.get("eligible_nonoverlap") == "1"]
    if len(eligible) < target:
        raise ValueError(f"only {len(eligible)} eligible HPRC subjects for target {target}")
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in eligible:
        groups[row.get("population", "").strip() or "UNSPECIFIED"].append(row)

    ideal = {group: target * len(rows) / len(eligible) for group, rows in groups.items()}
    quota = {group: min(len(groups[group]), math.floor(value)) for group, value in ideal.items()}
    remaining = target - sum(quota.values())
    order = sorted(groups, key=lambda group: (-(ideal[group] - quota[group]), group))
    for group in order:
        if remaining == 0:
            break
        if quota[group] < len(groups[group]):
            quota[group] += 1
            remaining -= 1
    if remaining:
        raise ValueError("population quota allocation failed")

    selected = []
    for group, rows in groups.items():
        ranked = sorted(rows, key=lambda row: (
            hashlib.sha256(f"{seed}|{row['subject']}".encode()).hexdigest(), row["subject"]
        ))
        selected.extend(ranked[:quota[group]])
    selected.sort(key=lambda row: (row.get("population", ""), row["subject"]))
    for index, row in enumerate(selected, start=1):
        row["selection_rank"] = index
        row["selection_seed"] = seed
        row["selection_rule"] = "population-proportional largest-remainder; seeded SHA-256 within stratum"
    fields = list(source[0]) + ["selection_rank", "selection_seed", "selection_rule"]
    write_tsv(output_tsv, selected, fields)
    summary = {
        "schema_version": "hprc-confirmation-selection-1",
        "truth_blind": True,
        "candidate_sha256": sha256(candidate_tsv),
        "target": target,
        "selected": len(selected),
        "seed": seed,
        "selection_rule": "population-proportional largest-remainder; seeded SHA-256 within stratum",
        "population_counts": dict(sorted((group, quota[group]) for group in quota)),
        "prediction_or_truth_used": False,
    }
    write_json(summary_json, summary)
    return summary
