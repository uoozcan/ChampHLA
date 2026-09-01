from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from champhla_confirmation.io import canonical_allele

from .io import sha256, write_json, write_tsv


GENES = ("A", "B", "C")
MODALITIES = ("wgs", "wes", "rnaseq")


def _unique_two_field(value: str, gene: str) -> str:
    options = set()
    for raw in value.strip().split("/"):
        token = raw.strip()
        if not token:
            continue
        if "*" not in token:
            token = f"{gene}*{token}"
        options.add(canonical_allele(token, gene))
    if len(options) != 1:
        raise ValueError(f"ambiguous two-field allele {gene} {value!r}")
    return next(iter(options))


def _read_2014(path: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter=" ", skipinitialspace=True))


def prepare_locked_truth(source: str, truth_output: str, registry_output: str,
                         manifest_output: str) -> dict:
    source_hash = sha256(source)
    parsed: dict[str, list[tuple[dict[str, tuple[str, str]], str]]] = defaultdict(list)
    for row in _read_2014(source):
        subject = row.get("id", "").strip().upper()
        if not subject:
            continue
        pairs = {}
        for gene in GENES:
            try:
                pair = sorted((_unique_two_field(row[gene], gene),
                               _unique_two_field(row[f"{gene}.1"], gene)))
                pairs[gene] = (pair[0], pair[1])
            except (KeyError, ValueError):
                continue
        parsed[subject].append((pairs, row.get("sbgroup", "").strip()))

    truth_rows = []
    registry_rows = []
    conflicts = 0
    for subject in sorted(parsed):
        resolved = {}
        subject_conflict = False
        for gene in GENES:
            unique = {entry[0][gene] for entry in parsed[subject] if gene in entry[0]}
            if len(unique) == 1:
                resolved[gene] = next(iter(unique))
            elif len(unique) > 1:
                subject_conflict = True
        if subject_conflict:
            conflicts += 1
        eligible = bool(resolved)
        complete = len(resolved) == len(GENES)
        population = next((entry[1] for entry in parsed[subject] if entry[1]), "")
        registry_rows.append({
            "subject": subject,
            "population": population,
            "primary_eligible": int(eligible),
            "complete_truth_available": int(complete),
            "exact_truth_loci": len(resolved),
            "source_record_count": len(parsed[subject]),
            "eligibility_reason": "" if eligible else "duplicate_conflict" if subject_conflict else "unresolved_or_ambiguous",
            "truth_source_sha256": source_hash,
        })
        if not eligible:
            continue
        for modality in MODALITIES:
            for gene in sorted(resolved):
                truth_rows.append({
                    "cohort": "1000G-2014-unseen",
                    "subject": subject,
                    "modality": modality,
                    "gene": gene,
                    "truth_allele1": resolved[gene][0],
                    "truth_allele2": resolved[gene][1],
                    "truth_status": "callable",
                    "truth_source": "2014 1000 Genomes laboratory HLA typing",
                    "source_sha256": source_hash,
                })

    write_tsv(truth_output, truth_rows, [
        "cohort", "subject", "modality", "gene", "truth_allele1", "truth_allele2",
        "truth_status", "truth_source", "source_sha256",
    ])
    write_tsv(registry_output, registry_rows, [
        "subject", "population", "primary_eligible", "complete_truth_available",
        "exact_truth_loci", "source_record_count", "eligibility_reason", "truth_source_sha256",
    ])
    result = {
        "schema_version": "locked-1000g-truth-1",
        "source_sha256": source_hash,
        "source_subjects": len(parsed),
        "eligible_any_exact_locus_subjects": sum(int(row["primary_eligible"]) for row in registry_rows),
        "eligible_complete_subjects": sum(int(row["complete_truth_available"]) for row in registry_rows),
        "duplicate_conflicts": conflicts,
        "truth_rows": len(truth_rows),
        "truth_output_sha256": sha256(truth_output),
        "registry_output_sha256": sha256(registry_output),
        "primary_truth": True,
        "inferred_2018_truth_allowed": False,
    }
    write_json(manifest_output, result)
    return result
