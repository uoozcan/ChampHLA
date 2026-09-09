from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

from .io import read_tsv, reject_truth_columns, sha256, write_json


DATASET_FIELDS = {
    "dataset_id", "cohort", "repository", "study_accession", "file_accessions",
    "assay", "organism", "reference_build", "donor_count", "library_count",
    "read_layout", "read_length", "platform", "library_prep", "tissue",
    "raw_data_status", "access_class", "estimated_size_gb", "publication_url",
    "evidence_url", "independence_stratum", "sample_link_status",
    "capacity_minimum", "capacity_met", "audit_status", "audit_notes",
    "access_terms", "tool_compatibility",
}
TRUTH_FIELDS = {
    "truth_id", "cohort", "method", "source_material", "loci", "resolution",
    "ambiguity_status", "typing_date", "reconciliation", "rna_informed",
    "truth_independence", "circularity", "accuracy_eligible", "reported_subjects",
    "usable_truth_count", "access_class", "source_accession", "evidence_url",
    "audit_status", "audit_notes",
}
CROSSWALK_FIELDS = {
    "cohort", "repository_sample", "donor", "aliases", "modality",
    "sequencing_accession", "biosample_accession", "experiment_accession",
    "file_accessions", "source_urls", "truth_identifier", "sample_link_status",
    "existing_champhla_overlap", "relatedness", "independence_stratum",
    "access_class", "gate_status", "exclusion_reason",
}
PILOT_FIELDS = {
    "cohort", "subject", "modality", "gene", "source", "source_accession",
    "input_type", "read1_url", "read2_url", "access_class",
    "independence_stratum", "gate_status", "exclusion_reason",
}

AUDIT_STATUSES = {"verified", "conditional", "rejected"}
ASSAYS = {"wgs", "wes", "rnaseq", "targeted_mhc", "scrnaseq", "truth_only"}
ACCESS_CLASSES = {"public", "controlled", "mixed", "unresolved"}
INDEPENDENCE_STRATA = {"donor_independent", "new_library_overlap", "not_applicable", "unresolved"}
PILOT_GATES = {
    "ready", "blocked_raw_accession", "blocked_truth_access",
    "blocked_truth_crosswalk", "blocked_truth_ambiguity", "excluded_related",
}
PILOT_SUBJECT_TARGET = 12
LANE_STATES = (
    "inventory",
    "metadata_verified",
    "truth_verified",
    "executable",
    "prediction_frozen",
    "evaluated",
)
PILOT_COHORT_ALIASES = {
    "NCI60_PUBLIC_PILOT": "NCI60",
    "AFGR_MKK_GATED_PILOT": "AFGR_MKK",
}


def _require_header(rows: list[dict[str, str]], fields: set[str], context: str) -> None:
    if not rows:
        raise ValueError(f"{context} is empty")
    missing = sorted(fields - set(rows[0]))
    if missing:
        raise ValueError(f"{context} missing columns: {missing}")


def _integer(value: str, field: str, context: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise ValueError(f"{context} has invalid {field}: {value!r}") from error
    if parsed < 0:
        raise ValueError(f"{context} has negative {field}: {parsed}")
    return parsed


def _yes(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def _split(value: str) -> list[str]:
    return [token.strip() for token in value.replace(";", ",").split(",") if token.strip()]


def _validate_datasets(rows: list[dict[str, str]]) -> None:
    _require_header(rows, DATASET_FIELDS, "dataset registry")
    seen = set()
    for row in rows:
        dataset_id = row["dataset_id"].strip()
        if not dataset_id or dataset_id in seen:
            raise ValueError(f"duplicate or empty dataset_id: {dataset_id!r}")
        seen.add(dataset_id)
        context = f"dataset {dataset_id}"
        if row["audit_status"] not in AUDIT_STATUSES:
            raise ValueError(f"{context} has invalid audit_status")
        if row["assay"] not in ASSAYS:
            raise ValueError(f"{context} has invalid assay")
        if row["access_class"] not in ACCESS_CLASSES:
            raise ValueError(f"{context} has invalid access_class")
        if row["independence_stratum"] not in INDEPENDENCE_STRATA:
            raise ValueError(f"{context} has invalid independence_stratum")
        mandatory = ("access_terms", "tool_compatibility")
        missing = [field for field in mandatory if not row[field].strip()]
        if missing:
            raise ValueError(f"{context} missing mandatory routing fields: {missing}")
        donors = _integer(row["donor_count"], "donor_count", context)
        _integer(row["library_count"], "library_count", context)
        minimum = _integer(row["capacity_minimum"], "capacity_minimum", context)
        if _yes(row["capacity_met"]) != (donors >= minimum if minimum else False):
            raise ValueError(f"{context} capacity_met disagrees with donor_count/minimum")
        if row["audit_status"] == "verified":
            required = ("repository", "study_accession", "file_accessions", "evidence_url",
                        "raw_data_status", "read_layout", "reference_build", "sample_link_status")
            missing = [field for field in required if not row[field].strip()]
            if missing:
                raise ValueError(f"{context} verified without {missing}")


def _validate_truth(rows: list[dict[str, str]]) -> None:
    _require_header(rows, TRUTH_FIELDS, "truth registry")
    seen = set()
    for row in rows:
        truth_id = row["truth_id"].strip()
        if not truth_id or truth_id in seen:
            raise ValueError(f"duplicate or empty truth_id: {truth_id!r}")
        seen.add(truth_id)
        context = f"truth {truth_id}"
        if row["audit_status"] not in AUDIT_STATUSES:
            raise ValueError(f"{context} has invalid audit_status")
        if row["access_class"] not in ACCESS_CLASSES:
            raise ValueError(f"{context} has invalid access_class")
        if row["circularity"] not in {"orthogonal", "partly_circular", "circular", "unresolved"}:
            raise ValueError(f"{context} has invalid circularity")
        _integer(row["reported_subjects"], "reported_subjects", context)
        usable = _integer(row["usable_truth_count"], "usable_truth_count", context)
        if _yes(row["accuracy_eligible"]) and (row["circularity"] != "orthogonal" or usable == 0):
            raise ValueError(f"{context} cannot be accuracy eligible")
        if _yes(row["accuracy_eligible"]):
            required = ("method", "source_material", "loci", "resolution", "source_accession",
                        "evidence_url", "truth_independence")
            missing = [field for field in required if not row[field].strip()]
            if missing:
                raise ValueError(f"{context} eligible without provenance fields: {missing}")


def _validate_crosswalk(rows: list[dict[str, str]]) -> None:
    _require_header(rows, CROSSWALK_FIELDS, "sample crosswalk")
    seen = set()
    for row in rows:
        key = (row["cohort"], row["repository_sample"], row["modality"],
               row["sequencing_accession"])
        if key in seen:
            raise ValueError(f"duplicate crosswalk row: {key}")
        seen.add(key)
        if row["modality"] not in {"wgs", "wes", "rnaseq", "targeted_mhc", "scrnaseq"}:
            raise ValueError(f"crosswalk row has invalid modality: {row['modality']!r}")
        if row["access_class"] not in ACCESS_CLASSES:
            raise ValueError(f"crosswalk row has invalid access_class: {row['access_class']!r}")
        if row["independence_stratum"] not in INDEPENDENCE_STRATA:
            raise ValueError("crosswalk row has invalid independence_stratum")
        if row["gate_status"] not in PILOT_GATES:
            raise ValueError(f"crosswalk row has invalid gate_status: {row['gate_status']!r}")
        if row["gate_status"] != "ready" and not row["exclusion_reason"].strip():
            raise ValueError(f"blocked crosswalk row lacks exclusion_reason: {key}")
        if row["gate_status"] == "excluded_related" and not row["relatedness"].strip():
            raise ValueError(f"related exclusion lacks relationship evidence: {key}")
        if row["gate_status"] == "ready":
            required = ("cohort", "repository_sample", "donor", "sequencing_accession",
                        "biosample_accession", "experiment_accession", "file_accessions",
                        "source_urls", "truth_identifier")
            missing = [field for field in required if not row[field].strip()]
            if missing:
                raise ValueError(f"ready crosswalk row missing {missing}: {key}")


def _validate_pilot(rows: list[dict[str, str]]) -> None:
    _require_header(rows, PILOT_FIELDS, "pilot manifest")
    reject_truth_columns(rows, "dataset-discovery pilot manifest")
    seen = set()
    for row in rows:
        key = (row["cohort"], row["subject"], row["modality"], row["gene"])
        if key in seen:
            raise ValueError(f"duplicate pilot locus row: {key}")
        seen.add(key)
        if row["modality"] not in {"wgs", "wes", "rnaseq"}:
            raise ValueError(f"pilot row has unsupported modality: {row['modality']!r}")
        if row["gene"] not in {"A", "B", "C"}:
            raise ValueError(f"pilot row has unsupported primary gene: {row['gene']!r}")
        if row["gate_status"] not in PILOT_GATES:
            raise ValueError(f"pilot row has invalid gate_status: {row['gate_status']!r}")
        if row["gate_status"] == "ready":
            if row["access_class"] != "public":
                raise ValueError(f"ready pilot row is not public: {key}")
            required = ("source", "source_accession", "input_type", "read1_url")
            missing = [field for field in required if not row[field].strip()]
            if missing:
                raise ValueError(f"ready pilot row missing {missing}: {key}")
            if row["input_type"] == "fastq" and not row["read2_url"].strip():
                raise ValueError(f"ready paired FASTQ row lacks read2_url: {key}")


def audit_dataset_discovery_registry(
    dataset_path: str,
    truth_path: str,
    crosswalk_path: str,
    pilot_path: str,
    output_json: str,
) -> dict:
    """Validate the discovery registries and emit fail-closed readiness gates."""
    datasets = read_tsv(dataset_path)
    truth = read_tsv(truth_path)
    crosswalk = read_tsv(crosswalk_path)
    pilot = read_tsv(pilot_path)
    _validate_datasets(datasets)
    _validate_truth(truth)
    _validate_crosswalk(crosswalk)
    _validate_pilot(pilot)

    dataset_status = Counter(row["audit_status"] for row in datasets)
    assay_status: dict[str, Counter] = defaultdict(Counter)
    for row in datasets:
        assay_status[row["assay"]][row["audit_status"]] += 1

    ready_rows = [row for row in pilot if row["gate_status"] == "ready"]
    ready_subject_modalities: dict[str, set[str]] = defaultdict(set)
    ready_loci: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in ready_rows:
        key = f"{row['cohort']}:{row['subject']}"
        ready_subject_modalities[key].add(row["modality"])
        ready_loci[(key, row["modality"])].add(row["gene"])

    complete_ready = {
        modality: sorted(key for (key, observed_modality), loci in ready_loci.items()
                         if observed_modality == modality and loci == {"A", "B", "C"})
        for modality in ("wgs", "wes", "rnaseq")
    }
    matched_public_wes_rna = sorted(
        key for key, modalities in ready_subject_modalities.items()
        if {"wes", "rnaseq"}.issubset(modalities)
        and key in complete_ready["wes"] and key in complete_ready["rnaseq"]
    )
    matched_resolved_wgs_rna = sorted({
        f"{row['cohort']}:{row['donor']}"
        for row in crosswalk
        if row["sample_link_status"] == "verified"
        and row["gate_status"] == "ready"
        and {"wgs", "rnaseq"}.issubset({
            candidate["modality"] for candidate in crosswalk
            if candidate["cohort"] == row["cohort"] and candidate["donor"] == row["donor"]
            and candidate["sample_link_status"] == "verified"
            and candidate["gate_status"] == "ready"
        })
    })

    verified_wes_capacity = [
        row["dataset_id"] for row in datasets
        if row["assay"] == "wes" and row["audit_status"] == "verified"
        and _yes(row["capacity_met"]) and row["sample_link_status"] == "verified"
    ]
    truth_by_cohort = defaultdict(list)
    for row in truth:
        truth_by_cohort[row["cohort"]].append(row)
    verified_wes_capacity = [dataset_id for dataset_id in verified_wes_capacity
                             if any(_yes(item["accuracy_eligible"])
                                    for item in truth_by_cohort[
                                        next(row["cohort"] for row in datasets
                                             if row["dataset_id"] == dataset_id)
                                    ])]

    lane_states = {}
    for dataset in datasets:
        dataset_id = dataset["dataset_id"]
        cohort = dataset["cohort"]
        modality = dataset["assay"]
        state = "inventory"
        reasons = []
        metadata_ready = (
            dataset["audit_status"] == "verified"
            and dataset["sample_link_status"] == "verified"
            and dataset["raw_data_status"] == "public_resolved"
        )
        if metadata_ready:
            state = "metadata_verified"
        else:
            reasons.append("metadata_or_sample_link_not_verified")
        eligible_truth = [row for row in truth_by_cohort[cohort]
                          if _yes(row["accuracy_eligible"])
                          and row["audit_status"] == "verified"]
        if metadata_ready and eligible_truth:
            state = "truth_verified"
        elif not eligible_truth:
            reasons.append("independent_truth_not_verified")
        executable_subjects = ([value for value in complete_ready.get(modality, [])
                                if PILOT_COHORT_ALIASES.get(value.split(":", 1)[0],
                                                            value.split(":", 1)[0]) == cohort]
                               if modality in complete_ready else [])
        if state == "truth_verified" and executable_subjects:
            state = "executable"
        elif modality in {"wgs", "wes", "rnaseq"} and not executable_subjects:
            reasons.append("no_complete_truth_blind_pilot_subjects")
        lane_states[dataset_id] = {
            "cohort": cohort,
            "assay": modality,
            "state": state,
            "state_index": LANE_STATES.index(state),
            "complete_executable_subjects": len(executable_subjects),
            "blocked_reasons": sorted(set(reasons)),
        }

    acceptance = {
        "public_rna_pilot_ready": len(complete_ready["rnaseq"]) >= PILOT_SUBJECT_TARGET,
        "public_wes_pilot_ready": len(complete_ready["wes"]) >= PILOT_SUBJECT_TARGET,
        "public_matched_wes_rna_pilot_ready": (
            len(matched_public_wes_rna) >= PILOT_SUBJECT_TARGET
        ),
        "matched_wgs_rna_candidate_accession_resolved": (
            len(matched_resolved_wgs_rna) >= PILOT_SUBJECT_TARGET
        ),
        "wes_minimum_closed_by_verified_cohort": bool(verified_wes_capacity),
        "overlap_report_complete": all(row["independence_stratum"] != "unresolved"
                                       and row["existing_champhla_overlap"]
                                       for row in crosswalk),
        "relatedness_report_complete": all(
            row["relatedness"]
            and "unresolved" not in row["relatedness"].lower()
            and "requires" not in row["relatedness"].lower()
            for row in crosswalk
        ),
        "truth_circularity_report_complete": all(row["truth_independence"]
                                                  and row["circularity"] for row in truth),
    }
    result = {
        "schema_version": "hla-ground-truth-dataset-discovery-1",
        "inventory_only": True,
        "raw_reads_downloaded": False,
        "controlled_access_requested": False,
        "dataset_rows": len(datasets),
        "dataset_status_counts": dict(sorted(dataset_status.items())),
        "assay_status_counts": {key: dict(sorted(value.items()))
                                for key, value in sorted(assay_status.items())},
        "truth_rows": len(truth),
        "crosswalk_rows": len(crosswalk),
        "pilot_rows": len(pilot),
        "pilot_subject_target": PILOT_SUBJECT_TARGET,
        "complete_ready_subjects": complete_ready,
        "matched_public_wes_rna_subjects": matched_public_wes_rna,
        "matched_resolved_wgs_rna_subjects": matched_resolved_wgs_rna,
        "verified_wes_capacity_datasets": verified_wes_capacity,
        "lane_states": lane_states,
        "lane_state_order": list(LANE_STATES),
        "acceptance": acceptance,
        "registry_valid": True,
        "all_readiness_gates_pass": all(acceptance.values()),
        "passed": all(acceptance.values()),
        "source_sha256": {
            "dataset_registry": sha256(dataset_path),
            "truth_registry": sha256(truth_path),
            "sample_crosswalk": sha256(crosswalk_path),
            "pilot_manifest": sha256(pilot_path),
        },
    }
    write_json(output_json, result)
    return result
