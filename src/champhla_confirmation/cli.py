from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from .audit import audit_wgs
from .cohorts import build_overlap_crosswalk
from .consensus import build_consensus, build_guarded_cc, build_mv_floored_cc
from .dataset_discovery import audit_dataset_discovery_registry
from .evaluation import capacity, evaluate, join_truth
from .external import build_hprc_release2_candidates, select_hprc_confirmation_roster
from .freeze import freeze_bundle, validate_freeze
from .hprc_truth import build_hprc_assembly_truth
from .io import read_json, read_tsv, reject_truth_columns, sha256, write_json, write_tsv
from .imgt_release import SOURCE_KINDS, load_release_index, resolve_source
from .manifests import (
    caller_release_map,
    validate_caller_reference_attestation,
    validate_comparator_manifest,
    validate_hprc_truth_protocol,
)
from .panels import (
    METHOD_GUARDED_CC,
    METHOD_PLURALITY,
    METHOD_RAW_CC,
    PANELS,
    PLURALITY_METHOD_VERSION,
    canonical_method,
)
from .power import simulate
from .public_reads import (
    build_ihwg_provenance_review_packet,
    finalize_ihwg_read_roster,
    sanitize_ihwg_registry,
    scout_ihwg_public_reads,
)
from .raw_cc import freeze_raw_cc_policy, predict_raw_cc
from .roihu import (
    assess_storage,
    audit_run_manifest,
    build_cleanup_plan,
    build_hprc_run_manifest,
    build_nci60_run_manifest,
    build_same_resource_run_manifest,
    resolve_same_resource_index_checksums,
    collect_run_outputs,
    freeze_workflow_lock,
    freeze_run_manifest,
    initialize_run_ledger,
    inventory_environment,
    transition_run_sample,
    validate_workflow_lock,
)
from .schema import (
    explode_candidate_rows,
    normalize_caller_row,
    normalize_method_row,
    normalize_raw_cc_row,
    validate_production_caller_matrix,
)
from .wgs_pilot import collect_full_cram_wgs_pilot


def audit_wgs_inputs_main() -> int:
    parser = argparse.ArgumentParser(description="Fail-closed audit of WGS calls against caller-native sources")
    parser.add_argument("--harmonized", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--manual-review")
    parser.add_argument("--environment-manifest", required=True)
    args = parser.parse_args()
    summary = audit_wgs(
        args.harmonized, args.output_dir, args.manual_review, args.environment_manifest,
    )
    print(f"WGS audit passed={summary['passed']} records={summary['locus_caller_records']}")
    return 0 if summary["passed"] else 2


def run_truth_blind_predictions_main() -> int:
    parser = argparse.ArgumentParser(description="Run frozen pair-level plurality without truth")
    parser.add_argument("--calls", required=True)
    raw_group = parser.add_mutually_exclusive_group(required=False)
    raw_group.add_argument("--raw-cc")
    raw_group.add_argument("--raw-cc-policy")
    parser.add_argument("--secondary")
    parser.add_argument("--mvfloor-policy", help="frozen truth-free MV-floor routing policy")
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--cohort-default", default="")
    parser.add_argument("--candidate-support-format", action="store_true")
    args = parser.parse_args()
    source_calls = read_tsv(args.calls)
    source_cc = read_tsv(args.raw_cc) if args.raw_cc else []
    reject_truth_columns(source_calls, "runtime caller calls")
    reject_truth_columns(source_cc, "runtime raw CC calls")
    calls = (explode_candidate_rows(source_calls, args.cohort_default) if args.candidate_support_format
             else [normalize_caller_row(row, args.cohort_default) for row in source_calls])
    validate_production_caller_matrix(calls)
    if args.raw_cc_policy:
        raw_cc = predict_raw_cc(calls, args.raw_cc_policy)
    elif args.raw_cc:
        raw_cc = [normalize_raw_cc_row(row, args.cohort_default) for row in source_cc
                  if row.get("method", METHOD_RAW_CC) == METHOD_RAW_CC]
    else:
        raw_cc = []
    predictions = build_guarded_cc(calls, raw_cc) if raw_cc else build_consensus(calls)
    if raw_cc:
        predictions.extend(raw_cc)
        if args.mvfloor_policy:
            predictions.extend(build_mv_floored_cc(calls, raw_cc, args.mvfloor_policy))
    elif args.mvfloor_policy:
        raise ValueError("--mvfloor-policy requires --raw-cc or --raw-cc-policy")
    for row in calls:
        predictions.append({
            "cohort": row["cohort"], "subject": row["subject"],
            "superpopulation": row.get("superpopulation", ""), "modality": row["modality"],
            "gene": row["gene"], "method": f"Caller:{row['caller']}",
            "allele1": row["allele1"], "allele2": row["allele2"],
            "call_status": row["call_status"], "decision_reason": "individual_intended_use_caller",
        })
    if args.secondary:
        secondary_source = read_tsv(args.secondary)
        reject_truth_columns(secondary_source, "runtime secondary method calls")
        predictions.extend(
            normalize_method_row(row, args.cohort_default)
            for row in secondary_source
            if row.get("method")
            and canonical_method(row["method"]) not in {METHOD_RAW_CC, METHOD_PLURALITY}
        )
    predictions.sort(key=lambda row: (row.get("cohort", ""), row["subject"], row["modality"],
                                      row["gene"], row["method"]))
    leading = ["cohort", "subject", "superpopulation", "modality", "gene", "method",
               "allele1", "allele2", "call_status", "decision_reason"]
    all_fields = {field for row in predictions for field in row}
    write_tsv(args.predictions, predictions, leading + sorted(all_fields - set(leading)))
    plurality_rows = [row for row in predictions if row["method"] == METHOD_PLURALITY]
    guarded_rows = [row for row in predictions if row["method"] == METHOD_GUARDED_CC]
    reasons = Counter(row.get("decision_reason", "") for row in predictions
                      if row["method"] == METHOD_GUARDED_CC)
    write_json(args.manifest, {
        "schema_version": "truth-blind-predictions-2", "truth_blind": True,
        "primary_method": METHOD_PLURALITY,
        "primary_method_version": PLURALITY_METHOD_VERSION,
        "production_caller_matrix_validated": True,
        "frozen_panels": {key: list(value) for key, value in PANELS.items()},
        "caller_partial_records": sum(row["call_status"] == "partial" for row in calls),
        "caller_missing_records": sum(row["call_status"] == "missing" for row in calls),
        "plurality_rows": len(plurality_rows),
        "plurality_no_evidence_loci": sum(row["call_status"] != "callable" for row in plurality_rows),
        "raw_cc_included": bool(raw_cc),
        "mvfloor_included": bool(raw_cc and args.mvfloor_policy),
        "mvfloor_policy_sha256": sha256(args.mvfloor_policy) if args.mvfloor_policy else "",
        "always_emit_guarded_cc": bool(guarded_rows) and all(row["call_status"] == "callable" for row in guarded_rows),
        "guarded_no_evidence_loci": sum(row["call_status"] != "callable" for row in guarded_rows),
        "homozygosity_guard_primary": False,
        "prediction_rows": len(predictions),
        "guarded_decision_counts": dict(sorted(reasons.items())),
        "input_sha256": {"calls": sha256(args.calls),
                         "raw_cc_or_policy": sha256(args.raw_cc or args.raw_cc_policy)
                         if (args.raw_cc or args.raw_cc_policy) else "not_supplied"},
        "predictions_sha256": sha256(args.predictions),
    })
    return 0


def freeze_confirmation_bundle_main() -> int:
    parser = argparse.ArgumentParser(description="Freeze code, policies, predictions, and truth-blind gates")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--input", action="append", default=[], metavar="NAME=PATH")
    args = parser.parse_args()
    inputs = {}
    for item in args.input:
        if "=" not in item:
            raise ValueError(f"--input must be NAME=PATH: {item}")
        name, path = item.split("=", 1)
        if not name or name in inputs:
            raise ValueError(f"invalid or duplicate freeze input name: {name!r}")
        inputs[name] = path
    required = {"predictions", "protocol", "comparators", "amendment", "wgs_audit"}
    if not required.issubset(inputs):
        raise ValueError(f"external freeze missing required inputs: {sorted(required - set(inputs))}")
    plurality = [row for row in read_tsv(inputs["predictions"])
                 if row.get("method") == METHOD_PLURALITY]
    if not plurality:
        raise ValueError("external freeze requires pair-level plurality prediction rows")
    amendment = read_json(inputs["amendment"])
    if amendment.get("status") != "SIGNED_BY_AUTHOR" or not amendment.get("signed_by"):
        raise ValueError("external freeze requires a signed consensus-primary amendment")
    if not read_json(inputs["wgs_audit"]).get("passed"):
        raise ValueError("external freeze requires a passed WGS audit")
    comparator_failures = validate_comparator_manifest(inputs["comparators"], require_frozen=True)
    if comparator_failures:
        raise ValueError(f"external freeze requires a frozen comparator manifest: {comparator_failures}")
    freeze_bundle(args.project_root, inputs, args.output)
    return 0


def validate_prediction_freeze_main() -> int:
    parser = argparse.ArgumentParser(description="Validate an immutable confirmation freeze")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = validate_freeze(args.manifest)
    write_json(args.output, result)
    print(f"freeze valid={result['valid']} failures={len(result['failures'])}")
    return 0 if result["valid"] else 2


def join_external_truth_main() -> int:
    parser = argparse.ArgumentParser(description="Join independently frozen truth to immutable predictions once")
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--freeze-manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--join-manifest", required=True)
    args = parser.parse_args()
    if Path(args.join_manifest).exists() or Path(args.output).exists():
        raise ValueError("truth join output already exists; join-once policy forbids overwrite")
    join_truth(args.predictions, args.truth, args.freeze_manifest, args.output, args.join_manifest)
    return 0


def evaluate_confirmation_main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate plurality against frozen comparator families")
    parser.add_argument("--joined", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--bootstrap", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20260831)
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--mode", choices=("discovery", "external"), default="discovery")
    parser.add_argument("--wgs-audit-summary")
    parser.add_argument("--evaluation-design")
    parser.add_argument("--comparator-manifest")
    args = parser.parse_args()
    result = evaluate(
        args.joined, args.output_dir, args.bootstrap, args.seed, args.allow_partial,
        args.mode, args.wgs_audit_summary, args.evaluation_design, args.comparator_manifest,
    )
    print(f"confirmation status={result['status']} "
          f"three_modality_claim_ready={result['three_modality_claim_ready']}")
    return 0


def check_confirmation_capacity_main() -> int:
    parser = argparse.ArgumentParser(description="Truth-blind mathematical capacity gate")
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = capacity(args.predictions, args.output)
    print(f"capacity passed={result['passed']} discordant_subjects={result['discordant_subjects']}")
    return 0 if result["passed"] else 2


def build_overlap_crosswalk_main() -> int:
    parser = argparse.ArgumentParser(description="Exclude external aliases and relatives overlapping development")
    parser.add_argument("--development", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()
    build_overlap_crosswalk(args.development, args.candidates, args.output, args.summary)
    return 0


def simulate_confirmation_power_main() -> int:
    parser = argparse.ArgumentParser(description="Discovery-based subject-cluster power projection")
    parser.add_argument("--paired", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--wgs", type=int, default=120)
    parser.add_argument("--wes", type=int, default=89)
    parser.add_argument("--rnaseq", type=int, default=150)
    parser.add_argument("--iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260831)
    parser.add_argument("--wgs-audit-summary")
    args = parser.parse_args()
    simulate(args.paired, args.output, {"wgs": args.wgs, "wes": args.wes, "rnaseq": args.rnaseq},
             args.iterations, args.seed, args.wgs_audit_summary)
    return 0


def build_hprc_release2_candidates_main() -> int:
    parser = argparse.ArgumentParser(description="Build a truth-free, non-overlapping HPRC R2 roster")
    parser.add_argument("--assemblies", required=True)
    parser.add_argument("--illumina", required=True)
    parser.add_argument("--sample-metadata", required=True)
    parser.add_argument("--development", required=True)
    parser.add_argument("--historical-exclusions", required=True)
    parser.add_argument("--source-commit", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()
    result = build_hprc_release2_candidates(
        args.assemblies, args.illumina, args.sample_metadata, args.development,
        args.historical_exclusions, args.output, args.summary, args.source_commit,
    )
    print(f"HPRC eligible non-overlap={result['eligible_nonoverlap']} "
          f"target120={result['target_120_available']}")
    return 0 if result["target_120_available"] else 2


def select_hprc_confirmation_roster_main() -> int:
    parser = argparse.ArgumentParser(description="Freeze a truth-free stratified HPRC WGS roster")
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--target", type=int, default=120)
    parser.add_argument("--seed", default="ChampHLA-CC-confirmation-20260831")
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()
    result = select_hprc_confirmation_roster(
        args.candidates, args.target, args.seed, args.output, args.summary,
    )
    print(f"HPRC selected={result['selected']} target={result['target']}")
    return 0


def build_frozen_raw_cc_policy_main() -> int:
    parser = argparse.ArgumentParser(description="Freeze raw CC weights from development labels only")
    parser.add_argument("--development-harmonized", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = read_tsv(args.development_harmonized)
    result = freeze_raw_cc_policy(rows, args.development_harmonized, args.output)
    print(f"raw CC policy weights={len(result['diagnostics'])} output={args.output}")
    return 0


def sanitize_ihwg_registry_main() -> int:
    parser = argparse.ArgumentParser(description="Remove genotype fields from an IHWG name registry")
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()
    result = sanitize_ihwg_registry(args.source, args.output, args.manifest)
    print(f"sanitized IHWG subjects={result['subjects']} truth_blind={result['truth_blind']}")
    return 0


def scout_ihwg_public_reads_main() -> int:
    parser = argparse.ArgumentParser(description="Truth-free ENA discovery for public IHWG WES/RNA reads")
    parser.add_argument("--registry", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--sleep-seconds", type=float, default=0.1)
    parser.add_argument("--limit", type=int, default=1000)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    result = scout_ihwg_public_reads(
        args.registry, args.output, args.summary, args.cache_dir,
        args.sleep_seconds, args.limit, args.offline,
    )
    print(f"ENA candidate subjects={result['candidate_subjects_by_modality']} "
          f"failures={len(result['query_failures'])}")
    return 0 if not result["query_failures"] else 2


def build_ihwg_provenance_review_packet_main() -> int:
    parser = argparse.ArgumentParser(description="Build a fail-closed IHWG run provenance review sheet")
    parser.add_argument("--matches", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()
    result = build_ihwg_provenance_review_packet(args.matches, args.output, args.summary)
    print(f"IHWG runs requiring provenance review={result['candidate_runs']}")
    return 0


def finalize_ihwg_read_roster_main() -> int:
    parser = argparse.ArgumentParser(description="Freeze provenance-verified IHWG WES/RNA run rosters")
    parser.add_argument("--review", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--min-wes", type=int, default=89)
    parser.add_argument("--min-rnaseq", type=int, default=130)
    args = parser.parse_args()
    result = finalize_ihwg_read_roster(
        args.review, args.output, args.summary, args.min_wes, args.min_rnaseq,
    )
    print(f"IHWG selected={result['selected_subjects_by_modality']} "
          f"minimum_met={result['minimum_met']}")
    return 0 if all(result["minimum_met"].values()) else 2


def collect_full_cram_wgs_pilot_main() -> int:
    parser = argparse.ArgumentParser(description="Collect strict truth-free full-CRAM WGS pilot calls")
    parser.add_argument("--caller-root", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()
    result = collect_full_cram_wgs_pilot(
        args.caller_root, args.manifest, args.sample, args.output, args.summary,
    )
    print(f"WGS pilot samples={result['samples']} records={result['observed_locus_caller_records']}")
    return 0


def audit_dataset_discovery_registry_main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate HLA-ground-truth dataset, truth, crosswalk, and pilot registries"
    )
    parser.add_argument("--datasets", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--crosswalk", required=True)
    parser.add_argument("--pilot", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = audit_dataset_discovery_registry(
        args.datasets, args.truth, args.crosswalk, args.pilot, args.output,
    )
    print(f"dataset discovery passed={result['passed']} acceptance={result['acceptance']}")
    return 0 if result["passed"] else 2


def audit_run_manifest_main() -> int:
    parser = argparse.ArgumentParser(description="Audit a truth-free Roihu run manifest")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = audit_run_manifest(args.manifest)
    write_json(args.output, result)
    print(f"run manifest passed={result['passed']} rows={result['rows']}")
    return 0 if result["passed"] else 2


def freeze_run_manifest_main() -> int:
    parser = argparse.ArgumentParser(description="Freeze a validated truth-free run manifest")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-counts", help="optional JSON modality-to-sample mapping")
    args = parser.parse_args()
    expected = read_json(args.expected_counts) if args.expected_counts else None
    freeze_run_manifest(args.manifest, args.output, expected)
    return 0


def initialize_run_ledger_main() -> int:
    parser = argparse.ArgumentParser(description="Create the caller-level Roihu run ledger")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--git-commit", default="")
    parser.add_argument("--job-id", default="")
    args = parser.parse_args()
    rows = initialize_run_ledger(args.manifest, args.output, args.git_commit, args.job_id)
    print(f"run ledger records={len(rows)}")
    return 0


def transition_run_ledger_main() -> int:
    parser = argparse.ArgumentParser(description="Atomically transition one Roihu sample in the ledger")
    parser.add_argument("--ledger", required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--modality", required=True)
    parser.add_argument("--state", required=True)
    parser.add_argument("--exit-code", default="")
    parser.add_argument("--runtime-seconds", default="")
    parser.add_argument("--peak-memory-bytes", default="")
    parser.add_argument("--peak-disk-bytes", default="")
    parser.add_argument("--retained-disk-bytes", default="")
    parser.add_argument("--output-sha256", default="")
    parser.add_argument("--job-id", default="")
    parser.add_argument("--attempt", default="")
    parser.add_argument("--supersedes-job-id", default="")
    args = parser.parse_args()
    updates = {
        key: value for key, value in {
            "exit_code": args.exit_code, "runtime_seconds": args.runtime_seconds,
            "peak_memory_bytes": args.peak_memory_bytes,
            "peak_disk_bytes": args.peak_disk_bytes,
            "retained_disk_bytes": args.retained_disk_bytes,
            "output_sha256": args.output_sha256,
            "job_id": args.job_id, "attempt": args.attempt,
            "supersedes_job_id": args.supersedes_job_id,
        }.items() if value != ""
    }
    transition_run_sample(
        args.ledger, args.cohort, args.sample, args.modality, args.state, **updates,
    )
    return 0


def audit_roihu_environment_main() -> int:
    parser = argparse.ArgumentParser(description="Inventory Roihu tools, storage, and frozen artifacts")
    parser.add_argument("--site-config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = inventory_environment(args.site_config, args.output)
    print(f"Roihu environment passed={result['passed']} failures={len(result['failures'])}")
    return 0 if result["passed"] else 2


def freeze_workflow_lock_main() -> int:
    parser = argparse.ArgumentParser(description="Freeze the repository-owned Roihu workflow")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--workflow-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = freeze_workflow_lock(args.project_root, args.workflow_root, args.output)
    print(f"workflow lock frozen files={len(result['files'])}")
    return 0


def resolve_caller_database_release_main() -> int:
    """Identify a caller database's IPD-IMGT/HLA release from its content.

    A filename is not provenance. This resolves a database source against a pinned cache of
    official per-release allele lists and reports NO_MATCH rather than guessing.
    """
    parser = argparse.ArgumentParser(
        description="Identify a caller database release from content")
    parser.add_argument("--component", required=True,
                        help="extracted database source to identify")
    parser.add_argument("--kind", required=True, choices=sorted(SOURCE_KINDS))
    parser.add_argument("--allele-list-cache",
                        help="pinned per-release allele lists; required unless kind is declared_header")
    parser.add_argument("--caller")
    parser.add_argument("--output")
    args = parser.parse_args()
    record = resolve_source(args.component, args.kind, args.allele_list_cache)
    if args.caller:
        record["caller"] = args.caller
    record["component"] = args.component
    if args.allele_list_cache:
        record["allele_list_releases"] = len(load_release_index(args.allele_list_cache))
    if args.output:
        write_json(args.output, record)
    print(f"{args.caller or args.component}: release={record['release']} "
          f"evidence={record.get('evidence')}")
    return 0 if record["release"] != "NO_MATCH" else 2


def validate_caller_reference_attestation_main() -> int:
    """Gate production on component-level database-release evidence."""
    parser = argparse.ArgumentParser(
        description="Validate the caller/database reference attestation")
    parser.add_argument("--attestation", required=True)
    parser.add_argument("--require-ready", action="store_true")
    parser.add_argument("--output")
    args = parser.parse_args()
    failures = validate_caller_reference_attestation(args.attestation, args.require_ready)
    result = {
        "attestation": args.attestation,
        "require_ready": args.require_ready,
        "passed": not failures,
        "failures": failures,
        "derived_releases": caller_release_map(args.attestation),
    }
    if args.output:
        write_json(args.output, result)
    for failure in failures:
        print(f"FAIL {failure}")
    print(f"caller-reference attestation passed={result['passed']}")
    return 0 if result["passed"] else 2


def validate_workflow_lock_main() -> int:
    parser = argparse.ArgumentParser(description="Validate the frozen Roihu workflow lock")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = validate_workflow_lock(args.lock, args.project_root)
    if args.output:
        write_json(args.output, result)
    print(f"workflow lock passed={result['passed']} files={result['files']}")
    return 0 if result["passed"] else 2


def assess_roihu_storage_main() -> int:
    parser = argparse.ArgumentParser(description="Apply the pilot-derived Roihu storage gate")
    parser.add_argument("--pilot-ledger", required=True)
    parser.add_argument("--targets", required=True, help="JSON mapping modality to sample count")
    parser.add_argument("--environment-inventory", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    inventory = read_json(args.environment_inventory)
    if not inventory.get("passed"):
        raise ValueError("storage gate requires a passed environment inventory")
    project_storage = inventory.get("project_storage", {})
    if not project_storage.get("filesystem_global_space_ignored"):
        raise ValueError("environment inventory does not contain project-allocation accounting")
    result = assess_storage(
        args.pilot_ledger, read_json(args.targets), project_storage.get("free_bytes", 0),
        args.output, "project_allocation", sha256(args.environment_inventory),
    )
    print(f"storage gate passed={result['passed']} required={result['required_available_bytes']}")
    return 0 if result["passed"] else 2


def plan_roihu_cleanup_main() -> int:
    parser = argparse.ArgumentParser(description="Write a dry-run cleanup plan for validated work files")
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--freeze-validation", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = build_cleanup_plan(args.run_root, args.freeze_validation, args.output)
    print(f"cleanup executable={result['executable']} paths={len(result['eligible_paths'])}")
    return 0 if result["executable"] else 2


def build_hprc_assembly_truth_main() -> int:
    parser = argparse.ArgumentParser(description="Build conservative dual-method HPRC assembly truth")
    parser.add_argument("--calls", required=True)
    parser.add_argument("--truth-output", required=True)
    parser.add_argument("--audit-output", required=True)
    args = parser.parse_args()
    result = build_hprc_assembly_truth(args.calls, args.truth_output, args.audit_output)
    print(f"HPRC truth resolved={result['resolved_loci']} unresolved={result['unresolved_loci']}")
    return 0


def collect_roihu_outputs_main() -> int:
    parser = argparse.ArgumentParser(description="Collect and reparse complete Roihu caller outputs")
    parser.add_argument("--caller-root", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()
    result = collect_run_outputs(args.caller_root, args.manifest, args.output, args.summary)
    print(f"caller collection passed={result['passed']} records={result['observed_records']}")
    return 0 if result["passed"] else 2


def build_nci60_run_manifest_main() -> int:
    parser = argparse.ArgumentParser(description="Build the strict 11-subject NCI-60 RNA manifest")
    parser.add_argument("--pilot", required=True)
    parser.add_argument("--ena-report", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = build_nci60_run_manifest(args.pilot, args.ena_report, args.output)
    print(f"NCI-60 run manifest passed={result['passed']} rows={result['rows']}")
    return 0


def build_same_resource_run_manifest_main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact truth-free 1000G rerun manifest")
    parser.add_argument("--roster", required=True)
    parser.add_argument("--assay-manifest", required=True)
    parser.add_argument("--ena-report", required=True)
    parser.add_argument("--index-checksums", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = build_same_resource_run_manifest(
        args.roster, args.assay_manifest, args.ena_report, args.index_checksums, args.output,
    )
    print(f"same-resource run manifest passed={result['passed']} rows={result['rows']}")
    return 0


def resolve_same_resource_index_checksums_main() -> int:
    parser = argparse.ArgumentParser(description="Stream and hash missing same-resource CRAI files")
    parser.add_argument("--roster", required=True)
    parser.add_argument("--assay-manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = resolve_same_resource_index_checksums(
        args.roster, args.assay_manifest, args.output, args.workers,
    )
    print(f"resolved public index checksums={result['records']}")
    return 0


def build_hprc_run_manifest_main() -> int:
    parser = argparse.ArgumentParser(description="Build the frozen 120-subject HPRC WGS manifest")
    parser.add_argument("--roster", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = build_hprc_run_manifest(args.roster, args.output)
    print(f"HPRC run manifest passed={result['passed']} rows={result['rows']}")
    return 0


def validate_hprc_truth_protocol_main() -> int:
    parser = argparse.ArgumentParser(description="Validate pinned HPRC assembly-truth execution inputs")
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--allow-draft", action="store_true")
    parser.add_argument("--output")
    args = parser.parse_args()
    failures = validate_hprc_truth_protocol(args.protocol, not args.allow_draft)
    result = {"passed": not failures, "failures": failures, "protocol_sha256": sha256(args.protocol)}
    if args.output:
        write_json(args.output, result)
    print(f"HPRC truth protocol passed={not failures} failures={len(failures)}")
    return 0 if not failures else 2
