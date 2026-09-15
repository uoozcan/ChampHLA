from __future__ import annotations

import argparse

from .assays import build_assay_manifest
from .exposure import build_exposure_ledger
from .figures import validate_figure_manifest
from .manuscript import audit_claims, write_docx_text, write_source_issue_register
from .provenance import build_originals_manifest
from .registry import render_tables, validate_registry
from .release import (
    build_compact_export_manifest,
    create_release_archive,
    freeze_release,
    verify_compact_export_manifest,
    verify_release_archive,
)
from .readiness import audit_goal_completion, audit_release_readiness
from .recount import independent_recount
from .rosters import freeze_rosters
from .submission import build_submission_package
from .truth import prepare_locked_truth


def build_subject_exposure_ledger_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Build a cross-project subject exposure ledger")
    parser.add_argument("--source", action="append", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args(argv)
    build_exposure_ledger(args.source, args.output, args.summary)
    return 0


def prepare_locked_1000g_truth_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Normalize and seal 2014 laboratory HLA truth")
    parser.add_argument("--source", required=True)
    parser.add_argument("--truth-output", required=True)
    parser.add_argument("--registry-output", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args(argv)
    prepare_locked_truth(args.source, args.truth_output, args.registry_output, args.manifest)
    return 0


def freeze_extension_rosters_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Freeze truth-blind 1000G extension rosters")
    parser.add_argument("--truth-registry", required=True)
    parser.add_argument("--assay-manifest", required=True)
    parser.add_argument("--exposure-ledger", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args(argv)
    result = freeze_rosters(args.truth_registry, args.assay_manifest, args.exposure_ledger,
                            args.config, args.output, args.summary)
    print(f"all_minimums_met={result['all_minimums_met']}")
    return 0 if result["all_minimums_met"] else 2


def freeze_release_bundle_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Checksum the compact release candidate")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--policy")
    args = parser.parse_args(argv)
    result = freeze_release(args.project_root, args.output, args.policy)
    print(
        f"release files={result['file_count']} git_dirty={result['git_dirty']} "
        f"valid={result['valid']} violations={len(result['violations'])}"
    )
    return 0 if result["valid"] else 2


def build_compact_export_manifest_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Build a path-aware compact Roihu export manifest")
    parser.add_argument("--root", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--list-output", required=True)
    args = parser.parse_args(argv)
    result = build_compact_export_manifest(
        args.root, args.run_id, args.policy, args.output, args.list_output,
    )
    print(f"compact files={result['file_count']} valid={result['valid']}")
    return 0 if result["valid"] else 2


def verify_compact_export_manifest_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Verify compact export files against their manifest")
    parser.add_argument("--root", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args(argv)
    failures = verify_compact_export_manifest(args.root, args.manifest)
    for failure in failures:
        print(failure)
    return 0 if not failures else 2


def create_release_archive_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Create a deterministic archive from a valid release freeze")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--freeze-manifest", required=True)
    parser.add_argument("--archive", required=True)
    args = parser.parse_args(argv)
    result = create_release_archive(args.project_root, args.freeze_manifest, args.archive)
    print(f"release archive files={result['file_count']} sha256={result['archive_sha256']}")
    return 0


def verify_release_archive_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Extract and verify a ChampHLA release archive")
    parser.add_argument("--archive", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    result = verify_release_archive(args.archive, args.output)
    print(f"archive passed={result['passed']} failures={len(result['failures'])}")
    return 0 if result["passed"] else 2


def audit_release_readiness_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Evaluate all plurality release gates")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    result = audit_release_readiness(args.project_root, args.config, args.output)
    print(f"release_ready={result['release_ready']} blockers={len(result['blockers'])}")
    return 0 if result["release_ready"] else 2


def audit_goal_completion_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Evaluate the full ChampHLA completion contract")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    result = audit_goal_completion(args.project_root, args.config, args.output)
    print(f"completion_ready={result['completion_ready']} blockers={len(result['blockers'])}")
    return 0 if result["completion_ready"] else 2


def validate_figure_manifest_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Validate reproducible ChampHLA figure artifacts")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    result = validate_figure_manifest(args.project_root, args.manifest, args.output)
    print(
        f"figures passed={result['passed']} generated={result['generated_figures']} "
        f"blocked={len(result['blocked_figures'])} failures={len(result['failures'])}"
    )
    return 0 if result["passed"] else 2


def build_submission_package_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Build deterministic journal-neutral Markdown/DOCX files")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--main", required=True)
    parser.add_argument("--supplement", required=True)
    parser.add_argument("--bibliography", required=True)
    parser.add_argument("--declarations", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--parity-audit", required=True)
    parser.add_argument("--bibliography-audit", required=True)
    parser.add_argument("--declarations-audit", required=True)
    args = parser.parse_args(argv)
    result = build_submission_package(
        args.project_root, args.main, args.supplement, args.bibliography, args.declarations,
        args.output_dir, args.manifest, args.parity_audit, args.bibliography_audit,
        args.declarations_audit,
    )
    print(
        f"submission parity={result['parity']['passed']} "
        f"bibliography={result['bibliography']['passed']} "
        f"declarations={result['declarations']['passed']}"
    )
    return 0 if (result["parity"]["passed"] and result["bibliography"]["passed"]
                 and result["declarations"]["passed"]) else 2


def independent_recount_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Independently recount exact plurality results")
    parser.add_argument("--joined", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--method", default="SimplePluralityLex")
    args = parser.parse_args(argv)
    result = independent_recount(args.joined, args.output, args.summary, args.method)
    print(f"independent recount strata={len(result['rows'])}")
    return 0


def render_manuscript_tables_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Render manuscript tables from the result registry")
    parser.add_argument("--registry", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--root")
    args = parser.parse_args(argv)
    failures = validate_registry(args.registry, args.root)
    if failures:
        raise ValueError("; ".join(failures))
    render_tables(args.registry, args.output)
    return 0


def audit_manuscript_claims_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Audit manuscript claims against frozen results")
    parser.add_argument("--manuscript", required=True)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--claims", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--root")
    parser.add_argument("--supplement")
    args = parser.parse_args(argv)
    result = audit_claims(
        args.manuscript, args.registry, args.claims, args.output, args.root, args.supplement,
    )
    print(f"submission_ready={result['submission_ready']} failures={len(result['failures'])} warnings={len(result['warnings'])}")
    return 0 if result["submission_ready"] else 2


def build_originals_manifest_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Verify immutable source copies")
    parser.add_argument("--source-map", required=True)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    result = build_originals_manifest(args.source_map, args.project_root, args.output)
    print(f"originals passed={result['passed']} records={len(result['records'])}")
    return 0 if result["passed"] else 2


def extract_source_manuscript_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Extract and audit the immutable source DOCX")
    parser.add_argument("--source", required=True)
    parser.add_argument("--text-output", required=True)
    parser.add_argument("--issues-output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args(argv)
    write_docx_text(args.source, args.text_output)
    write_source_issue_register(args.source, args.issues_output, args.summary)
    return 0


def build_1000g_assay_manifest_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Build truth-blind official WGS/WES/RNA input manifest")
    parser.add_argument("--wgs-index", required=True)
    parser.add_argument("--wes-index", required=True)
    parser.add_argument("--rna-sdrf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args(argv)
    build_assay_manifest(args.wgs_index, args.wes_index, args.rna_sdrf, args.output, args.manifest)
    return 0


def _main() -> int:
    parser = argparse.ArgumentParser(description="ChampHLA publication recovery utilities")
    subparsers = parser.add_subparsers(dest="command", required=True)
    commands = {
        "build-subject-exposure-ledger": build_subject_exposure_ledger_main,
        "prepare-locked-1000g-truth": prepare_locked_1000g_truth_main,
        "freeze-extension-rosters": freeze_extension_rosters_main,
        "freeze-release-bundle": freeze_release_bundle_main,
        "build-compact-export-manifest": build_compact_export_manifest_main,
        "verify-compact-export-manifest": verify_compact_export_manifest_main,
        "create-release-archive": create_release_archive_main,
        "verify-release-archive": verify_release_archive_main,
        "audit-release-readiness": audit_release_readiness_main,
        "audit-goal-completion": audit_goal_completion_main,
        "validate-figure-manifest": validate_figure_manifest_main,
        "build-submission-package": build_submission_package_main,
        "independent-recount": independent_recount_main,
        "render-manuscript-tables": render_manuscript_tables_main,
        "audit-manuscript-claims": audit_manuscript_claims_main,
        "build-originals-manifest": build_originals_manifest_main,
        "extract-source-manuscript": extract_source_manuscript_main,
        "build-1000g-assay-manifest": build_1000g_assay_manifest_main,
    }
    for name in commands:
        subparsers.add_parser(name, add_help=False)
    args, remaining = parser.parse_known_args()
    return commands[args.command](remaining)


if __name__ == "__main__":
    raise SystemExit(_main())
