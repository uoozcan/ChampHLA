from __future__ import annotations

import argparse

from .assays import build_assay_manifest
from .exposure import build_exposure_ledger
from .manuscript import audit_claims, write_docx_text, write_source_issue_register
from .provenance import build_originals_manifest
from .registry import render_tables, validate_registry
from .release import freeze_release
from .readiness import audit_release_readiness
from .recount import independent_recount
from .rosters import freeze_rosters
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
    args = parser.parse_args(argv)
    result = freeze_release(args.project_root, args.output)
    print(f"release files={result['file_count']} git_dirty={result['git_dirty']}")
    return 0


def audit_release_readiness_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Evaluate all plurality release gates")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    result = audit_release_readiness(args.project_root, args.config, args.output)
    print(f"release_ready={result['release_ready']} blockers={len(result['blockers'])}")
    return 0 if result["release_ready"] else 2


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
        "audit-release-readiness": audit_release_readiness_main,
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
