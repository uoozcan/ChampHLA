import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_roihu_preflight_enables_nounset_after_csc_environment():
    text = (ROOT / "scripts/roihu_preflight.sh").read_text(encoding="utf-8")
    assert text.index("source /etc/profile.d/zz-csc-env.sh") < text.index("set -u")


def test_no_script_sources_the_csc_environment_under_nounset():
    """zz-csc-env.sh reads PS1, which is unbound in a batch shell.

    Sourcing it with `set -u` already active aborts the job in under a second with
    "PS1: unbound variable". Both production sbatch scripts carried this defect and it
    was never caught, because it only fires once a job is actually submitted.
    """
    offenders = []
    for path in sorted((ROOT / "scripts").iterdir()):
        if path.suffix not in {".sh", ".sbatch"}:
            continue
        text = path.read_text(encoding="utf-8")
        if "source /etc/profile.d/zz-csc-env.sh" not in text:
            continue
        source_at = text.index("source /etc/profile.d/zz-csc-env.sh")
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("set -") and "u" in stripped.split()[1].lstrip("-"):
                if text.index(line) < source_at:
                    offenders.append(f"{path.name}: `{stripped}` precedes the CSC source")
                break
    assert not offenders, offenders


def test_wgs_scripts_fail_closed():
    run = (ROOT / "scripts/roihu_run_full_cram_wgs_pilot.sbatch").read_text(encoding="utf-8")
    repair = (ROOT / "scripts/roihu_repair_kourami_full_bam.sbatch").read_text(encoding="utf-8")
    assert "set -eo pipefail" in run
    assert "Expected exactly one nonempty" in run
    assert "set -eo pipefail" in repair
    assert "Kourami produced no nonempty result" in repair


def test_roihu_environment_loads_samtools_before_creating_venv():
    text = (ROOT / "scripts/setup_roihu_test_env.sh").read_text(encoding="utf-8")
    assert text.index("module load samtools/1.21") < text.index("python3 -m venv")
    assert "EXPECTED_PYTHON=3.11.15" in text
    assert '"${ROOT}/.venv/bin/pytest" -q' in text
    assert '"${ROOT}/.venv/bin/python" -m pytest' in text


def test_imported_wgs_pilot_artifact_hashes():
    root = ROOT / "artifacts" / "wgs_pilot_review"
    manifest = json.loads((root / "import_manifest.json").read_text(encoding="utf-8"))
    for relative, expected in manifest["files"].items():
        observed = hashlib.sha256((root / relative).read_bytes()).hexdigest()
        assert observed == expected


def test_supported_python_and_cross_platform_ci_are_pinned():
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    assert 'requires-python = ">=3.10"' in pyproject
    assert "ubuntu-latest" in workflow
    assert "windows-latest" in workflow
    assert 'python-version: ["3.10", "3.11"]' in workflow
    assert "pytest -q" in workflow
    assert "python -m pytest -q" in workflow


def test_canonical_workflow_is_fail_closed_and_submission_validates_lock():
    workflow_root = ROOT / "workflow"
    production_text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted(workflow_root.rglob("*"))
        if path.is_file() and path.suffix in {".nf", ".config"}
    )
    for forbidden in ("errorStrategy 'ignore'", "? 'retry' : 'ignore'", "|| true",
                      "No results generated", "(STUB)", "stub:"):
        assert forbidden not in production_text
    submit = (ROOT / "scripts/roihu_submit_wave.sh").read_text(encoding="utf-8")
    runner = (ROOT / "scripts/roihu_run_sample.sbatch").read_text(encoding="utf-8")
    assert "validate_workflow_lock" in submit
    assert "validate_comparator_manifest" in submit
    assert "configs/roihu_workflow_lock.json" in runner
    assert "placeholder|stub|no results generated" in runner
    stage = (ROOT / "scripts" / "roihu_stage_inputs.sbatch").read_text(encoding="utf-8")
    assert ".staging.${SLURM_JOB_ID}" in stage
    assert "quarantine/staging" in stage
    assert "mv \"${stage_root}\" \"${sample_root}\"" in stage
    assert "sequential_low_storage" in submit
    assert 'afterok:${previous_job}' in submit
    assert "pilot_retained" not in submit
    retry = (ROOT / "scripts" / "roihu_retry_sample.sh").read_text(encoding="utf-8")
    assert "--hold" in retry
    assert "supersedes-job-id" in retry
    assert "scontrol release" in retry
    prune = (ROOT / "scripts" / "roihu_prune_sample_work.sh").read_text(encoding="utf-8")
    assert "CALLERS_COMPLETE" in prune
    assert "caller_output_validation.tsv" in prune
    assert "find \"${work_root}\" -depth -type f -delete" in prune


def test_wgs_staging_streams_and_never_writes_the_source_cram():
    """A 30x CRAM is ~18 GB and no caller sees it; 137 of them do not fit the allocation."""
    text = (ROOT / "scripts/roihu_stage_inputs.sbatch").read_text(encoding="utf-8")
    wgs = text[text.index('elif [[ "${modality}" == wgs ]]'):text.index('elif [[ "${input_type}" == cram')]
    # The source is read straight from its URI, never fetched to the staging root.
    assert '-X "${input_uri}"' in wgs
    assert 'fetch_one "${input_uri}"' not in wgs
    # Only the index is fetched, and it is still verified.
    assert 'fetch_one "${index_uri}"' in wgs
    assert 'verify_one "${index_checksum}"' in wgs
    # Remote sources only; a local path would silently defeat the point.
    assert '"${input_uri}" == https://*' in wgs


def test_streamed_wgs_records_that_the_source_checksum_was_not_reverified():
    """The CRAM is never on disk, so claiming verification would be false."""
    text = (ROOT / "scripts/roihu_stage_inputs.sbatch").read_text(encoding="utf-8")
    assert '"source_checksum_verified": false' in text
    assert '"source_retained": false' in text
    assert '"source_checksum_declared"' in text
    assert "not independently reverified" in text


def test_wes_staging_still_downloads_and_verifies_its_source():
    """WES callers consume the CRAM itself, so it cannot be replaced by an extract."""
    text = (ROOT / "scripts/roihu_stage_inputs.sbatch").read_text(encoding="utf-8")
    wes = text[text.index('elif [[ "${input_type}" == cram'):text.index("else\n  echo \"Unsupported input type")]
    assert 'fetch_one "${input_uri}" "${raw}"' in wes
    assert 'verify_one "${source_checksum}" "${raw}"' in wes


def test_retained_disk_is_measured_after_the_release_not_before():
    """The gate multiplies retained bytes by the cohort size, so it must be post-release."""
    text = (ROOT / "scripts/roihu_run_sample.sbatch").read_text(encoding="utf-8")
    peak = text.index("peak_disk=$(du -sb")
    prune = text.index("roihu_prune_sample_work.sh")
    retained = text.index("retained_disk=$(du -sb")
    assert peak < prune < retained


def test_prune_releases_wes_and_rna_input_but_never_the_wgs_bam():
    text = (ROOT / "scripts/roihu_prune_sample_work.sh").read_text(encoding="utf-8")
    assert '[[ "${modality}" == wes || "${modality}" == rnaseq ]]' in text
    release = text[text.index('if [[ "${modality}" == wes'):]
    assert "mate_aware_hla" not in release
    # The record of what was staged must survive the release.
    assert 'test -s "${input_root}/staged_files_sha256.txt"' in text
    assert "! -name staged_files_sha256.txt" in text
    assert "INPUT_RELEASED_AFTER_VALIDATION" in text


def test_prune_refuses_outside_the_project_input_root():
    text = (ROOT / "scripts/roihu_prune_sample_work.sh").read_text(encoding="utf-8")
    assert '/scratch/project_2008084/champhla_plurality_inputs*' in text
    assert "CHAMPHLA_INPUT_ROOT:?" in text


def test_pilot_wave_runs_in_the_mode_production_will_use():
    """Otherwise the pilot's retained figure describes a state that never persists."""
    text = (ROOT / "scripts/roihu_submit_wave.sh").read_text(encoding="utf-8")
    pilot = text[text.index("  pilot)"):text.index("  capacity)")]
    assert "execution_mode=sequential_low_storage" in pilot
