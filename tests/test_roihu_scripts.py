import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_roihu_preflight_enables_nounset_after_csc_environment():
    text = (ROOT / "scripts/roihu_preflight.sh").read_text(encoding="utf-8")
    assert text.index("source /etc/profile.d/zz-csc-env.sh") < text.index("set -u")


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
