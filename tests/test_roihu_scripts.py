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
