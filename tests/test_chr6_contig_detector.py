from __future__ import annotations

from pathlib import Path
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "bin" / "detect_chr6_contig.sh"


def run_detector(tmp_path: Path, header: str) -> subprocess.CompletedProcess[str]:
    path = tmp_path / "header.sam"
    path.write_text(header, encoding="utf-8")
    return subprocess.run(
        ["bash", str(SCRIPT), str(path)], capture_output=True, text=True, check=False,
    )


@pytest.mark.parametrize("contig", ["6", "chr6"])
def test_detector_accepts_exact_primary_contig_with_real_sam_tabs(tmp_path: Path, contig: str):
    result = run_detector(
        tmp_path,
        "@HD\tVN:1.6\n"
        f"@SQ\tSN:{contig}\tLN:170805979\tUR:/reference.fa\n"
        "@SQ\tSN:chr6_GL000250v2_alt\tLN:4672374\n",
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == f"{contig}\n"


@pytest.mark.parametrize(
    "header",
    [
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n",
        "@SQ\tSN:6\tLN:10\n@SQ\tSN:chr6\tLN:10\n",
        "@SQ\tSN:chr6_GL000250v2_alt\tLN:10\n",
    ],
)
def test_detector_rejects_missing_ambiguous_or_alt_only_headers(tmp_path: Path, header: str):
    result = run_detector(tmp_path, header)
    assert result.returncode == 2
    assert "expected exactly one" in result.stderr


def test_every_region_extract_uses_the_frozen_detector():
    modules = {
        "hlahd.nf": 1,
        "kourami.nf": 1,
        "spechla.nf": 1,
        "bam_to_fastq.nf": 1,
        "polysolver.nf": 1,
    }
    for name, expected in modules.items():
        text = (ROOT / "workflow" / "modules" / name).read_text(encoding="utf-8")
        assert text.count("bin/detect_chr6_contig.sh") == expected
        assert "SN:chr6[[:space:]]" not in text
