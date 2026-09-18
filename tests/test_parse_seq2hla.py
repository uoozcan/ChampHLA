"""seq2HLA's "Confidence" column is a p-value, and the filter must read it that way.

On a real Roihu run seq2HLA typed HG00096 correctly --

    A   A*29:02  0.00603534   A*01:01  0.02059562
    B   B*08:01  2.53597e-06  B*44:03' 9.12599e-06
    C   C*07:01' 0.02720768   C*16:01  0.04395955

-- and the parser discarded every one of them, emitting only DQA1 (p 0.135, 0.217)
and DPA1 (p 0.28, 0.023). It skipped rows whose values were <= 0.1, treating the
column as a confidence score.

seq2HLA's README describes its output as "a p-value for each call", so small is
strong. The filter had it backwards: it dropped the best calls and kept the worst.
These tests use the real numbers from that run.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "bin" / "parse_seq2hla_results.py"

CLASS_I = (
    "#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence\n"
    "A\tA*29:02\t0.00603534\tA*01:01\t0.02059562\n"
    "B\tB*08:01\t2.53597e-06\tB*44:03'\t9.12599e-06\n"
    "C\tC*07:01'\t0.02720768\tC*16:01\t0.04395955\n"
)
CLASS_II = (
    "#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence\n"
    "DQA1\tDQA1*02:01\t0.1351333\tDQA1*05:01\t0.2173955\n"
    "DQB1\tDQB1*02:01'\t0.0\tDQB1*02:01\tNA\n"
    "DRB1\tDRB1*03:01\t0.0009694682\tDRB1*07:01\t0.009087313\n"
    "DRA\tDRA*01:01\tNA\tDRA*01:01\tNA\n"
)
NONCLASS = (
    "#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence\n"
    "E\tE*01:01\tNA\tE*01:01\tNA\n"
    "P\tno\tNA\tno\tNA\n"
)


def run_parser():
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        (d / "S.-ClassI-class.HLAgenotype4digits").write_text(CLASS_I, encoding="utf-8")
        (d / "S.-ClassII.HLAgenotype4digits").write_text(CLASS_II, encoding="utf-8")
        (d / "S.-ClassI-nonclass.HLAgenotype4digits").write_text(NONCLASS, encoding="utf-8")
        out = d / "out.txt"
        subprocess.run([sys.executable, str(SCRIPT), "--sample", "S",
                        "--prefix", "S.", "--output", str(out)],
                       check=True, cwd=str(d), capture_output=True)
        rows = {}
        for line in out.read_text(encoding="utf-8").splitlines():
            if line.startswith("#") or line.startswith("Gene"):
                continue
            gene, a1, a2 = line.split("\t")[:3]
            rows[gene] = (a1, a2)
        return rows


class PValueDirectionTests(unittest.TestCase):
    def setUp(self):
        self.rows = run_parser()

    def test_class_one_survives_and_is_correct(self):
        """The regression: these were being discarded for being too significant."""
        self.assertEqual(self.rows.get("A"), ("A*29:02", "A*01:01"))
        self.assertEqual(self.rows.get("B"), ("B*08:01", "B*44:03"))
        self.assertEqual(self.rows.get("C"), ("C*07:01", "C*16:01"))

    def test_a_genuinely_unsupported_locus_is_dropped(self):
        """DQA1: both p above 0.1, so neither allele is supported."""
        self.assertNotIn("DQA1", self.rows)

    def test_a_strongly_supported_class_two_locus_survives(self):
        self.assertEqual(self.rows.get("DRB1"), ("DRB1*03:01", "DRB1*07:01"))

    def test_na_is_unscored_and_does_not_exclude(self):
        """DQB1 has p=0.0 and NA; NA means unscored, not unsupported."""
        self.assertIn("DQB1", self.rows)

    def test_the_ambiguity_apostrophe_is_stripped(self):
        for gene, (a1, a2) in self.rows.items():
            self.assertNotIn("'", a1, gene)
            self.assertNotIn("'", a2, gene)

    def test_non_classical_loci_are_excluded(self):
        for gene in ("E", "P", "DRA"):
            self.assertNotIn(gene, self.rows)


if __name__ == "__main__":
    unittest.main()
