from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from champhla_confirmation.audit import audit_wgs
from champhla_confirmation.cohorts import build_overlap_crosswalk
from champhla_confirmation.io import read_tsv, write_tsv
from champhla_confirmation.panels import PANELS


class OverlapAuditTests(unittest.TestCase):
    def test_alias_and_relative_exclusion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            development = root / "dev.tsv"
            candidates = root / "candidates.tsv"
            write_tsv(development, [{"subject": "HG1", "aliases": "COR1"}])
            write_tsv(candidates, [
                {"cohort": "X", "subject": "NEW1", "coriell_id": "COR1", "development_relative": ""},
                {"cohort": "X", "subject": "NEW2", "coriell_id": "", "development_relative": "true"},
                {"cohort": "X", "subject": "NEW3", "coriell_id": "", "development_relative": ""},
            ])
            result = build_overlap_crosswalk(str(development), str(candidates), str(root / "out.tsv"),
                                             str(root / "summary.json"))
            self.assertEqual(1, result["eligible_nonoverlap"])

    def test_wgs_audit_fails_closed_without_all_expected_records_and_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "harm.tsv"
            write_tsv(path, [{"sample": "S", "modality": "wgs", "gene": "A", "tool": "OptiType",
                              "allele1_raw": "A*01:01:01", "allele2_raw": "A*02:01:01",
                              "allele1": "A*01:01", "allele2": "A*02:01", "call_status": "callable",
                              "source_file": str(root / "missing.txt")}])
            result = audit_wgs(path, root / "audit")
            self.assertFalse(result["passed"])
            self.assertFalse(result["automated_passed"])
            self.assertFalse(result["manual_review_complete"])
            self.assertTrue(any("absent_expected_record" in row["audit_status"]
                                for row in read_tsv(root / "audit" / "wgs_audit.tsv")))

    def test_wgs_audit_can_pass_after_manual_review(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = []
            for index in range(10):
                for gene, pair in {"A": ("A*01:01", "A*02:01"),
                                   "B": ("B*07:02", "B*08:01"),
                                   "C": ("C*07:01", "C*07:02")}.items():
                    for caller in PANELS["wgs"]:
                        source = root / f"{index}_{gene}_{caller}.txt"
                        if caller == "OptiType":
                            source.write_text(f"{gene}1\t{gene}2\n{pair[0]}\t{pair[1]}\n", encoding="utf-8")
                        else:
                            source.write_text(f"HLA-{gene}\t{pair[0]}\t{pair[1]}\n", encoding="utf-8")
                        rows.append({"sample": f"S{index}", "modality": "wgs", "gene": gene,
                                     "tool": caller, "allele1_raw": pair[0], "allele2_raw": pair[1],
                                     "allele1": pair[0], "allele2": pair[1], "call_status": "callable",
                                     "source_file": str(source)})
            harmonized = root / "harm.tsv"
            write_tsv(harmonized, rows)
            first = audit_wgs(harmonized, root / "first")
            self.assertFalse(first["passed"])
            self.assertTrue(first["automated_passed"])
            self.assertFalse(first["manual_review_complete"])
            review = read_tsv(root / "first" / "wgs_manual_review.tsv")
            for row in review:
                row["review_status"] = "pass"
                row["reviewer"] = "fixture-reviewer"
            reviewed = root / "reviewed.tsv"
            write_tsv(reviewed, review)
            second = audit_wgs(harmonized, root / "second", reviewed)
            self.assertTrue(second["passed"])
            self.assertTrue(second["automated_passed"])
            self.assertTrue(second["manual_review_complete"])


if __name__ == "__main__":
    unittest.main()
