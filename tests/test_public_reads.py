import csv
import tempfile
import unittest
from pathlib import Path

from champhla_confirmation.io import read_json, read_tsv
from champhla_confirmation.public_reads import (
    _classify,
    _parse_ena_tsv,
    finalize_ihwg_read_roster,
    sanitize_ihwg_registry,
)


class PublicReadTests(unittest.TestCase):
    def test_registry_sanitizer_discards_genotype_derived_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.tsv"
            output = Path(tmp) / "clean.tsv"
            manifest = Path(tmp) / "manifest.json"
            with source.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle, delimiter="\t",
                    fieldnames=["ihw_number", "primary_name", "ancestry", "homozygous"],
                )
                writer.writeheader()
                writer.writerow({"ihw_number": "IHW1", "primary_name": "JY",
                                 "ancestry": "A", "homozygous": "Y"})
            result = sanitize_ihwg_registry(str(source), str(output), str(manifest))
            self.assertEqual(set(read_tsv(output)[0]), {"ihw_number", "primary_name", "ancestry"})
            self.assertEqual(read_json(manifest)["discarded_columns"], ["homozygous"])
            self.assertTrue(result["truth_blind"])

    def test_ena_parser_and_classification_fail_closed(self):
        text = ("run_accession\tsample_title\tscientific_name\tlibrary_strategy\t"
                "library_source\tlibrary_layout\tinstrument_platform\n"
                "SRR1\tJY\tHomo sapiens\tRNA-Seq\tTRANSCRIPTOMIC\tPAIRED\tILLUMINA\n"
                "SRR2\tJY\triver metagenome\tWGS\tMETAGENOMIC\tSINGLE\tILLUMINA\n")
        rows = _parse_ena_tsv(text)
        self.assertEqual(_classify(rows[0], "JY")["provenance_status"],
                         "requires_study_level_review")
        self.assertEqual(_classify(rows[1], "JY")["provenance_status"], "ineligible")

    def test_final_roster_requires_review_and_uses_largest_library(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            review = tmp / "review.tsv"
            rows = []
            for run, size, verified in (("SRR1", "100", "yes"), ("SRR2", "200", "yes"),
                                        ("SRR3", "300", "")):
                rows.append({
                    "ihw_number": "IHW1", "run_accession": run,
                    "study_accession": "STUDY", "modality": "rnaseq", "base_count": size,
                    "technical_match": "1", "exact_sample_title": "1",
                    "cell_line_provenance_verified": verified,
                    "orthogonal_dna_truth_available": verified,
                    "development_or_relative_overlap": "0", "excluded_study": "0",
                    "evidence_url": "https://example.org/evidence" if verified else "",
                    "reviewer": "reviewer" if verified else "", "review_notes": "",
                })
            with review.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            result = finalize_ihwg_read_roster(
                str(review), str(tmp / "roster.tsv"), str(tmp / "summary.json"), 1, 1,
            )
            selected = read_tsv(tmp / "roster.tsv")
            self.assertEqual(selected[0]["run_accession"], "SRR2")
            self.assertFalse(result["minimum_met"]["wes"])
            self.assertTrue(result["minimum_met"]["rnaseq"])


if __name__ == "__main__":
    unittest.main()
