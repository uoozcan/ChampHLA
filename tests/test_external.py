import csv
import json
import tempfile
import unittest
from pathlib import Path

from champhla_confirmation.external import (
    build_hprc_release2_candidates,
    select_hprc_confirmation_roster,
)


def write_csv(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class ExternalRosterTests(unittest.TestCase):
    def test_hprc_requires_both_haplotypes_and_excludes_relatives(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            assemblies = []
            for sample in ("KEEP", "DEV", "REL", "ONEHAP"):
                for hap in (("1",) if sample == "ONEHAP" else ("1", "2")):
                    assemblies.append({"sample_id": sample, "haplotype": hap,
                                       "assembly": f"s3://assembly/{sample}/{hap}"})
            illumina = [{"sample_id": sample, "filetype": "cram", "library_strategy": "WGS",
                         "library_layout": "PAIRED", "coverage": "30", "read_length": "150",
                         "path": f"s3://reads/{sample}.cram"}
                        for sample in ("KEEP", "DEV", "REL", "ONEHAP")]
            metadata = [
                {"sample_id": "KEEP", "alternative_id": "ALTKEEP", "biosample_id": "BS1",
                 "family_id": "F1", "paternal_id": "", "maternal_id": "", "siblings": "",
                 "population_abbreviation": "P1"},
                {"sample_id": "DEV", "alternative_id": "", "biosample_id": "BS2",
                 "family_id": "F2", "paternal_id": "", "maternal_id": "", "siblings": "",
                 "population_abbreviation": "P2"},
                {"sample_id": "REL", "alternative_id": "", "biosample_id": "BS3",
                 "family_id": "F2", "paternal_id": "DEV", "maternal_id": "", "siblings": "",
                 "population_abbreviation": "P2"},
                {"sample_id": "ONEHAP", "alternative_id": "", "biosample_id": "BS4",
                 "family_id": "F4", "paternal_id": "", "maternal_id": "", "siblings": "",
                 "population_abbreviation": "P4"},
            ]
            write_csv(tmp / "assemblies.csv", assemblies)
            write_csv(tmp / "illumina.csv", illumina)
            write_csv(tmp / "metadata.csv", metadata)
            (tmp / "development.tsv").write_text("subject\nDEV\n")
            (tmp / "excluded.txt").write_text("# none\n")
            result = build_hprc_release2_candidates(
                str(tmp / "assemblies.csv"), str(tmp / "illumina.csv"),
                str(tmp / "metadata.csv"), str(tmp / "development.tsv"),
                str(tmp / "excluded.txt"), str(tmp / "out.tsv"), str(tmp / "summary.json"), "abc",
            )
            self.assertEqual(result["eligible_nonoverlap"], 1)
            with (tmp / "out.tsv").open() as handle:
                rows = {row["subject"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["KEEP"]["eligible_nonoverlap"], "1")
            self.assertIn("known_development_relative", rows["REL"]["exclusion_reason"])
            self.assertIn("missing_phased_assembly", rows["ONEHAP"]["exclusion_reason"])
            self.assertTrue(json.loads((tmp / "summary.json").read_text())["truth_blind"])

    def test_selection_is_deterministic_and_stratified(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fields = ["subject", "population", "eligible_nonoverlap"]
            with (tmp / "candidates.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(
                    [{"subject": f"A{i}", "population": "A", "eligible_nonoverlap": "1"}
                     for i in range(6)]
                    + [{"subject": f"B{i}", "population": "B", "eligible_nonoverlap": "1"}
                       for i in range(4)]
                )
            first = select_hprc_confirmation_roster(
                str(tmp / "candidates.tsv"), 5, "seed", str(tmp / "one.tsv"), str(tmp / "one.json")
            )
            second = select_hprc_confirmation_roster(
                str(tmp / "candidates.tsv"), 5, "seed", str(tmp / "two.tsv"), str(tmp / "two.json")
            )
            self.assertEqual((tmp / "one.tsv").read_bytes(), (tmp / "two.tsv").read_bytes())
            self.assertEqual(first["population_counts"], {"A": 3, "B": 2})
            self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
