import csv
import hashlib
import tempfile
import unittest
from pathlib import Path

from champhla_confirmation.io import read_tsv
from champhla_confirmation.panels import PANELS
from champhla_confirmation.wgs_pilot import CALLER_SUFFIX, collect_full_cram_wgs_pilot


class WgsPilotTests(unittest.TestCase):
    def test_collection_requires_five_validated_native_summaries(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            manifest = tmp / "pilot.tsv"
            manifest.write_text("sample_id\tcram_url\nS1\thttps://example.org/S1.cram\n")
            root = tmp / "callers" / "S1"
            results = root / "results"
            results.mkdir(parents=True)
            validations = []
            text = "A\tA*01:01\tA*02:01\nB\tB*07:02\tB*08:01\nC\tC*07:01\tC*07:02\n"
            for caller in PANELS["wgs"]:
                suffix = CALLER_SUFFIX[caller]
                path = results / suffix / f"S1_{suffix}.txt"
                path.parent.mkdir()
                path.write_text(text)
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                validations.append({"caller": suffix, "native_summary": str(path), "sha256": digest})
            with (root / "caller_output_validation.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(validations[0]))
                writer.writeheader()
                writer.writerows(validations)
            (root / "CALLERS_COMPLETE").touch()
            result = collect_full_cram_wgs_pilot(
                str(tmp / "callers"), str(manifest), ["S1"],
                str(tmp / "calls.tsv"), str(tmp / "summary.json"),
            )
            self.assertEqual(result["observed_locus_caller_records"], 15)
            self.assertTrue(result["truth_blind"])
            self.assertEqual(sum(row["call_status"] == "callable" for row in read_tsv(tmp / "calls.tsv")), 15)


if __name__ == "__main__":
    unittest.main()
