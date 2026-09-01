import tempfile
import unittest
from pathlib import Path

from champhla_confirmation.io import write_json
from champhla_confirmation.raw_cc import FROZEN_POLICY, predict_raw_cc


class RawCCTests(unittest.TestCase):
    def test_truth_free_runtime_and_order_invariance(self):
        with tempfile.TemporaryDirectory() as tmp:
            bundle = {
                "schema_version": "truth-free-raw-cc-policy-1", "policies": FROZEN_POLICY,
                "weights": {"wes": {caller: {gene: 1.0 for gene in "ABC"}
                                      for caller in ("HLA-HD", "OptiType", "POLYSOLVER", "SpecHLA", "T1K")}},
            }
            path = Path(tmp) / "policy.json"
            write_json(path, bundle)
            calls = [
                {"cohort": "c", "subject": "s", "modality": "wes", "gene": "A",
                 "caller": caller, "allele1": pair[0], "allele2": pair[1], "call_status": "callable"}
                for caller, pair in [
                    ("OptiType", ("A*01:01", "A*02:01")),
                    ("HLA-HD", ("A*03:01", "A*24:02")),
                    ("POLYSOLVER", ("A*03:01", "A*24:02")),
                    ("SpecHLA", ("A*03:01", "A*24:02")),
                    ("T1K", ("A*03:01", "A*24:02")),
                ]
            ]
            one = predict_raw_cc(calls, str(path))[0]
            two = predict_raw_cc(list(reversed(calls)), str(path))[0]
            self.assertEqual((one["allele1"], one["allele2"]), ("A*03:01", "A*24:02"))
            self.assertEqual(one, two)
            self.assertNotIn("truth", "".join(one))

            contaminated = [dict(calls[0], truth_allele1="A*01:01")]
            with self.assertRaises(ValueError):
                predict_raw_cc(contaminated, str(path))

    def test_missing_champion_can_fall_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            bundle = {
                "schema_version": "truth-free-raw-cc-policy-1", "policies": FROZEN_POLICY,
                "weights": {"wgs": {caller: {gene: 1.0 for gene in "ABC"}
                                      for caller in ("HLA-HD", "Kourami", "OptiType", "SpecHLA", "T1K")}},
            }
            path = Path(tmp) / "policy.json"
            write_json(path, bundle)
            calls = [{"cohort": "c", "subject": "s", "modality": "wgs", "gene": "A",
                      "caller": caller, "allele1": "A*01:01", "allele2": "A*02:01",
                      "call_status": "callable"} for caller in ("HLA-HD", "Kourami", "SpecHLA")]
            result = predict_raw_cc(calls, str(path))[0]
            self.assertEqual(result["decision_reason"], "champion_missing_weighted_consensus")
            self.assertEqual(result["call_status"], "callable")


if __name__ == "__main__":
    unittest.main()
