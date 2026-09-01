import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "evaluate_real_reference_ranker.py"
SPEC = importlib.util.spec_from_file_location("real_eval", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
SPEC.loader.exec_module(MODULE)


class RealEvaluationTests(unittest.TestCase):
    def test_mcnemar_identical(self):
        row = {"sample": "S1", "modality": "wgs", "gene": "A", "is_correct_2field": "1"}
        candidate = {MODULE.key(row): row}
        self.assertEqual(MODULE.mcnemar_exact(candidate, candidate), (1.0, 0, 0))

    def test_holm_is_monotone_in_sorted_order(self):
        rows = [{"p_value": 0.01}, {"p_value": 0.03}, {"p_value": 0.02}]
        adjusted = MODULE.holm(rows)
        ordered = sorted(adjusted, key=lambda row: row["p_value"])
        self.assertEqual([row["holm_p"] for row in ordered], sorted(row["holm_p"] for row in ordered))


if __name__ == "__main__":
    unittest.main()
