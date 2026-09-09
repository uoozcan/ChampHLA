from __future__ import annotations

import unittest

from champhla_confirmation.io import canonical_allele, canonical_pair
from champhla_confirmation.statistics import (
    exact_cluster_signflip,
    holm_adjust,
    simultaneous_cluster_bootstrap_ci,
    subject_deltas,
)


class IoStatisticsTests(unittest.TestCase):
    def test_two_field_normalization_and_symmetry(self):
        self.assertEqual("A*01:01", canonical_allele("HLA-A*01:01:01", "A"))
        self.assertEqual(("A*01:01", "A*02:01"), canonical_pair("A*02:01", "01:01:01", "A"))

    def test_gene_mismatch_and_one_field_fail(self):
        with self.assertRaises(ValueError):
            canonical_allele("B*07:02", "A")
        with self.assertRaises(ValueError):
            canonical_allele("A*01", "A")

    def test_exact_signflip(self):
        self.assertEqual(1.0, exact_cluster_signflip([]))
        self.assertAlmostEqual(0.03125, exact_cluster_signflip([1] * 6))
        self.assertAlmostEqual(0.015625, exact_cluster_signflip([1] * 7))

    def test_holm_is_monotone(self):
        rows = holm_adjust([{"p_value": 0.01}, {"p_value": 0.02}, {"p_value": 0.5}])
        self.assertEqual([0.03, 0.04, 0.5], [row["holm_adjusted_p"] for row in rows])

    def test_donor_cluster_overrides_library_subject(self):
        rows = [
            {"subject": "L1", "cluster_id": "D1", "a": 0, "b": 1},
            {"subject": "L2", "cluster_id": "D1", "a": 0, "b": 1},
        ]
        self.assertEqual({"D1": 2}, subject_deltas(rows, "a", "b"))

    def test_simultaneous_intervals_cover_each_comparator(self):
        rows = [{"subject": f"S{i}", "a": 0, "b": 1} for i in range(8)]
        intervals = simultaneous_cluster_bootstrap_ci(
            {"m1": rows, "m2": rows}, "a", "b", iterations=200, seed=4,
        )
        self.assertEqual({"m1", "m2"}, set(intervals))
        self.assertEqual((1.0, 1.0), intervals["m1"])


if __name__ == "__main__":
    unittest.main()
