from __future__ import annotations

import unittest

from champhla_confirmation.parsers import parse_caller_call, parse_caller_pair


class CallerParserTests(unittest.TestCase):
    def test_optitype_matrix(self):
        text = "A1\tA2\tB1\tB2\tC1\tC2\nA*01:01\tA*02:01\tB*07:02\tB*08:01\tC*07:01\tC*07:02\n"
        self.assertEqual(("A*01:01", "A*02:01"), parse_caller_pair("OptiType", text, "A"))

    def test_optitype_dash_is_partial(self):
        text = "A1\tA2\tB1\tB2\tC1\tC2\nA*01:01\t-\tB*07:02\tB*08:01\tC*07:01\tC*07:02\n"
        self.assertEqual(
            {"allele1": "A*01:01", "allele2": "", "call_status": "partial"},
            parse_caller_call("OptiType", text, "A"),
        )

    def test_line_wrappers_for_all_wgs_callers(self):
        for caller in ("HLA-HD", "Kourami", "SpecHLA", "T1K"):
            with self.subTest(caller=caller):
                text = "HLA-A\tA*03:01:01\tA*24:02:01\n"
                self.assertEqual(("A*03:01", "A*24:02"), parse_caller_pair(caller, text, "A"))

    def test_homozygous_duplicate(self):
        self.assertEqual(("C*07:02", "C*07:02"),
                         parse_caller_pair("HLA-HD", "C\tC*07:02\tC*07:02\n", "C"))

    def test_ambiguous_many_pairs_fails(self):
        text = "A*01:01 A*02:01 A*03:01"
        self.assertIsNone(parse_caller_pair("T1K", text, "A"))

    def test_pseudogene_suffix_is_not_misparsed_as_class_i(self):
        text = ("A\tHLA-A*03:01:01\t-\n"
                "DMA\tHLA-DMA*01:01:01\tHLA-DMA*01:02:01\n"
                "DOB\tHLA-DOB*01:01:01\tHLA-DOB*01:01:03\n")
        self.assertIsNone(parse_caller_pair("HLA-HD", text, "A"))
        self.assertEqual(
            {"allele1": "A*03:01", "allele2": "", "call_status": "partial"},
            parse_caller_call("HLA-HD", text, "A"),
        )
        self.assertIsNone(parse_caller_pair("HLA-HD", text, "B"))

    def test_single_allele_is_partial_not_homozygous(self):
        self.assertEqual(
            {"allele1": "B*07:02", "allele2": "", "call_status": "partial"},
            parse_caller_call("T1K", "B\tB*07:02\t-\n", "B"),
        )

    def test_explicit_homozygous_single_allele_is_duplicated(self):
        self.assertEqual(
            ("C*07:02", "C*07:02"),
            parse_caller_pair("HLA-HD", "C C*07:02 homozygous\n", "C"),
        )


if __name__ == "__main__":
    unittest.main()
