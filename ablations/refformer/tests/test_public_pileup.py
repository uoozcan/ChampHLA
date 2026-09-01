import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "generate_public_pileup_tensors.py"
SPEC = importlib.util.spec_from_file_location("public_pileup", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
SPEC.loader.exec_module(MODULE)


class PublicPileupTests(unittest.TestCase):
    def test_indel_parser_skips_payload(self):
        counts = MODULE.parse_bases(".,A+3TTT,-2ac*^F,$")
        self.assertEqual(counts["ref"], 4)
        self.assertEqual(counts["A"], 1)
        self.assertEqual(counts["ins"], 1)
        self.assertEqual(counts["del"], 2)

    def test_tensor_shape_and_symptom_channels(self):
        text = "chr6\t1\tA\t2\t.,\tII\nchr6\t2\tC\t0\t*\t*\n"
        tensor = MODULE.tensor_from_text(text)
        self.assertEqual(tensor.shape, (2, 44))
        self.assertEqual(tensor[0, 1], 1.0)
        self.assertEqual(tensor[0, 18], 1.0)
        self.assertEqual(tensor[1, 2], 1.0)
        self.assertEqual(tensor[1, 18], 0.0)

    def test_regions_have_expected_lengths(self):
        for build, regions in MODULE.REGIONS.items():
            for region in regions.values():
                start, end = map(int, region.split(":", 1)[1].split("-"))
                expected = 6001 if build == "GRCh38DH" else 5001
                self.assertEqual(end - start + 1, expected)


if __name__ == "__main__":
    unittest.main()
