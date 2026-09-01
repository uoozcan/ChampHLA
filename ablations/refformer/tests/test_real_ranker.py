import unittest

import torch

from champhla_refformer.real_ranker import bin_pileup


class RealRankerTests(unittest.TestCase):
    def test_bin_pileup_shape_and_extrema(self):
        pileup = torch.zeros(100, 44)
        pileup[37, 5] = 1.0
        binned = bin_pileup(pileup, 10)
        self.assertEqual(tuple(binned.shape), (10, 88))
        self.assertEqual(float(binned[:, 44 + 5].max()), 1.0)

    def test_bin_pileup_rejects_wrong_channels(self):
        with self.assertRaises(ValueError):
            bin_pileup(torch.zeros(10, 43), 5)


if __name__ == "__main__":
    unittest.main()
