from __future__ import annotations

import random
import unittest

from champhla_confirmation.consensus import build_guarded_cc
from champhla_confirmation.schema import normalize_caller_row


def call(subject, caller, pair, modality="wgs", gene="A", status="callable"):
    return {"cohort": "C", "subject": subject, "modality": modality, "gene": gene,
            "caller": caller, "allele1": pair[0] if pair else "", "allele2": pair[1] if pair else "",
            "call_status": status, "caller_version": "1", "source_path": "", "source_sha256": ""}


def cc(subject, pair=None, modality="wgs", gene="A"):
    return {"cohort": "C", "subject": subject, "modality": modality, "gene": gene,
            "method": "ChampionChallenger", "allele1": pair[0] if pair else "",
            "allele2": pair[1] if pair else "", "call_status": "callable" if pair else "missing"}


class ConsensusTests(unittest.TestCase):
    def setUp(self):
        self.p1 = ("A*01:01", "A*02:01")
        self.p2 = ("A*03:01", "A*24:02")

    def _method(self, rows, method):
        return next(row for row in rows if row["method"] == method)

    def test_exact_two_thirds_is_protected(self):
        callers = [call("S", "HLA-HD", self.p1), call("S", "Kourami", self.p1),
                   call("S", "OptiType", self.p2)]
        rows = build_guarded_cc(callers, [cc("S", self.p2)])
        baseline = self._method(rows, "SimpleTwoThirdsConsensus")
        guarded = self._method(rows, "TwoThirdsGuardedCC")
        self.assertEqual("callable", baseline["call_status"])
        self.assertEqual(self.p1, (guarded["allele1"], guarded["allele2"]))
        self.assertEqual("protected_two_thirds_consensus", guarded["decision_reason"])

    def test_below_threshold_invokes_cc(self):
        callers = [call("S", "HLA-HD", self.p1), call("S", "Kourami", self.p1),
                   call("S", "OptiType", self.p2), call("S", "SpecHLA", self.p2),
                   call("S", "T1K", ("A*11:01", "A*26:01"))]
        rows = build_guarded_cc(callers, [cc("S", self.p2)])
        self.assertEqual("no_consensus", self._method(rows, "SimpleTwoThirdsConsensus")["call_status"])
        guarded = self._method(rows, "TwoThirdsGuardedCC")
        self.assertEqual(self.p2, (guarded["allele1"], guarded["allele2"]))
        self.assertEqual("champion_challenger_resolved_disagreement", guarded["decision_reason"])

    def test_missing_cc_uses_deterministic_plurality(self):
        callers = [call("S", "HLA-HD", self.p2), call("S", "Kourami", self.p1)]
        rows = build_guarded_cc(callers, [cc("S")])
        guarded = self._method(rows, "TwoThirdsGuardedCC")
        self.assertEqual(min(self.p1, self.p2), (guarded["allele1"], guarded["allele2"]))
        self.assertEqual("plurality_fallback_missing_cc", guarded["decision_reason"])

    def test_homozygous_pair_and_single_callable(self):
        homo = ("A*01:01", "A*01:01")
        rows = build_guarded_cc([call("S", "OptiType", homo)], [cc("S", self.p2)])
        self.assertEqual(homo, tuple(self._method(rows, "TwoThirdsGuardedCC")[field]
                                    for field in ("allele1", "allele2")))

    def test_caller_order_invariance(self):
        callers = [call("S", "HLA-HD", self.p1), call("S", "Kourami", self.p1),
                   call("S", "OptiType", self.p2)]
        expected = build_guarded_cc(callers, [cc("S", self.p2)])
        random.Random(7).shuffle(callers)
        self.assertEqual(expected, build_guarded_cc(callers, [cc("S", self.p2)]))

    def test_truth_column_is_rejected(self):
        row = call("S", "HLA-HD", self.p1)
        row["truth_allele1"] = "A*01:01"
        with self.assertRaises(ValueError):
            build_guarded_cc([row], [cc("S", self.p1)])

    def test_nottyped_is_explicit_missing(self):
        row = call("S", "HLA-HD", ("NOTTYPED", "NOTTYPED"))
        normalized = normalize_caller_row(row)
        self.assertEqual("missing", normalized["call_status"])
        self.assertEqual("", normalized["allele1"])

    def test_called_status_is_normalized_callable(self):
        row = call("S", "HLA-HD", self.p1)
        row["call_status"] = "called"
        self.assertEqual("callable", normalize_caller_row(row)["call_status"])


if __name__ == "__main__":
    unittest.main()
