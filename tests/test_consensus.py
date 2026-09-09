from __future__ import annotations

import random
import unittest

from champhla_confirmation.consensus import build_consensus, build_guarded_cc, build_mv_floored_cc
from champhla_confirmation.io import write_json
from champhla_confirmation.schema import (
    normalize_caller_row,
    normalize_method_row,
    validate_production_caller_matrix,
)
from champhla_confirmation.panels import GENES, PANELS


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

    def test_partial_call_is_retained_but_does_not_vote(self):
        partial = call("S", "HLA-HD", ("A*01:01", "-"))
        normalized = normalize_caller_row(partial)
        self.assertEqual("partial", normalized["call_status"])
        self.assertEqual("A*01:01", normalized["allele1"])
        rows = build_guarded_cc(
            [normalized, call("S", "Kourami", self.p2)], [cc("S", self.p1)]
        )
        plurality = self._method(rows, "SimplePluralityLex")
        self.assertEqual(self.p2, (plurality["allele1"], plurality["allele2"]))
        self.assertEqual(1, plurality["partial_tools"])
        self.assertEqual(1, plurality["complete_tools"])

    def test_tie_audit_and_duplicate_caller_rejection(self):
        callers = [call("S", "HLA-HD", self.p2), call("S", "Kourami", self.p1)]
        plurality = self._method(build_guarded_cc(callers, [cc("S", self.p1)]), "SimplePluralityLex")
        self.assertEqual(1, plurality["tie_at_top"])
        self.assertIn("|", plurality["tied_pairs"])
        with self.assertRaisesRegex(ValueError, "duplicate caller row"):
            build_guarded_cc(callers + [call("S", "HLA-HD", self.p2)], [cc("S", self.p1)])

    def test_complete_missingness_returns_no_evidence(self):
        rows = build_consensus([
            call("S", "HLA-HD", None, status="missing"),
            call("S", "Kourami", None, status="missing"),
        ])
        plurality = self._method(rows, "SimplePluralityLex")
        self.assertEqual("no_evidence", plurality["call_status"])
        self.assertEqual(0, plurality["complete_tools"])
        self.assertEqual(5, plurality["missing_tools"])

    def test_majority_vote_is_a_legacy_input_alias(self):
        row = normalize_method_row({
            "cohort": "C", "subject": "S", "modality": "wgs", "gene": "A",
            "method": "MajorityVote", "allele1": self.p1[0], "allele2": self.p1[1],
            "call_status": "callable",
        })
        self.assertEqual("SimplePluralityLex", row["method"])

    def test_explicit_native_homozygous_status_can_duplicate_one_allele(self):
        row = call("S", "HLA-HD", ("A*01:01", "-"), status="homozygous")
        normalized = normalize_caller_row(row)
        self.assertEqual("callable", normalized["call_status"])
        self.assertEqual(("A*01:01", "A*01:01"),
                         (normalized["allele1"], normalized["allele2"]))

    def test_production_matrix_requires_explicit_missing_rows(self):
        pairs = {
            "A": ("A*01:01", "A*02:01"),
            "B": ("B*07:02", "B*08:01"),
            "C": ("C*07:01", "C*07:02"),
        }
        rows = [
            normalize_caller_row(call("S", caller, pairs[gene], gene=gene))
            for gene in GENES for caller in PANELS["wgs"]
        ]
        validate_production_caller_matrix(rows)
        with self.assertRaisesRegex(ValueError, "incomplete production caller matrix"):
            validate_production_caller_matrix(rows[:-1])

    def test_mv_floor_applies_frozen_truth_free_routes(self):
        import tempfile
        from pathlib import Path

        callers = [call("S", "HLA-HD", self.p1), call("S", "Kourami", self.p1),
                   call("S", "OptiType", self.p2), call("S", "SpecHLA", self.p2),
                   call("S", "T1K", ("A*11:01", "A*26:01"))]
        policy = {"default_route": "mv", "policy": {
            "wgs": {"unanimous": "mv", "clear_majority": "cc",
                    "split": "cc", "single_tool": "mv"}}}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "policy.json"
            write_json(path, policy)
            result = build_mv_floored_cc(callers, [cc("S", self.p2)], path)[0]
        self.assertEqual(self.p2, (result["allele1"], result["allele2"]))
        self.assertEqual("split", result["vote_stratum"])
        self.assertEqual("cc", result["selected_route"])
        self.assertNotIn("truth", "".join(result))


if __name__ == "__main__":
    unittest.main()
