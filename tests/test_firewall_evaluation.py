from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from champhla_confirmation.evaluation import capacity, evaluate, join_truth
from champhla_confirmation.freeze import freeze_bundle, validate_freeze
from champhla_confirmation.io import read_json, write_tsv


FIELDS = ["cohort", "subject", "modality", "gene", "method", "allele1", "allele2",
          "call_status", "decision_reason"]


def prediction(subject, modality, gene, method, pair, status="callable"):
    return {"cohort": "EXT", "subject": subject, "modality": modality, "gene": gene,
            "method": method, "allele1": pair[0] if pair else "", "allele2": pair[1] if pair else "",
            "call_status": status, "decision_reason": "fixture"}


class FirewallEvaluationTests(unittest.TestCase):
    def test_freeze_detects_mutation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "source.txt").write_text("code\n", encoding="utf-8")
            predictions = root / "pred.tsv"
            predictions.write_text("a\n1\n", encoding="utf-8")
            protocol = root / "protocol.json"
            protocol.write_text("{}\n", encoding="utf-8")
            manifest = root / "freeze.json"
            freeze_bundle(root, {"predictions": predictions, "protocol": protocol}, manifest)
            self.assertTrue(validate_freeze(manifest)["valid"])
            predictions.write_text("a\n2\n", encoding="utf-8")
            self.assertFalse(validate_freeze(manifest)["valid"])

    def test_freeze_hash_scope_excludes_truth_data(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "src").mkdir()
            (root / "src" / "code.py").write_text("VALUE = 1\n", encoding="utf-8")
            (root / "discovery").mkdir()
            (root / "discovery" / "truth.tsv").write_text("truth_allele1\nA*01:01\n")
            predictions = root / "pred.tsv"
            predictions.write_text("a\n1\n", encoding="utf-8")
            protocol = root / "protocol.json"
            protocol.write_text("{}\n", encoding="utf-8")
            manifest = root / "freeze.json"
            payload = freeze_bundle(root, {"predictions": predictions, "protocol": protocol}, manifest)
            frozen_paths = {row["path"] for row in payload["project_files"]}
            self.assertIn("src/code.py", frozen_paths)
            self.assertNotIn("discovery/truth.tsv", frozen_paths)
            self.assertEqual("pred.tsv", payload["inputs"]["predictions"]["path"])
            self.assertNotIn("\\", payload["project_root_locator"])
            self.assertTrue(validate_freeze(manifest)["valid"])

    def test_join_once_checksum_and_partial_evaluation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pred_path, truth_path = root / "pred.tsv", root / "truth.tsv"
            pair1, pair2 = ("A*01:01", "A*02:01"), ("A*03:01", "A*24:02")
            rows = []
            truths = []
            for index in range(8):
                subject = f"S{index}"
                rows.extend([
                    prediction(subject, "wes", "A", "SimpleTwoThirdsConsensus", None, "no_consensus"),
                    prediction(subject, "wes", "A", "SimplePluralityLex", pair2),
                    prediction(subject, "wes", "A", "TwoThirdsGuardedCC", pair1),
                    prediction(subject, "wes", "A", "ChampionChallenger", pair1),
                ])
                truths.append({"cohort": "EXT", "subject": subject, "modality": "wes", "gene": "A",
                               "truth_allele1": pair1[0], "truth_allele2": pair1[1], "truth_status": "resolved",
                               "truth_source": "orthogonal", "source_sha256": "abc"})
            write_tsv(pred_path, rows, FIELDS)
            write_tsv(truth_path, truths)
            protocol = root / "protocol.json"
            protocol.write_text("{}\n", encoding="utf-8")
            freeze_path = root / "freeze.json"
            freeze_bundle(root, {"predictions": pred_path, "protocol": protocol}, freeze_path)
            joined, join_manifest = root / "joined.tsv", root / "join.json"
            join_truth(str(pred_path), str(truth_path), str(freeze_path), str(joined), str(join_manifest))
            result = evaluate(str(joined), str(root / "eval"), bootstrap=1000, allow_partial=True,
                              mode="discovery")
            self.assertFalse(result["three_modality_claim_ready"])
            self.assertEqual("SimplePluralityLex", result["reference_method"])
            primary = (root / "eval" / "primary_modality_results.tsv").read_text(encoding="utf-8")
            self.assertIn("\t0\t0.0\t8\t1.0\t", primary)
            self.assertEqual(8, read_json(join_manifest)["loci"])

    def test_capacity_requires_three_modalities(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = []
            pair = ("A*01:01", "A*02:01")
            for modality in ("wgs", "wes", "rnaseq"):
                for index in range(7):
                    subject = f"{modality}{index}"
                    rows.extend([
                        prediction(subject, modality, "A", "SimpleTwoThirdsConsensus", None, "no_consensus"),
                        prediction(subject, modality, "A", "TwoThirdsGuardedCC", pair),
                    ])
            path = root / "pred.tsv"
            write_tsv(path, rows, FIELDS)
            result = capacity(str(path), str(root / "capacity.json"))
            self.assertTrue(result["passed"])


if __name__ == "__main__":
    unittest.main()
