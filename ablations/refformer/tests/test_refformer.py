import copy
import csv
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from champhla_refformer.alleles import canonical_pair, normalize_allele
from champhla_refformer.data import build_loci, exclude_subjects, fit_runtime_context
from champhla_refformer.evaluation import assign_grouped_folds, fixed_universe, nested_run
from champhla_refformer.fallback import PileupEvidenceEncoder, StandaloneReferenceRanker
from champhla_refformer.model import RefFormerModel, ReferenceEmbeddingStore
from champhla_refformer.reference import ReferenceIndex, build_index, parse_nucleotide_alignment
from champhla_refformer.synthetic import (SyntheticPileupConfig, sample_candidate_pairs,
                                           sequence_pair_to_pileup, stable_pair_split)
from champhla_refformer.tokenizer import KmerTokenizer
from champhla_refformer.training import (coverage_thresholds, fit_temperature, load_bundle,
                                         predict_loci, save_bundle, train_model)


def call(sample, tool, pair, truth=("A*01:01", "A*02:01"), modality="wgs", score="0.8",
         callable_value="1", population="EUR"):
    return {"sample": sample, "superpopulation": population, "modality": modality, "gene": "A",
            "tool": tool, "allele1": pair[0], "allele2": pair[1], "truth_allele1": truth[0],
            "truth_allele2": truth[1], "is_callable": callable_value,
            "is_correct_2field": str(int(canonical_pair(*pair, "A") == canonical_pair(*truth, "A"))),
            "confidence_score": score}


def rows_fixture():
    tools = ["HLA-HD", "Kourami", "OptiType", "SpecHLA", "T1K"]
    rows = []
    for index in range(8):
        sample = f"S{index}"
        truth = ("A*01:01", "A*02:01") if index % 2 == 0 else ("A*03:01", "A*24:02")
        wrong = ("A*11:01", "A*26:01")
        for tool_index, tool in enumerate(tools):
            pair = truth if tool_index < 2 + (index % 3) else wrong
            rows.append(call(sample, tool, pair, truth, score=str(0.9 - tool_index * 0.1)))
    return rows


def embedding_store(tmp: Path):
    alleles = ["A*01:01", "A*02:01", "A*03:01", "A*24:02", "A*11:01", "A*26:01"]
    generator = torch.Generator().manual_seed(7)
    payload = {"embeddings": {allele: torch.randn((2, 16), generator=generator) for allele in alleles},
               "reference_sha256": "reference-test", "encoder_sha256": "encoder-test",
               "embedding_dim": 16, "metadata": {}}
    path = tmp / "embeddings.pt"; torch.save(payload, path)
    return ReferenceEmbeddingStore(path)


class StubReference:
    checksum = "reference-test"
    def resolve(self, allele):
        return normalize_allele(allele), {"alias_resolved": False, "deleted": False,
            "unresolved": False, "partial": False}


class AlleleAndReferenceTests(unittest.TestCase):
    def test_normalization_and_unordered_pair(self):
        self.assertEqual(normalize_allele("HLA-A*02:01:01G"), "A*02:01")
        self.assertEqual(normalize_allele("A*02010101"), "A*02:01")
        self.assertEqual(canonical_pair("A*02:01", "A*01:01"), ("A*01:01", "A*02:01"))

    def test_index_retains_multiple_sequences_and_reference_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source"
            (source / "fasta").mkdir(parents=True); (source / "alignments").mkdir()
            for gene in ("A", "B", "C"):
                (source / "fasta" / f"{gene}_nuc.fasta").write_text(
                    f">HLA:X1 {gene}*01:01:01:01 8 bp\nACGTACGT\n"
                    f">HLA:X2 {gene}*01:01:02:01 8 bp\nACGTTCGT\n")
                (source / "alignments" / f"{gene}_nuc.txt").write_text("alignment\n")
            (source / "Allelelist_history.txt").write_text(
                "HLA_ID,3590,3580\nX1,A*01:01:01:01,A*01020101\n")
            (source / "Deleted_alleles.txt").write_text(
                "AlleleID,Allele,Description\nX3,A*02:01:08,Sequence extended and renamed A*02:1040\n")
            for name in ("release_version.txt", "LICENCE.md", "Allelelist.txt"):
                (source / name).write_text("fixture\n")
            index_path = Path(tmp) / "index.json"
            payload = build_index(source, index_path)
            self.assertEqual(len(payload["groups"]["A*01:01"]["sequences"]), 2)
            reference = ReferenceIndex(index_path)
            resolved, flags = reference.resolve("A*01020101")
            self.assertTrue(flags["alias_resolved"])
            self.assertEqual(resolved, "A*01:01")
            _, deleted = reference.resolve("A*02:01:08")
            self.assertTrue(deleted["deleted"])
            _, missing = reference.resolve("A*99:99")
            self.assertTrue(missing["unresolved"])
            tampered = json.loads(index_path.read_text())
            tampered["release"] = "moving-latest"
            index_path.write_text(json.dumps(tampered))
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                ReferenceIndex(index_path)

    def test_actual_pinned_index(self):
        path = ROOT / "resources" / "imgthla_3.59.0" / "imgt_index.json"
        reference = ReferenceIndex(path)
        self.assertEqual(reference.payload["release"], "IPD-IMGT/HLA 3.59.0")
        self.assertEqual(reference.payload["git_commit"], "f03bbe951d3106ff00f1801407a7a13278c8866b")
        self.assertGreater(len(reference.payload["groups"]), 10000)

    def test_tokenizer_is_deterministic(self):
        tokenizer = KmerTokenizer(k=3, max_tokens=12)
        self.assertEqual(tokenizer.encode("ACGTAC"), tokenizer.encode("ACGTAC"))
        tokens, mask = tokenizer.encode("ACNT")
        self.assertEqual(len(tokens), 12); self.assertEqual(len(mask), 12)

    def test_imgt_alignment_dash_expansion_and_masking(self):
        with tempfile.TemporaryDirectory() as tmp:
            alignment = Path(tmp) / "A_nuc.txt"
            alignment.write_text(
                " cDNA 1\n"
                " A*01:01:01:01  AC.G|T\n"
                " A*02:01:01:01  --T-*\n"
                " A*05:01:01:01  -----\n"
                " A*06:01:01:01  -----\n"
                " A*07:01:01:01  -----\n"
                " A*08:01:01:01  -----\n"
                " A*09:01:01:01  -----\n"
                " A*10:01:01:01  -----\n"
                " A*11:01:01:01  -----\n"
                " A*12:01:01:01  -----\n"
                " cDNA 5\n"
                " A*01:01:01:01  AAC\n"
                " A*02:01:01:01  -G-\n"
                " A*05:01:01:01  ---\n"
                " A*06:01:01:01  ---\n"
                " A*07:01:01:01  ---\n"
                " A*08:01:01:01  ---\n"
                " A*09:01:01:01  ---\n"
                " A*10:01:01:01  ---\n"
                " A*11:01:01:01  ---\n"
                " A*12:01:01:01  ---\n"
                " cDNA 8\n"
                " A*03:01:01:01  TT\n"
                " A*04:01:01:01  AAAAA\n"
            )
            records, metadata = parse_nucleotide_alignment(alignment, "A")
            self.assertEqual(records["A*01:01:01:01"]["aligned_sequence"], "AC-GTAAC")
            self.assertEqual(records["A*02:01:01:01"]["aligned_sequence"], "ACTGNAGC")
            self.assertEqual(records["A*02:01:01:01"]["alignment_mask"], "11110111")
            self.assertEqual(metadata["shared_columns"], 8)
            self.assertEqual(metadata["ignored_allele_specific_blocks"], 1)


class CandidateAndLeakageTests(unittest.TestCase):
    def test_truth_poison_does_not_change_runtime_features(self):
        rows = rows_fixture()
        context = fit_runtime_context(rows)
        reference = StubReference()
        original = build_loci(rows, context, reference, include_labels=False)
        poisoned = copy.deepcopy(rows)
        for row in poisoned:
            row["truth_allele1"], row["truth_allele2"] = "A*98:98", "A*99:99"
            row["is_correct_2field"] = "0"
        runtime = build_loci(poisoned, context, reference, include_labels=False)
        self.assertEqual(original, runtime)

    def test_missing_callers_tie_homozygous_and_oracle(self):
        rows = [call("S", "HLA-HD", ("A*01:01", "A*01:01"), ("A*01:01", "A*01:01")),
                call("S", "Kourami", ("A*02:01", "A*03:01"), ("A*01:01", "A*01:01"))]
        context = fit_runtime_context(rows)
        locus = build_loci(rows, context, StubReference(), include_labels=True)[0]
        self.assertEqual(locus["candidate_count"], 2)
        self.assertEqual(locus["candidate_set_oracle"], 1)
        hom = next(item for item in locus["candidates"] if item["pair"][0] == item["pair"][1])
        self.assertEqual(hom["global_features"][6], 1.0)
        self.assertEqual(sum(token["callable"] for token in hom["tool_tokens"]), 2)

    def test_subjects_are_globally_grouped(self):
        rows = rows_fixture()
        for row in rows[:5]:
            duplicate = dict(row); duplicate["modality"] = "wes"; rows.append(duplicate)
        folds = assign_grouped_folds(rows, 3, 9)
        self.assertEqual(len(folds), 8)
        self.assertTrue(all(isinstance(value, int) for value in folds.values()))

    def test_reserved_subjects_removed_across_modalities(self):
        rows = rows_fixture()
        extra = dict(rows[0]); extra["modality"] = "wes"; rows.append(extra)
        retained, overlap = exclude_subjects(rows, {"S0", "HPRC_NOT_PRESENT"})
        self.assertEqual(overlap, {"S0"})
        self.assertNotIn("S0", {row["sample"] for row in retained})

    def test_heldout_truth_cannot_change_fold_statistics(self):
        rows = rows_fixture()
        training = [row for row in rows if row["sample"] != "S0"]
        heldout = [dict(row) for row in rows if row["sample"] == "S0"]
        original = fit_runtime_context(training)
        for row in heldout:
            row["truth_allele1"], row["truth_allele2"] = "A*98:98", "A*99:99"
        self.assertEqual(original, fit_runtime_context(training))

    def test_fixed_universe_does_not_shrink(self):
        rows = [call("S1", "HLA-HD", ("A*01:01", "A*02:01")),
                call("S2", "HLA-HD", ("A*01:01", "A*02:01"))]
        calls = [{"sample": "S1", "modality": "wgs", "gene": "A", "is_correct": "1"}]
        universe = fixed_universe(calls, rows, "wgs", "test")
        self.assertEqual(len(universe), 2)
        self.assertEqual(sum(row["is_correct"] == "1" for row in universe), 1)

    def test_outer_fold_shard_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "at least one outer fold"):
                nested_run([], {}, None, None, {}, 3, 2, 7, Path(tmp), [])
            with self.assertRaisesRegex(ValueError, "must be unique"):
                nested_run([], {}, None, None, {}, 3, 2, 7, Path(tmp), [1, 1])
            with self.assertRaisesRegex(ValueError, "outside"):
                nested_run([], {}, None, None, {}, 3, 2, 7, Path(tmp), [3])


class ModelAndBundleTests(unittest.TestCase):
    def test_pair_symmetry_and_caller_permutation_invariance(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = embedding_store(Path(tmp))
            context = fit_runtime_context(rows_fixture())
            locus = build_loci(rows_fixture(), context, StubReference(), include_labels=True)[0]
            model = RefFormerModel(reference_dim=16, hidden_dim=16, layers=1, heads=4, dropout=0.0)
            model.eval()
            first = model._pair(("A*01:01", "A*02:01"), store, torch.device("cpu"))
            second = model._pair(("A*02:01", "A*01:01"), store, torch.device("cpu"))
            torch.testing.assert_close(first, second)
            scores1, _ = model.forward_locus(locus, store)
            scores2, _ = model.forward_locus(locus, store, list(reversed(range(5))))
            torch.testing.assert_close(scores1, scores2, rtol=1e-5, atol=1e-6)

    def test_training_prediction_checksum_and_abstention(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); store = embedding_store(tmp)
            context = fit_runtime_context(rows_fixture())
            loci = build_loci(rows_fixture(), context, StubReference(), include_labels=True)
            model, summary = train_model(loci[:6], loci[6:], store,
                {"hidden_dim": 16, "set_layers": 1, "heads": 4, "dropout": 0.0,
                 "epochs": 2, "patience": 3, "learning_rate": 0.001, "device": "cpu"}, 3)
            manifest = save_bundle(tmp / "bundle", model, context, store, 1.0,
                                   {"95": 0.9, "90": 0.8, "80": 0.7}, summary)
            self.assertIn("tokenizer", manifest)
            self.assertIn("reference_manifest", manifest)
            restored, restored_manifest = load_bundle(tmp / "bundle", store)
            calls, audit = predict_loci(restored, loci[:1], store, threshold=1.01,
                                         model_sha256=restored_manifest["bundle_sha256"])
            self.assertEqual(calls[0]["call_status"], "no_call")
            self.assertTrue(audit)
            weights = tmp / "bundle" / "model.pt"
            weights.write_bytes(weights.read_bytes() + b"tamper")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                load_bundle(tmp / "bundle", store)

    def test_calibration_and_coverage_operating_points(self):
        temperature = fit_temperature([([4.0, 1.0], 0), ([1.0, 3.0], 1)])
        self.assertGreaterEqual(temperature, 0.2)
        self.assertLessEqual(temperature, 4.0)
        thresholds = coverage_thresholds([{"confidence": value} for value in (0.1, 0.2, 0.3, 0.4, 0.5)])
        self.assertLessEqual(thresholds["95"], thresholds["90"])
        self.assertLessEqual(thresholds["90"], thresholds["80"])


class StandaloneFallbackTests(unittest.TestCase):
    def test_modality_specific_adapters_and_input_validation(self):
        encoder = PileupEvidenceEncoder(hidden_dim=16, layers=1, heads=4, dropout=0.0)
        self.assertEqual(set(encoder.adapters), {"wgs", "wes", "rnaseq"})
        self.assertIsNot(encoder.adapters["wgs"][1].weight, encoder.adapters["wes"][1].weight)
        pileup = torch.zeros((1, 12, 44))
        mask = torch.ones((1, 12), dtype=torch.bool)
        self.assertEqual(tuple(encoder(pileup, mask, "wgs").shape), (1, 16))
        with self.assertRaisesRegex(ValueError, "unsupported modality"):
            encoder(pileup, mask, "longread")

    def test_synthetic_pileup_is_deterministic_and_modality_aware(self):
        config = SyntheticPileupConfig(max_positions=32)
        first = sequence_pair_to_pileup(
            "ACGT" * 20, "AGGT" * 20, "wgs", torch.Generator().manual_seed(9), config
        )
        second = sequence_pair_to_pileup(
            "ACGT" * 20, "AGGT" * 20, "wgs", torch.Generator().manual_seed(9), config
        )
        torch.testing.assert_close(first[0], second[0])
        torch.testing.assert_close(first[1], second[1])
        self.assertEqual(tuple(first[0].shape), (1, 32, 44))
        self.assertTrue(bool(first[1].any()))
        wes, _ = sequence_pair_to_pileup(
            "ACGT" * 20, "AGGT" * 20, "wes", torch.Generator().manual_seed(9), config
        )
        self.assertFalse(torch.equal(first[0], wes))

    def test_pair_sampling_and_pair_scoring_are_unordered(self):
        candidates = sample_candidate_pairs(
            ["A*01:01", "A*02:01", "A*03:01", "A*24:02"],
            ("A*02:01", "A*01:01"), 4, random.Random(5),
        )
        self.assertEqual(candidates[0], ("A*01:01", "A*02:01"))
        self.assertEqual(len(candidates), len(set(candidates)))
        self.assertEqual(stable_pair_split(candidates[0], 7),
                         stable_pair_split(tuple(reversed(candidates[0])), 7))
        with tempfile.TemporaryDirectory() as tmp:
            store = embedding_store(Path(tmp))
            model = StandaloneReferenceRanker(
                reference_dim=16, hidden_dim=16, evidence_layers=1,
                evidence_heads=4, dropout=0.0,
            )
            evidence = torch.randn(16)
            first = model.score_pairs(evidence, [("A*01:01", "A*02:01")], store)
            second = model.score_pairs(evidence, [("A*02:01", "A*01:01")], store)
            torch.testing.assert_close(first, second)
            with self.assertRaisesRegex(ValueError, "duplicate unordered"):
                model.score_pairs(
                    evidence,
                    [("A*01:01", "A*02:01"), ("A*02:01", "A*01:01")],
                    store,
                )


class BaselineCompatibilityTests(unittest.TestCase):
    def test_existing_strict_baselines_reproduce(self):
        script = ROOT / "vendor" / "metaconsensus" / "verify_strict_baselines.py"
        bundled = ROOT / "inputs" / "baseline_oof_calls.tsv"
        analysis = Path(os.environ.get(
            "CHAMPHLA_STRICT_OOF",
            r"D:\Users\uonur\Downloads\ChampHLA_MetaConsensus_implementation\analysis\strict_reanalysis\oof_calls.tsv"
            if os.name == "nt"
            else "/mnt/d/Users/uonur/Downloads/ChampHLA_MetaConsensus_implementation/analysis/strict_reanalysis/oof_calls.tsv",
        ))
        oof = bundled if bundled.exists() else analysis
        result = subprocess.run(
            [sys.executable, str(script), "--oof", str(oof)],
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
