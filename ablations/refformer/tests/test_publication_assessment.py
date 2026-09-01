import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "build_publication_assessment.py"
SPEC = importlib.util.spec_from_file_location("build_publication_assessment", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class PublicationAssessmentTests(unittest.TestCase):
    def test_failed_gates_select_negative_ablation_route(self):
        candidate = {"passed": False}
        meta = {"passed": False}
        synthetic = {
            "synthetic_pretraining": {
                "alignment_aware_pilot": {
                    "best_validation_accuracy": 0.07,
                    "random_eight_candidate_expectation": 0.125,
                }
            }
        }
        imgt = {
            "available": True,
            "missing_tensor_count": 0,
            "test": {"accuracy": 0.10},
        }
        no_imgt = {"available": True, "test": {"accuracy": 0.12}}
        result = MODULE.determine_decision(candidate, meta, {}, synthetic, imgt, no_imgt)
        self.assertFalse(result["retain_refformer_as_method_headline"])
        self.assertFalse(result["scale_to_ten_outer_by_five_inner_folds"])
        self.assertEqual(result["recommendation"], "negative_ablation_in_benchmark_software_paper")

    def test_candidate_gate_can_retain_method(self):
        candidate = {"passed": True}
        meta = {"passed": False}
        synthetic = {
            "synthetic_pretraining": {
                "alignment_aware_pilot": {
                    "best_validation_accuracy": 0.0,
                    "random_eight_candidate_expectation": 0.125,
                }
            }
        }
        result = MODULE.determine_decision(candidate, meta, {}, synthetic, {}, {})
        self.assertTrue(result["retain_refformer_as_method_headline"])

    def test_meta_gate_is_retained_for_external_validation(self):
        candidate = {"passed": False}
        meta = {"passed": True}
        synthetic = {
            "synthetic_pretraining": {
                "alignment_aware_pilot": {
                    "best_validation_accuracy": 0.07,
                    "random_eight_candidate_expectation": 0.125,
                }
            }
        }
        result = MODULE.determine_decision(candidate, meta, {}, synthetic, {}, {})
        self.assertTrue(result["retain_metaconsensus_for_frozen_external_validation"])
        self.assertFalse(result["submission_ready_method_claim"])

    def test_external_nonregression_failure_selects_benchmark_route(self):
        candidate = {"passed": False}
        meta = {"passed": True}
        external = {
            "rna": {"available": True, "report": {
                "delta": -0.024, "delta_cluster_ci_lo": -0.06,
                "within_2pp_nonregression": False,
            }}
        }
        synthetic = {"synthetic_pretraining": {"alignment_aware_pilot": {
            "best_validation_accuracy": 0.07, "random_eight_candidate_expectation": 0.125,
        }}}
        result = MODULE.determine_decision(candidate, meta, external, synthetic, {}, {})
        self.assertFalse(result["retain_metaconsensus_as_method_headline"])
        self.assertIn("benchmark_software_paper", result["recommendation"])


if __name__ == "__main__":
    unittest.main()
