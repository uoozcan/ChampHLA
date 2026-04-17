import importlib.util
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VOTING_SCRIPT = REPO / "bin" / "majority_voting.py"


spec = importlib.util.spec_from_file_location("majority_voting", VOTING_SCRIPT)
majority_voting = importlib.util.module_from_spec(spec)
spec.loader.exec_module(majority_voting)


class MajorityVotingRuntimeTest(unittest.TestCase):
    def test_runtime_schema_tool_weight_lookup(self):
        weights = {
            "tool_weights": {
                "optitype": {
                    "wes": {"final_weight": 0.82}
                }
            },
            "gene_weights": {}
        }
        row = {"tool": "optitype", "modality": "wes", "gene": "A"}
        self.assertAlmostEqual(majority_voting.lookup_weight(weights, row, use_gene_specific=False), 0.82)

    def test_runtime_schema_gene_weight_lookup(self):
        weights = {
            "tool_weights": {
                "optitype": {
                    "wes": {"final_weight": 0.50}
                }
            },
            "gene_weights": {
                "optitype": {
                    "wes": {
                        "A": {"final_weight": 0.91}
                    }
                }
            }
        }
        row = {"tool": "optitype", "modality": "wes", "gene": "A"}
        self.assertAlmostEqual(majority_voting.lookup_weight(weights, row, use_gene_specific=True), 0.91)

    def test_raw_accuracy_gene_specific_lookup(self):
        weights = {
            "raw_accuracy": {
                "optitype": {"A": 0.90, "B": 0.80, "C": 0.70, "DRB1": 0.60, "DQB1": 0.50}
            }
        }
        row = {"tool": "optitype", "modality": "wgs", "gene": "B"}
        self.assertAlmostEqual(majority_voting.lookup_weight(weights, row, use_gene_specific=True), 0.80)

    def test_raw_accuracy_tool_mean_lookup(self):
        weights = {
            "raw_accuracy": {
                "optitype": {"A": 0.90, "B": 0.80, "C": 0.70, "DRB1": 0.60, "DQB1": 0.50}
            }
        }
        row = {"tool": "optitype", "modality": "wgs", "gene": "A"}
        expected = (0.90 + 0.80 + 0.70 + 0.60 + 0.50) / 5.0
        self.assertAlmostEqual(majority_voting.lookup_weight(weights, row, use_gene_specific=False), expected)

    def test_raw_accuracy_missing_gene_defaults_zero(self):
        weights = {
            "raw_accuracy": {
                "optitype": {"A": 0.90}
            }
        }
        row = {"tool": "optitype", "modality": "wgs", "gene": "B"}
        self.assertEqual(majority_voting.lookup_weight(weights, row, use_gene_specific=True), 0.0)

    def test_equal_mode_bypasses_weight_schemas(self):
        weights = {"raw_accuracy": {"optitype": {"A": 0.1}}}
        row = {"tool": "missing_tool", "modality": "wes", "gene": "A"}
        self.assertEqual(majority_voting.lookup_weight(weights, row, use_gene_specific=False, equal_mode=True), 1.0)


class ModalityWiringTest(unittest.TestCase):
    def test_main_resolves_modality_from_effective_input(self):
        text = (REPO / "main.nf").read_text(encoding="utf-8")
        self.assertIn("def resolved_modality = params.run_modality", text)
        self.assertIn("(eff_input_type == 'fastq' ? 'wes' : 'wgs')", text)
        self.assertIn("(seq_type == 'rna') ? 'rnaseq'", text)

    def test_consensus_workflow_receives_resolved_modality(self):
        text = (REPO / "main.nf").read_text(encoding="utf-8")
        self.assertIn("MAJORITY_VOTING_WORKFLOW(", text)
        self.assertIn("resolved_modality", text)

    def test_majority_voting_module_uses_passed_modality(self):
        text = (REPO / "modules" / "majority_voting.nf").read_text(encoding="utf-8")
        self.assertIn("val modality", text)
        self.assertIn("AGGREGATE_RESULTS(ch_all_results.collect(), ch_modality)", text)
        self.assertIn("MAJORITY_VOTING(AGGREGATE_RESULTS.out.calls, ch_modality)", text)

    def test_default_consensus_genes_are_abc(self):
        text = (REPO / "nextflow.config").read_text(encoding="utf-8")
        self.assertIn("mv_genes                = 'A,B,C'", text)


if __name__ == "__main__":
    unittest.main()
