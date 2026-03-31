import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path("/users/ozcanumu/scratch/project_2008084/ozcanumu/repos/hla-typing-pipeline")
SCRIPT = REPO / "bin" / "hla_benchmark.py"
FIXTURES = REPO / "tests" / "fixtures"


class BenchmarkWorkflowTest(unittest.TestCase):
    def test_end_to_end_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "out"
            config = FIXTURES / "benchmark_config.yaml"
            subprocess.run(["python3", str(SCRIPT), "--config", str(config), "--output-dir", str(outdir)], check=True, cwd=str(REPO))

            summary_path = outdir / "tables" / "summary_full_cohort.tsv"
            self.assertTrue(summary_path.exists())
            with summary_path.open("r", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            summary = {(row["tool"], row["modality"]): row for row in rows}

            self.assertEqual(summary[("SpecHLA", "wes")]["overall_correct_call_rate"], "1.0")
            self.assertEqual(summary[("SpecHLA", "wgs")]["overall_correct_call_rate"], "1.0")
            self.assertEqual(summary[("OptiType", "wes")]["overall_correct_call_rate"], "1.0")
            self.assertEqual(summary[("OptiType", "wgs")]["overall_correct_call_rate"], "0.8333")
            self.assertEqual(summary[("OptiType", "rnaseq")]["overall_correct_call_rate"], "0.6667")
            self.assertEqual(summary[("OptiType", "rnaseq")]["callable_rate"], "0.8333")
            self.assertEqual(summary[("ArcasHLA", "rnaseq")]["overall_correct_call_rate"], "0.8333")
            self.assertEqual(summary[("ArcasHLA", "rnaseq")]["median_runtime_hours"], "")
            self.assertEqual(summary[("ArcasHLA", "rnaseq")]["median_max_ram_gb"], "")

            harmonized_path = outdir / "tables" / "harmonized_benchmark_rows.tsv"
            with harmonized_path.open("r", encoding="utf-8") as handle:
                harmonized = list(csv.DictReader(handle, delimiter="\t"))
            self.assertIn("rnaseq", {row["modality"] for row in harmonized})
            self.assertIn("confidence_score", harmonized[0])
            self.assertIn("confidence_source", harmonized[0])
            self.assertIn("is_correct_3field", harmonized[0])
            self.assertIn("match_grade", harmonized[0])
            self.assertEqual({row["imgt_hla_version"] for row in harmonized}, {"3.59.0"})
            optitype_conf_rows = [row for row in harmonized if row["tool"] == "OptiType" and row["modality"] == "wes"]
            self.assertTrue(all(row["confidence_source"] == "optitype_result_objective" for row in optitype_conf_rows))

            modality_gene_path = outdir / "tables" / "summary_modality_gene_coverage.tsv"
            self.assertTrue(modality_gene_path.exists())
            with modality_gene_path.open("r", encoding="utf-8") as handle:
                modality_gene_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertIn(("rnaseq", "A"), {(row["modality"], row["gene"]) for row in modality_gene_rows})

            ambiguity_path = outdir / "tables" / "ambiguity_summary.tsv"
            ambiguity_gene_path = outdir / "tables" / "ambiguity_summary_by_gene.tsv"
            reference_metadata_path = outdir / "tables" / "reference_metadata.tsv"
            self.assertTrue(ambiguity_path.exists())
            self.assertTrue(ambiguity_gene_path.exists())
            self.assertTrue(reference_metadata_path.exists())
            with ambiguity_path.open("r", encoding="utf-8") as handle:
                ambiguity_rows = list(csv.DictReader(handle, delimiter="	"))
            ambiguity_map = {(row["tool"], row["modality"]): row for row in ambiguity_rows}
            self.assertEqual(ambiguity_map[("OptiType", "wgs")]["exact_2field_rate"], "0.8333")
            self.assertEqual(ambiguity_map[("ArcasHLA", "rnaseq")]["exact_3field_rate"], "0.8333")
            self.assertEqual(ambiguity_map[("SpecHLA", "wes")]["imgt_hla_version"], "3.59.0")
            with reference_metadata_path.open("r", encoding="utf-8") as handle:
                reference_rows = list(csv.DictReader(handle, delimiter="	"))
            self.assertEqual(reference_rows[0]["imgt_hla_version"], "3.59.0")
            self.assertEqual(reference_rows[0]["secondary_resolutions"], "2,3")

            method_comparison_path = outdir / "tables" / "method_comparison.tsv"
            self.assertTrue(method_comparison_path.exists())
            with method_comparison_path.open("r", encoding="utf-8") as handle:
                method_rows = list(csv.DictReader(handle, delimiter="\t"))
            method_summary = {(row["method"], row["modality"]): row for row in method_rows}
            self.assertEqual(method_summary[("MajorityVote", "wgs")]["overall_correct_call_rate"], "0.8333")
            self.assertEqual(method_summary[("WeightedConsensus", "wgs")]["overall_correct_call_rate"], "1.0")
            self.assertEqual(method_summary[("WeightedConsensus", "rnaseq")]["overall_correct_call_rate"], "0.8333")

            weights_path = outdir / "tables" / "tool_confidence_weights.tsv"
            runtime_weights_path = outdir / "tables" / "consensus_runtime_weights.json"
            self.assertTrue(weights_path.exists())
            self.assertTrue(runtime_weights_path.exists())
            with weights_path.open("r", encoding="utf-8") as handle:
                weight_rows = list(csv.DictReader(handle, delimiter="\t"))
            weight_map = {(row["tool"], row["modality"]): row for row in weight_rows}
            self.assertEqual(weight_map[("SpecHLA", "wes")]["base_reliability"], "1.0")
            self.assertEqual(weight_map[("ArcasHLA", "rnaseq")]["confidence_coverage_rate"], "1.0")
            self.assertEqual(weight_map[("SpecHLA", "wgs")]["confidence_coverage_rate"], "0.0")

            calibration_path = outdir / "tables" / "confidence_calibration_summary.tsv"
            self.assertTrue(calibration_path.exists())
            with calibration_path.open("r", encoding="utf-8") as handle:
                calibration_rows = list(csv.DictReader(handle, delimiter="\t"))
            calibration_map = {(row["tool"], row["modality"]): row for row in calibration_rows}
            self.assertEqual(calibration_map[("OptiType", "wes")]["observed_accuracy"], "1")
            self.assertEqual(calibration_map[("OptiType", "wgs")]["observed_accuracy"], "0.8333")
            self.assertEqual(calibration_map[("ArcasHLA", "rnaseq")]["n_rows"], "6")

            discordance_path = outdir / "tables" / "discordance_summary.tsv"
            self.assertTrue(discordance_path.exists())
            with discordance_path.open("r", encoding="utf-8") as handle:
                discordance_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(discordance_rows, [{"scope": "cross_modality", "tag": "dna_rna_discordance", "n_events": "1"}])

            self.assertTrue((outdir / "figures" / "figure_2_accuracy_comparison.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_3_per_gene_gains.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_4_confidence_calibration.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_5_abstention_tradeoff.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_6_discordance_taxonomy.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_7_confidence_weights.svg").exists())


if __name__ == "__main__":
    unittest.main()
