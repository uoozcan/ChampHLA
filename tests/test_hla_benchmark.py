import csv
import importlib.util
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "bin" / "hla_benchmark.py"
FIXTURES = REPO / "tests" / "fixtures"


spec = importlib.util.spec_from_file_location("hb_module", str(SCRIPT))
hb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hb)


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
            self.assertEqual(method_summary[("WeightedConsensus", "rnaseq")]["overall_correct_call_rate"], "1.0")

            method_ambiguity_path = outdir / "tables" / "method_ambiguity_summary.tsv"
            method_ambiguity_gene_path = outdir / "tables" / "method_ambiguity_summary_by_gene.tsv"
            majority_calls_path = outdir / "tables" / "majority_vote_baseline.tsv"
            weighted_calls_path = outdir / "tables" / "weighted_consensus_calls.tsv"
            self.assertTrue(method_ambiguity_path.exists())
            self.assertTrue(method_ambiguity_gene_path.exists())
            self.assertTrue(majority_calls_path.exists())
            self.assertTrue(weighted_calls_path.exists())
            with method_ambiguity_path.open("r", encoding="utf-8") as handle:
                method_ambiguity_rows = list(csv.DictReader(handle, delimiter="	"))
            method_ambiguity_map = {(row["method"], row["modality"]): row for row in method_ambiguity_rows}
            self.assertEqual(method_ambiguity_map[("MajorityVote", "wgs")]["exact_3field_rate"], "0.8333")
            self.assertEqual(method_ambiguity_map[("WeightedConsensus", "wgs")]["exact_3field_rate"], "1.0")
            with weighted_calls_path.open("r", encoding="utf-8") as handle:
                weighted_rows = list(csv.DictReader(handle, delimiter="	"))
            weighted_map = {(row["sample"], row["modality"], row["gene"]): row for row in weighted_rows}
            self.assertIn("is_correct_3field", weighted_rows[0])
            self.assertEqual(weighted_map[("S2", "wgs", "C")]["is_correct_3field"], "1")
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
            confidence_error_path = outdir / "tables" / "confidence_error_summary.tsv"
            confidence_error_gene_path = outdir / "tables" / "confidence_error_summary_by_gene.tsv"
            self.assertTrue(confidence_error_path.exists())
            self.assertTrue(confidence_error_gene_path.exists())
            with confidence_error_path.open("r", encoding="utf-8") as handle:
                confidence_error_rows = list(csv.DictReader(handle, delimiter="\t"))
            confidence_error_map = {(row["tool"], row["modality"], row["bin_index"]): row for row in confidence_error_rows}
            self.assertEqual(confidence_error_map[("OptiType", "wgs", "5")]["error_rate"], "0.1667")
            self.assertEqual(confidence_error_map[("ArcasHLA", "rnaseq", "2")]["observed_accuracy"], "0.0")
            with confidence_error_gene_path.open("r", encoding="utf-8") as handle:
                confidence_error_gene_rows = list(csv.DictReader(handle, delimiter="\t"))
            confidence_error_gene_map = {
                (row["tool"], row["modality"], row["gene"], row["bin_index"]): row for row in confidence_error_gene_rows
            }
            self.assertEqual(confidence_error_gene_map[("OptiType", "wgs", "C", "5")]["error_rate"], "0.5")
            self.assertEqual(confidence_error_gene_map[("ArcasHLA", "rnaseq", "C", "2")]["error_count"], "1")

            discordance_path = outdir / "tables" / "discordance_summary.tsv"
            self.assertTrue(discordance_path.exists())
            with discordance_path.open("r", encoding="utf-8") as handle:
                discordance_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(discordance_rows, [{"scope": "rnaseq", "tag": "possible_expression_bias", "n_events": "1"}])

            self.assertTrue((outdir / "figures" / "figure_2_accuracy_comparison.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_3_per_gene_gains.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_4_confidence_calibration.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_5_abstention_tradeoff.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_6_discordance_taxonomy.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_7_confidence_weights.svg").exists())

    def test_calibration_guardrail_blocks_overconfident_boost(self):
        poor_entries = [
            {"is_callable": "1", "is_correct": "1", "confidence_score": "1.0", "sample": "S1", "gene": "A", "tool": "PoorTool", "modality": "wgs"},
            {"is_callable": "1", "is_correct": "0", "confidence_score": "1.0", "sample": "S2", "gene": "A", "tool": "PoorTool", "modality": "wgs"},
            {"is_callable": "1", "is_correct": "0", "confidence_score": "1.0", "sample": "S3", "gene": "A", "tool": "PoorTool", "modality": "wgs"},
            {"is_callable": "1", "is_correct": "0", "confidence_score": "1.0", "sample": "S4", "gene": "A", "tool": "PoorTool", "modality": "wgs"},
        ]
        settings = hb.confidence_guardrail_settings({})
        poor_record = hb.compute_weight_record(poor_entries, "PoorTool", "wgs", settings=settings)
        self.assertEqual(poor_record["base_reliability"], 0.25)
        self.assertEqual(poor_record["calibrated_confidence"], 1.0)
        self.assertEqual(poor_record["guardrail_status"], "poor_calibration")
        self.assertEqual(poor_record["guardrail_factor"], 0.0)
        self.assertEqual(poor_record["effective_confidence"], 0.25)
        self.assertEqual(poor_record["final_weight"], 0.25)

        good_entries = [
            {"is_callable": "1", "is_correct": "1", "confidence_score": "0.8", "sample": "S1", "gene": "A", "tool": "GoodTool", "modality": "wgs"},
            {"is_callable": "1", "is_correct": "1", "confidence_score": "0.7", "sample": "S2", "gene": "A", "tool": "GoodTool", "modality": "wgs"},
            {"is_callable": "1", "is_correct": "0", "confidence_score": "0.4", "sample": "S3", "gene": "A", "tool": "GoodTool", "modality": "wgs"},
            {"is_callable": "1", "is_correct": "0", "confidence_score": "0.3", "sample": "S4", "gene": "A", "tool": "GoodTool", "modality": "wgs"},
        ]
        good_record = hb.compute_weight_record(good_entries, "GoodTool", "wgs", settings=settings)
        self.assertEqual(good_record["guardrail_status"], "applied")
        self.assertGreater(good_record["guardrail_factor"], 0.0)
        self.assertGreater(good_record["effective_confidence"], good_record["base_reliability"])
        self.assertGreater(good_record["final_weight"], good_record["base_reliability"])

    def test_native_t1k_and_arcashla_confidence_parsers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            t1k_path = tmpdir / "S1_genotype.tsv"
            t1k_path.write_text("HLA-A\t2\tHLA-A*01:01:01\t55.0\t30\tHLA-A*02:01:01\t45.0\t20\nHLA-B\t1\tHLA-B*07:02:01\t50.0\t25\t.\t0\t-1\n", encoding="utf-8")
            arcas_path = tmpdir / "S1.genes.json"
            arcas_path.write_text('{"A": [66.0, 23, 0.25], "B": [103.0, 64, 0.5]}', encoding="utf-8")

            t1k_entries = hb.parse_confidence_file(t1k_path, "t1k_genotype_confidence", "S1", {"target_reads": 50})
            arcas_entries = hb.parse_confidence_file(arcas_path, "arcashla_genes_confidence", "S1", {"target_reads": 50})

            t1k_map = {entry["gene"]: entry for entry in t1k_entries}
            arcas_map = {entry["gene"]: entry for entry in arcas_entries}

            self.assertEqual(t1k_map["A"]["sample"], "S1")
            self.assertEqual(t1k_map["A"]["read_support"], 50.0)
            self.assertEqual(t1k_map["A"]["confidence_score"], 1.0)
            self.assertEqual(t1k_map["B"]["read_support"], 25.0)
            self.assertEqual(t1k_map["B"]["confidence_source"], "t1k_genotype_confidence")

            self.assertEqual(arcas_map["A"]["sample"], "S1")
            self.assertEqual(arcas_map["A"]["read_support"], 23.0)
            self.assertEqual(arcas_map["A"]["raw_confidence"], 0.25)
            self.assertEqual(arcas_map["A"]["confidence_source"], "arcashla_genes_confidence")
            self.assertEqual(arcas_map["B"]["confidence_score"], 0.5)

            hlahd_path = tmpdir / "S1_A.read.txt"
            hlahd_path.write_text("HLA-A*01:01:01:01\t251\nR1 only\t92\n", encoding="utf-8")
            hlahd_entries = hb.parse_confidence_file(hlahd_path, "hlahd_read_confidence", "S1", {"target_reads": 50})
            hlahd_map = {entry["gene"]: entry for entry in hlahd_entries}
            self.assertEqual(hlahd_map["A"]["sample"], "S1")
            self.assertEqual(hlahd_map["A"]["read_support"], 251.0)
            self.assertEqual(hlahd_map["A"]["confidence_score"], 1.0)
            self.assertEqual(hlahd_map["A"]["confidence_source"], "hlahd_read_confidence")

            kourami_path = tmpdir / "S1.kourami.result"
            kourami_path.write_text("A*01:01:01G\t540\t0.98\t546\t546\t15.0\t7.0\t8.0\nA*02:01:01G\t530\t0.94\t546\t546\t15.0\t7.0\t8.0\nB*07:02:01G\t300\t0.75\t320\t320\t10.0\t5.0\t5.0\n", encoding="utf-8")
            kourami_entries = hb.parse_confidence_file(kourami_path, "kourami_result_confidence", "S1", {"target_reads": 50})
            kourami_map = {entry["gene"]: entry for entry in kourami_entries}
            self.assertEqual(kourami_map["A"]["sample"], "S1")
            self.assertEqual(kourami_map["A"]["read_support"], 540.0)
            self.assertAlmostEqual(kourami_map["A"]["raw_confidence"], 0.96, places=6)
            self.assertAlmostEqual(kourami_map["A"]["confidence_score"], 0.96, places=6)
            self.assertEqual(kourami_map["B"]["confidence_source"], "kourami_result_confidence")



if __name__ == "__main__":
    unittest.main()
