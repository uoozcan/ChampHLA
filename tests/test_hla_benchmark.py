import csv
import json
import importlib.util
import subprocess
import tempfile
import unittest
from pathlib import Path

import yaml

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
            self.assertIn("raw_score_family", harmonized[0])
            self.assertIn("raw_score_value", harmonized[0])
            self.assertIn("calibrated_probability", harmonized[0])
            self.assertIn("cv_calibrated_probability", harmonized[0])
            self.assertIn("is_ambiguity_compatible", harmonized[0])
            self.assertIn("is_resolution_compatible", harmonized[0])
            self.assertIn("compatibility_grade", harmonized[0])
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
            self.assertEqual(method_summary[("MajorityVote", "wgs")]["overall_correct_call_rate"], "1.0")
            self.assertEqual(method_summary[("WeightedConsensus", "wgs")]["overall_correct_call_rate"], "1.0")
            self.assertEqual(method_summary[("WeightedConsensus", "rnaseq")]["overall_correct_call_rate"], "1.0")
            meta_method_path = outdir / "tables" / "meta_method_comparison.tsv"
            self.assertTrue(meta_method_path.exists())
            with meta_method_path.open("r", encoding="utf-8") as handle:
                meta_method_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["method"] == "MetaConsensus" for row in meta_method_rows))

            method_ambiguity_path = outdir / "tables" / "method_ambiguity_summary.tsv"
            method_ambiguity_gene_path = outdir / "tables" / "method_ambiguity_summary_by_gene.tsv"
            majority_calls_path = outdir / "tables" / "majority_vote_baseline.tsv"
            weighted_calls_path = outdir / "tables" / "weighted_consensus_calls.tsv"
            self.assertTrue(method_ambiguity_path.exists())
            self.assertTrue(method_ambiguity_gene_path.exists())
            self.assertTrue(majority_calls_path.exists())
            self.assertTrue(weighted_calls_path.exists())
            self.assertTrue((outdir / "tables" / "meta_consensus_calls.tsv").exists())
            self.assertTrue((outdir / "tables" / "meta_consensus_decision_trace.tsv").exists())
            with method_ambiguity_path.open("r", encoding="utf-8") as handle:
                method_ambiguity_rows = list(csv.DictReader(handle, delimiter="	"))
            method_ambiguity_map = {(row["method"], row["modality"]): row for row in method_ambiguity_rows}
            self.assertEqual(method_ambiguity_map[("MajorityVote", "wgs")]["exact_3field_rate"], "1.0")
            self.assertEqual(method_ambiguity_map[("WeightedConsensus", "wgs")]["exact_3field_rate"], "1.0")
            with weighted_calls_path.open("r", encoding="utf-8") as handle:
                weighted_rows = list(csv.DictReader(handle, delimiter="	"))
            weighted_map = {(row["sample"], row["modality"], row["gene"]): row for row in weighted_rows}
            self.assertIn("is_correct_3field", weighted_rows[0])
            self.assertEqual(weighted_map[("S2", "wgs", "C")]["is_correct_3field"], "1")
            self.assertEqual(weighted_map[("S2", "wgs", "C")]["agreeing_tools"], "2")
            self.assertEqual(weighted_map[("S2", "wgs", "C")]["contributing_tools"], "5")
            self.assertIn("compatibility_grade", weighted_rows[0])
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
            self.assertEqual(weight_map[("SpecHLA", "wes")]["benchmark_mode"], "legacy_heuristic")
            self.assertIn("calibration_method", weight_rows[0])
            self.assertIn("mean_calibrated_probability", weight_rows[0])
            self.assertIn("calibrated_brier_score", weight_rows[0])

            calibration_path = outdir / "tables" / "confidence_calibration_summary.tsv"
            self.assertTrue(calibration_path.exists())
            with calibration_path.open("r", encoding="utf-8") as handle:
                calibration_rows = list(csv.DictReader(handle, delimiter="\t"))
            calibration_map = {(row["tool"], row["modality"]): row for row in calibration_rows}
            self.assertEqual(calibration_map[("OptiType", "wes")]["observed_accuracy"], "1")
            self.assertEqual(calibration_map[("OptiType", "wgs")]["observed_accuracy"], "0.8333")
            self.assertEqual(calibration_map[("ArcasHLA", "rnaseq")]["n_rows"], "6")
            self.assertEqual(calibration_map[("OptiType", "wes")]["benchmark_mode"], "legacy_heuristic")
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

            flags_path = outdir / "tables" / "sample_discordance_flags.tsv"
            self.assertTrue(flags_path.exists())
            with flags_path.open("r", encoding="utf-8") as handle:
                flags_rows = list(csv.DictReader(handle, delimiter="\t"))
            # Fixture has only possible_expression_bias (no dna_rna_discordance), so flags table should be empty
            self.assertEqual(flags_rows, [])

            self.assertTrue((outdir / "figures" / "figure_2_accuracy_comparison.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_3_per_gene_gains.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_4_confidence_calibration.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_5_abstention_tradeoff.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_6_discordance_taxonomy.svg").exists())
            self.assertTrue((outdir / "figures" / "figure_7_confidence_weights.svg").exists())
            self.assertTrue((outdir / "tables" / "tool_pairwise_agreement.tsv").exists())
            self.assertTrue((outdir / "tables" / "tool_pairwise_agreement_by_gene.tsv").exists())
            self.assertTrue((outdir / "tables" / "tool_vs_truth_error_taxonomy.tsv").exists())
            self.assertTrue((outdir / "tables" / "tool_disagreement_events.tsv").exists())
            self.assertTrue((outdir / "tables" / "consensus_decision_trace.tsv").exists())
            self.assertTrue((outdir / "tables" / "population_method_comparison.tsv").exists())
            self.assertTrue((outdir / "tables" / "locus_difficulty_summary.tsv").exists())
            self.assertTrue((outdir / "tables" / "summary_full_cohort_multiresolution.tsv").exists())
            self.assertTrue((outdir / "tables" / "summary_per_gene_multiresolution.tsv").exists())
            self.assertTrue((outdir / "tables" / "method_comparison_multiresolution.tsv").exists())
            self.assertTrue((outdir / "tables" / "population_resolution_summary.tsv").exists())
            self.assertTrue((outdir / "tables" / "population_conflict_summary.tsv").exists())

    def test_probabilistic_recalibration_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "out"
            config = yaml.safe_load((FIXTURES / "benchmark_config.yaml").read_text(encoding="utf-8"))
            config = hb.resolve_config_paths(config, FIXTURES)
            config.setdefault("benchmark", {})
            config["benchmark"]["mode"] = "probabilistic_recalibrated"
            config["benchmark"]["probabilistic_calibration"] = {"method": "platt", "cv_strategy": "loo"}
            config_path = tmpdir / "benchmark_prob.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--config", str(config_path), "--output-dir", str(outdir)], check=True, cwd=str(REPO))

            with (outdir / "tables" / "harmonized_benchmark_rows.tsv").open("r", encoding="utf-8") as handle:
                harmonized = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["calibrated_probability"] != "" for row in harmonized))
            self.assertTrue(any(row["cv_calibrated_probability"] != "" for row in harmonized))

            with (outdir / "tables" / "tool_confidence_weights.tsv").open("r", encoding="utf-8") as handle:
                weights = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(all(row["benchmark_mode"] == "probabilistic_recalibrated" for row in weights))
            self.assertTrue(any(row["mean_calibrated_probability"] != "" for row in weights))
            self.assertTrue(any(row["calibrated_brier_score"] != "" for row in weights))

            cv_summary = outdir / "tables" / "cross_validation_weight_summary.tsv"
            self.assertTrue(cv_summary.exists())
            with cv_summary.open("r", encoding="utf-8") as handle:
                cv_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["cv_mean_calibrated_probability"] != "" for row in cv_rows))

            runtime_weights = json.loads((outdir / "tables" / "consensus_runtime_weights.json").read_text(encoding="utf-8"))
            self.assertEqual(runtime_weights["benchmark_mode"], "probabilistic_recalibrated")

    def test_ensemble_ablation_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "out"
            config = yaml.safe_load((FIXTURES / "benchmark_config.yaml").read_text(encoding="utf-8"))
            config = hb.resolve_config_paths(config, FIXTURES)
            config.setdefault("benchmark", {})
            config["benchmark"]["ensemble_ablation"] = {
                "enabled": True,
                "tool_subsets": [
                    {"name": "all_tools", "tools": ["SpecHLA", "OptiType", "ArcasHLA"]},
                    {"name": "optitype_only", "tools": ["OptiType"]},
                ],
            }
            config_path = tmpdir / "benchmark_ablation.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--config", str(config_path), "--output-dir", str(outdir)], check=True, cwd=str(REPO))

            comparison_path = outdir / "tables" / "ablation_method_comparison.tsv"
            per_gene_path = outdir / "tables" / "ablation_method_per_gene.tsv"
            summary_path = outdir / "tables" / "ablation_summary.json"
            self.assertTrue(comparison_path.exists())
            self.assertTrue(per_gene_path.exists())
            self.assertTrue(summary_path.exists())

            with comparison_path.open("r", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["ablation_name"] == "all_tools" and row["method"] == "WeightedConsensus" for row in rows))
            self.assertTrue(any(row["ablation_name"] == "optitype_only" and row["method"] == "MetaConsensus" for row in rows))
            self.assertTrue(all("tool_subset" in row for row in rows))

            with per_gene_path.open("r", encoding="utf-8") as handle:
                per_gene_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["ablation_name"] == "all_tools" and row["gene"] == "A" for row in per_gene_rows))

            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertTrue(any(row["ablation_name"] == "all_tools" for row in payload))
            self.assertTrue(any(row["method"] == "MajorityVote" for row in payload))

    def test_locus_expert_consensus_settings_validation(self):
        rows = [
            {"tool": "OptiType", "modality": "wgs", "gene": "A"},
            {"tool": "T1K", "modality": "wgs", "gene": "A"},
            {"tool": "OptiType", "modality": "wgs", "gene": "B"},
            {"tool": "T1K", "modality": "wgs", "gene": "C"},
        ]
        valid = {
            "enabled": True,
            "parent_method": "weighted_consensus",
            "panel_sets": [{"name": "ok", "tool_subsets_by_gene": {"A": ["T1K"], "B": ["OptiType"], "C": ["T1K"]}}],
        }
        settings = hb.validate_locus_expert_consensus_settings(valid, rows, ["A", "B", "C"])
        self.assertEqual(settings["parent_method"], "weighted_consensus")

        invalid_parent = dict(valid, parent_method="bad_parent")
        with self.assertRaises(ValueError):
            hb.validate_locus_expert_consensus_settings(invalid_parent, rows, ["A", "B", "C"])

        missing_gene = {
            "enabled": True,
            "parent_method": "weighted_consensus",
            "panel_sets": [{"name": "missing", "tool_subsets_by_gene": {"A": ["T1K"], "B": ["OptiType"]}}],
        }
        with self.assertRaises(ValueError):
            hb.validate_locus_expert_consensus_settings(missing_gene, rows, ["A", "B", "C"])

        unknown_tool = {
            "enabled": True,
            "parent_method": "weighted_consensus",
            "panel_sets": [{"name": "unknown", "tool_subsets_by_gene": {"A": ["UnknownTool"], "B": ["OptiType"], "C": ["T1K"]}}],
        }
        with self.assertRaises(ValueError):
            hb.validate_locus_expert_consensus_settings(unknown_tool, rows, ["A", "B", "C"])

    def test_locus_expert_consensus_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "out"
            config = yaml.safe_load((FIXTURES / "benchmark_config.yaml").read_text(encoding="utf-8"))
            config = hb.resolve_config_paths(config, FIXTURES)
            config.setdefault("benchmark", {})
            config["benchmark"]["locus_expert_consensus"] = {
                "enabled": True,
                "parent_method": "weighted_consensus",
                "panel_sets": [
                    {
                        "name": "fixture_panel",
                        "tool_subsets_by_gene": {
                            "A": ["SpecHLA", "OptiType"],
                            "B": ["SpecHLA", "OptiType"],
                            "C": ["SpecHLA", "OptiType"],
                        },
                    }
                ],
            }
            config_path = tmpdir / "benchmark_locus_expert.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--config", str(config_path), "--output-dir", str(outdir)], check=True, cwd=str(REPO))

            calls_path = outdir / "tables" / "locus_expert_consensus_calls.tsv"
            trace_path = outdir / "tables" / "locus_expert_decision_trace.tsv"
            comparison_path = outdir / "tables" / "locus_expert_method_comparison.tsv"
            per_gene_path = outdir / "tables" / "locus_expert_method_per_gene.tsv"
            summary_path = outdir / "tables" / "locus_expert_summary.json"
            self.assertTrue(calls_path.exists())
            self.assertTrue(trace_path.exists())
            self.assertTrue(comparison_path.exists())
            self.assertTrue(per_gene_path.exists())
            self.assertTrue(summary_path.exists())

            with calls_path.open("r", encoding="utf-8") as handle:
                call_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["panel_name"] == "fixture_panel" and row["method"] == "LocusExpertConsensus" for row in call_rows))
            self.assertTrue(all("selected_tools" in row for row in call_rows))

            with comparison_path.open("r", encoding="utf-8") as handle:
                comparison_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["panel_name"] == "fixture_panel" and row["parent_method"] == "weighted_consensus" for row in comparison_rows))

    def test_champion_challenger_settings_validation(self):
        rows = [
            {"tool": "OptiType", "modality": "wgs", "gene": "A"},
            {"tool": "T1K", "modality": "wgs", "gene": "A"},
            {"tool": "OptiType", "modality": "wgs", "gene": "B"},
            {"tool": "OptiType", "modality": "wgs", "gene": "C"},
        ]
        valid = {
            "enabled": True,
            "champion_by_gene": {"A": "T1K", "B": "OptiType", "C": "OptiType"},
            "fallback_method": "weighted_consensus",
            "override_policy": {"min_challenger_support_fraction": 0.65, "min_challenger_margin": 0.20, "min_supporting_tools": 2, "require_non_ambiguity_override": True},
        }
        settings = hb.validate_champion_challenger_settings(valid, rows, ["A", "B", "C"])
        self.assertEqual(settings["fallback_method"], "weighted_consensus")

        invalid_fallback = dict(valid, fallback_method="majority_vote")
        with self.assertRaises(ValueError):
            hb.validate_champion_challenger_settings(invalid_fallback, rows, ["A", "B", "C"])

        missing_gene = dict(valid, champion_by_gene={"A": "T1K", "B": "OptiType"})
        with self.assertRaises(ValueError):
            hb.validate_champion_challenger_settings(missing_gene, rows, ["A", "B", "C"])

        unknown_tool = dict(valid, champion_by_gene={"A": "UnknownTool", "B": "OptiType", "C": "OptiType"})
        with self.assertRaises(ValueError):
            hb.validate_champion_challenger_settings(unknown_tool, rows, ["A", "B", "C"])

    def test_champion_challenger_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "out"
            config = yaml.safe_load((FIXTURES / "benchmark_config.yaml").read_text(encoding="utf-8"))
            config = hb.resolve_config_paths(config, FIXTURES)
            config.setdefault("benchmark", {})
            config["benchmark"]["champion_challenger"] = {
                "enabled": True,
                "champion_by_gene": {"A": "OptiType", "B": "OptiType", "C": "OptiType"},
                "fallback_method": "weighted_consensus",
                "override_policy": {
                    "min_challenger_support_fraction": 0.65,
                    "min_challenger_margin": 0.20,
                    "min_supporting_tools": 2,
                    "require_non_ambiguity_override": True,
                },
            }
            config_path = tmpdir / "benchmark_champion.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--config", str(config_path), "--output-dir", str(outdir)], check=True, cwd=str(REPO))

            calls_path = outdir / "tables" / "champion_challenger_calls.tsv"
            trace_path = outdir / "tables" / "champion_challenger_trace.tsv"
            comparison_path = outdir / "tables" / "champion_challenger_method_comparison.tsv"
            per_gene_path = outdir / "tables" / "champion_challenger_method_per_gene.tsv"
            summary_path = outdir / "tables" / "champion_challenger_summary.json"
            self.assertTrue(calls_path.exists())
            self.assertTrue(trace_path.exists())
            self.assertTrue(comparison_path.exists())
            self.assertTrue(per_gene_path.exists())
            self.assertTrue(summary_path.exists())

            with calls_path.open("r", encoding="utf-8") as handle:
                call_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["method"] == "ChampionChallenger" for row in call_rows))
            self.assertTrue(all("champion_tool" in row for row in call_rows))

            with comparison_path.open("r", encoding="utf-8") as handle:
                comparison_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["method"] == "ChampionChallenger" for row in comparison_rows))

    def test_consensus_gating_settings_validation(self):
        config = {
            "benchmark": {
                "champion_challenger": {
                    "enabled": True,
                    "champion_by_gene": {"A": "T1K", "B": "OptiType", "C": "OptiType"},
                    "fallback_method": "weighted_consensus",
                },
                "consensus_gating": {
                    "enabled": True,
                    "easy_policy": "majority_vote",
                    "hard_policy": "champion_challenger",
                    "hard_locus_rule": {
                        "min_distinct_pairs_for_hard": 4,
                        "min_support_margin_for_easy": 0.35,
                        "min_support_fraction_for_easy": 0.50,
                        "trigger_on_majority_weighted_disagreement": True,
                        "trigger_on_duplicated_top_pair_with_alternative": True,
                    },
                },
            }
        }
        settings = hb.validate_consensus_gating_settings(hb.consensus_gating_settings(config), config)
        self.assertEqual(settings["easy_policy"], "majority_vote")

        bad_easy = yaml.safe_load(yaml.safe_dump(config))
        bad_easy["benchmark"]["consensus_gating"]["easy_policy"] = "weighted_consensus"
        with self.assertRaises(ValueError):
            hb.validate_consensus_gating_settings(hb.consensus_gating_settings(bad_easy), bad_easy)

        bad_hard = yaml.safe_load(yaml.safe_dump(config))
        bad_hard["benchmark"]["consensus_gating"]["hard_policy"] = "weighted_consensus"
        with self.assertRaises(ValueError):
            hb.validate_consensus_gating_settings(hb.consensus_gating_settings(bad_hard), bad_hard)

        missing_champion = {"benchmark": {"consensus_gating": {"enabled": True, "easy_policy": "majority_vote", "hard_policy": "champion_challenger"}}}
        with self.assertRaises(ValueError):
            hb.validate_consensus_gating_settings(hb.consensus_gating_settings(missing_champion), missing_champion)

    def test_hard_locus_classification_and_duplicate_pair(self):
        self.assertTrue(hb.is_duplicated_pair(("A*01:01", "A*01:01")))
        self.assertFalse(hb.is_duplicated_pair(("A*01:01", "A*02:01")))
        settings = {
            "hard_locus_rule": {
                "min_distinct_pairs_for_hard": 4,
                "min_support_margin_for_easy": 0.35,
                "min_support_fraction_for_easy": 0.50,
                "trigger_on_majority_weighted_disagreement": True,
                "trigger_on_duplicated_top_pair_with_alternative": True,
            }
        }
        hard_features = {
            "distinct_pair_count": 4,
            "weighted_support_margin": 0.60,
            "weighted_support_fraction": 0.80,
            "majority_weighted_disagree": False,
            "weighted_top_pair_is_duplicated": False,
            "alternative_nonduplicated_exists": False,
        }
        difficulty, reasons = hb.classify_hard_locus(hard_features, settings)
        self.assertEqual(difficulty, "hard")
        self.assertEqual(reasons, ["distinct_pairs"])
        mixed_features = dict(hard_features, distinct_pair_count=2, weighted_support_margin=0.10, weighted_support_fraction=0.20, majority_weighted_disagree=True)
        difficulty, reasons = hb.classify_hard_locus(mixed_features, settings)
        self.assertEqual(difficulty, "hard")
        self.assertEqual(reasons, ["low_margin", "low_support", "majority_weighted_disagree"])

    def test_gated_consensus_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "out"
            config = yaml.safe_load((FIXTURES / "benchmark_config.yaml").read_text(encoding="utf-8"))
            config = hb.resolve_config_paths(config, FIXTURES)
            config.setdefault("benchmark", {})
            config["benchmark"]["champion_challenger"] = {
                "enabled": True,
                "champion_by_gene": {"A": "OptiType", "B": "OptiType", "C": "OptiType"},
                "fallback_method": "weighted_consensus",
                "override_policy": {
                    "min_challenger_support_fraction": 0.65,
                    "min_challenger_margin": 0.20,
                    "min_supporting_tools": 2,
                    "require_non_ambiguity_override": True,
                },
            }
            config["benchmark"]["consensus_gating"] = {
                "enabled": True,
                "easy_policy": "majority_vote",
                "hard_policy": "champion_challenger",
                "hard_locus_rule": {
                    "min_distinct_pairs_for_hard": 4,
                    "min_support_margin_for_easy": 0.35,
                    "min_support_fraction_for_easy": 0.50,
                    "trigger_on_majority_weighted_disagreement": True,
                    "trigger_on_duplicated_top_pair_with_alternative": True,
                },
            }
            config_path = tmpdir / "benchmark_gated.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--config", str(config_path), "--output-dir", str(outdir)], check=True, cwd=str(REPO))

            calls_path = outdir / "tables" / "gated_consensus_calls.tsv"
            trace_path = outdir / "tables" / "gated_consensus_trace.tsv"
            comparison_path = outdir / "tables" / "gated_method_comparison.tsv"
            per_gene_path = outdir / "tables" / "gated_method_per_gene.tsv"
            summary_path = outdir / "tables" / "gated_summary.json"
            self.assertTrue(calls_path.exists())
            self.assertTrue(trace_path.exists())
            self.assertTrue(comparison_path.exists())
            self.assertTrue(per_gene_path.exists())
            self.assertTrue(summary_path.exists())

            with calls_path.open("r", encoding="utf-8") as handle:
                call_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["method"] == "GatedConsensus" for row in call_rows))
            self.assertTrue(all(row["difficulty_class"] in {"easy", "hard"} for row in call_rows))
            self.assertTrue(all(row["selected_policy"] in {"majority_vote", "champion_challenger"} for row in call_rows))

    def test_threshold_sweep_settings_and_override_summary(self):
        config = {
            "benchmark_analysis": {
                "threshold_sweeps": {
                    "enabled": True,
                    "weighted_consensus": {
                        "min_support_values": [0.35, 0.55],
                        "min_margin_values": [0.0, 0.15],
                    },
                    "champion_challenger": {
                        "champion_by_gene": {"A": "OptiType", "B": "OptiType", "C": "OptiType"},
                        "min_challenger_support_fraction_values": [0.20, 0.65],
                        "min_challenger_margin_values": [0.0, 0.20],
                        "min_supporting_tools_values": [1, 2],
                        "require_non_ambiguity_override": True,
                    },
                }
            }
        }
        weighted = hb.weighted_threshold_sweep_settings(config)
        champion = hb.champion_override_sweep_settings(config)
        self.assertTrue(weighted["enabled"])
        self.assertEqual(len(weighted["min_support_values"]) * len(weighted["min_margin_values"]), 4)
        self.assertEqual(len(champion["min_challenger_support_fraction_values"]) * len(champion["min_challenger_margin_values"]) * len(champion["min_supporting_tools_values"]), 8)

        truth_index = {("S1", "wes", "A"): ("A*01:01", "A*02:01")}
        trace_rows = [
            {"override_triggered": "1", "sample": "S1", "modality": "wes", "gene": "A", "champion_pair": "A*01:01+A*01:01", "challenger_pair": "A*01:01+A*02:01"},
            {"override_triggered": "1", "sample": "S1", "modality": "wes", "gene": "A", "champion_pair": "A*01:01+A*02:01", "challenger_pair": "A*01:01+A*01:01"},
            {"override_triggered": "1", "sample": "S1", "modality": "wes", "gene": "A", "champion_pair": "A*01:01+A*01:01", "challenger_pair": "A*03:01+A*03:01"},
        ]
        summary = hb.summarize_override_effects(trace_rows, truth_index)
        self.assertEqual(summary["override_count"], 3)
        self.assertEqual(summary["corrective_override_count"], 1)
        self.assertEqual(summary["harmful_override_count"], 1)
        self.assertEqual(summary["neutral_override_count"], 1)

    def test_threshold_sweep_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "out"
            config = yaml.safe_load((FIXTURES / "benchmark_config.yaml").read_text(encoding="utf-8"))
            config = hb.resolve_config_paths(config, FIXTURES)
            config.setdefault("benchmark_analysis", {})
            config["benchmark_analysis"]["threshold_sweeps"] = {
                "enabled": True,
                "weighted_consensus": {
                    "min_support_values": [0.35, 0.55],
                    "min_margin_values": [0.0, 0.15],
                },
                "champion_challenger": {
                    "champion_by_gene": {"A": "OptiType", "B": "OptiType", "C": "OptiType"},
                    "min_challenger_support_fraction_values": [0.20, 0.65],
                    "min_challenger_margin_values": [0.0, 0.20],
                    "min_supporting_tools_values": [1, 2],
                    "require_non_ambiguity_override": True,
                },
            }
            config_path = tmpdir / "benchmark_sweeps.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--config", str(config_path), "--output-dir", str(outdir)], check=True, cwd=str(REPO))

            weighted_path = outdir / "sweeps" / "wes_weighted_threshold_sweep.tsv"
            champion_path = outdir / "sweeps" / "wes_champion_override_sweep.tsv"
            summary_path = outdir / "sweeps" / "wes_sweep_summary.json"
            self.assertTrue(weighted_path.exists())
            self.assertTrue(champion_path.exists())
            self.assertTrue(summary_path.exists())

            with weighted_path.open("r", encoding="utf-8") as handle:
                weighted_rows = list(csv.DictReader(handle, delimiter="\t"))
            with champion_path.open("r", encoding="utf-8") as handle:
                champion_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(weighted_rows), 4)
            self.assertEqual(len(champion_rows), 8)

    def test_threshold_sweep_outputs_rnaseq_with_forced_spechla(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "out"
            config = yaml.safe_load((FIXTURES / "benchmark_config.yaml").read_text(encoding="utf-8"))
            config = hb.resolve_config_paths(config, FIXTURES)
            config["runs"] = [
                dict(run)
                for run in config["runs"]
                if run.get("modality") == "rnaseq" or (run.get("tool") == "SpecHLA" and run.get("modality") == "wes")
            ]
            for run in config["runs"]:
                if run.get("tool") == "SpecHLA":
                    run["modality"] = "rnaseq"
                    run.pop("coverage_only", None)
            config.setdefault("benchmark_analysis", {})
            config["benchmark_analysis"]["threshold_sweeps"] = {
                "enabled": True,
                "weighted_consensus": {
                    "min_support_values": [0.35, 0.55],
                    "min_margin_values": [0.0, 0.15],
                },
                "champion_challenger": {
                    "champion_by_gene": {"A": "OptiType", "B": "ArcasHLA", "C": "ArcasHLA"},
                    "min_challenger_support_fraction_values": [0.20, 0.65],
                    "min_challenger_margin_values": [0.0, 0.20],
                    "min_supporting_tools_values": [1, 2],
                    "require_non_ambiguity_override": True,
                },
            }
            config_path = tmpdir / "benchmark_sweeps_rna.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

            subprocess.run(["python3", str(SCRIPT), "--config", str(config_path), "--output-dir", str(outdir)], check=True, cwd=str(REPO))

            weighted_path = outdir / "sweeps" / "rna_weighted_threshold_sweep.tsv"
            champion_path = outdir / "sweeps" / "rna_champion_override_sweep.tsv"
            summary_path = outdir / "sweeps" / "rna_sweep_summary.json"
            self.assertTrue(weighted_path.exists())
            self.assertTrue(champion_path.exists())
            self.assertTrue(summary_path.exists())

            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(summary["modality"], "rnaseq")
            self.assertTrue(summary["forced_spechla_accuracy"])
            self.assertTrue(summary["coverage_notes"]["spechla_forced_into_accuracy"])

            with weighted_path.open("r", encoding="utf-8") as handle:
                weighted_rows = list(csv.DictReader(handle, delimiter="\t"))
            with champion_path.open("r", encoding="utf-8") as handle:
                champion_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(weighted_rows), 4)
            self.assertEqual(len(champion_rows), 8)

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


class SampleIntegrityTest(unittest.TestCase):
    """Unit tests for build_sample_discordance_flags()."""

    def _discord(self, sample, gene):
        return {"sample": sample, "gene": gene, "scope": "cross_modality",
                "tag": "dna_rna_discordance", "detail": "rna_pair_differs_from_dna"}

    def _bias(self, sample, gene):
        return {"sample": sample, "gene": gene, "scope": "rnaseq",
                "tag": "possible_expression_bias", "detail": "rna_missing_dna_present"}

    def test_critical_flag_three_loci(self):
        rows = [self._discord("NA07000", g) for g in ["A", "B", "C"]]
        result = hb.build_sample_discordance_flags(rows)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["integrity_flag"], "critical")
        self.assertEqual(result[0]["n_discordant_loci"], 3)

    def test_warning_flag_two_loci(self):
        rows = [self._discord("S1", "A"), self._discord("S1", "B")]
        result = hb.build_sample_discordance_flags(rows)
        self.assertEqual(result[0]["integrity_flag"], "warning")
        self.assertEqual(result[0]["n_discordant_loci"], 2)

    def test_nominal_flag_one_locus(self):
        rows = [self._discord("S1", "A")]
        result = hb.build_sample_discordance_flags(rows)
        self.assertEqual(result[0]["integrity_flag"], "nominal")
        self.assertEqual(result[0]["n_discordant_loci"], 1)

    def test_expression_bias_only_excluded(self):
        rows = [self._bias("S1", g) for g in ["A", "B", "C"]]
        result = hb.build_sample_discordance_flags(rows)
        self.assertEqual(result, [])

    def test_discordant_loci_field_sorted(self):
        # Insert genes out of canonical order; expect sorted output A,B,C not C,A,B
        rows = [self._discord("S1", g) for g in ["C", "A", "B"]]
        result = hb.build_sample_discordance_flags(rows)
        self.assertEqual(result[0]["discordant_loci"], "A,B,C")

    def test_multiple_samples(self):
        rows = [
            self._discord("NA07000", "A"), self._discord("NA07000", "B"),
            self._discord("NA07000", "C"), self._bias("NA07000", "DRB1"),
            self._discord("HG00096", "B"), self._bias("HG00096", "A"),
        ]
        result = hb.build_sample_discordance_flags(rows)
        self.assertEqual(len(result), 2)
        by_sample = {row["sample"]: row for row in result}
        self.assertEqual(by_sample["NA07000"]["integrity_flag"], "critical")
        self.assertEqual(by_sample["NA07000"]["discordant_loci"], "A,B,C")
        self.assertEqual(by_sample["NA07000"]["n_expression_bias_loci"], 1)
        self.assertEqual(by_sample["NA07000"]["expression_bias_loci"], "DRB1")
        self.assertEqual(by_sample["HG00096"]["integrity_flag"], "nominal")
        self.assertEqual(by_sample["HG00096"]["discordant_loci"], "B")
        self.assertEqual(by_sample["HG00096"]["n_expression_bias_loci"], 1)
        self.assertEqual(by_sample["HG00096"]["expression_bias_loci"], "A")


class ArbitrationTest(unittest.TestCase):
    """Unit tests for arbitrate_dna_rna_discordance() and related helpers."""

    def _dna(self, callable_tools=4, support_fraction=0.80, support_margin=0.30, pair=("A*01:01", "A*02:01")):
        return {"pair": pair, "callable_tools": callable_tools, "support_fraction": support_fraction,
                "support_margin": support_margin, "mean_confidence": None, "mean_read_support": None}

    def _rna(self, callable_tools=3, support_fraction=0.75, support_margin=0.25, pair=("A*03:01", "A*02:01"),
             mean_read_support=25.0):
        return {"pair": pair, "callable_tools": callable_tools, "support_fraction": support_fraction,
                "support_margin": support_margin, "mean_confidence": None, "mean_read_support": mean_read_support}

    def test_rule1_rna_no_callable_tools(self):
        result = hb.arbitrate_dna_rna_discordance(self._dna(), self._rna(callable_tools=0))
        self.assertEqual(result["arbitration_rule"], "rna_no_callable_tools")
        self.assertEqual(result["arbitration_outcome"], "dna_wins")
        self.assertEqual(result["arbitrated_pair"], ("A*01:01", "A*02:01"))

    def test_rule2_rna_low_read_support(self):
        result = hb.arbitrate_dna_rna_discordance(self._dna(), self._rna(mean_read_support=5.0))
        self.assertEqual(result["arbitration_rule"], "rna_low_read_support")
        self.assertEqual(result["arbitration_outcome"], "dna_wins")

    def test_rule2_skipped_when_read_support_none(self):
        # Rule 2 must not fire when mean_read_support is None
        result = hb.arbitrate_dna_rna_discordance(
            self._dna(callable_tools=4, support_fraction=0.80),
            self._rna(callable_tools=3, support_fraction=0.75, mean_read_support=None),
        )
        self.assertNotEqual(result["arbitration_rule"], "rna_low_read_support")

    def test_rule3_dna_insufficient_evidence(self):
        result = hb.arbitrate_dna_rna_discordance(
            self._dna(callable_tools=1, support_fraction=0.60),
            self._rna(callable_tools=3, support_fraction=0.80),
        )
        self.assertEqual(result["arbitration_rule"], "dna_insufficient_evidence")
        self.assertEqual(result["arbitration_outcome"], "rna_wins")

    def test_rule4_dna_stronger_consensus(self):
        result = hb.arbitrate_dna_rna_discordance(
            self._dna(callable_tools=4, support_fraction=0.85),
            self._rna(callable_tools=3, support_fraction=0.40),
        )
        self.assertEqual(result["arbitration_rule"], "dna_stronger_consensus")
        self.assertEqual(result["arbitration_outcome"], "dna_wins")

    def test_rule5_rna_stronger_consensus(self):
        result = hb.arbitrate_dna_rna_discordance(
            self._dna(callable_tools=4, support_fraction=0.45),
            self._rna(callable_tools=3, support_fraction=0.75),
        )
        self.assertEqual(result["arbitration_rule"], "rna_stronger_consensus")
        self.assertEqual(result["arbitration_outcome"], "rna_wins")

    def test_rule6_abstain_high_confidence_conflict(self):
        result = hb.arbitrate_dna_rna_discordance(
            self._dna(callable_tools=4, support_fraction=0.60),
            self._rna(callable_tools=3, support_fraction=0.62, mean_read_support=30.0),
        )
        self.assertEqual(result["arbitration_rule"], "high_confidence_conflict")
        self.assertEqual(result["arbitration_outcome"], "abstain")
        self.assertIsNone(result["arbitrated_pair"])

    def test_rule1_priority_over_rule3(self):
        # Even when rule 3's condition is met (dna_callable=1, rna_frac=0.90),
        # rule 1 fires first because rna_callable=0.
        result = hb.arbitrate_dna_rna_discordance(
            self._dna(callable_tools=1, support_fraction=0.60),
            self._rna(callable_tools=0, support_fraction=0.90),
        )
        self.assertEqual(result["arbitration_rule"], "rna_no_callable_tools")

    def test_build_arbitration_accuracy_table_correct_call(self):
        discordance_rows = [{
            "sample": "S1", "gene": "A", "scope": "cross_modality",
            "tag": "dna_rna_discordance", "detail": "rna_pair_differs_from_dna",
            "arbitrated_pair": "A*01:01+A*02:01",
            "arbitration_rule": "dna_stronger_consensus",
            "arbitration_outcome": "dna_wins",
        }]
        harmonized_rows = [{
            "sample": "S1", "gene": "A", "modality": "wes",
            "truth_allele1": "A*01:01", "truth_allele2": "A*02:01",
        }]
        result = hb.build_arbitration_accuracy_table(discordance_rows, harmonized_rows)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["is_correct"], "1")

    def test_build_arbitration_accuracy_table_incorrect_call(self):
        discordance_rows = [{
            "sample": "S1", "gene": "A", "scope": "cross_modality",
            "tag": "dna_rna_discordance", "detail": "rna_pair_differs_from_dna",
            "arbitrated_pair": "A*03:01+A*02:01",
            "arbitration_rule": "rna_stronger_consensus",
            "arbitration_outcome": "rna_wins",
        }]
        harmonized_rows = [{
            "sample": "S1", "gene": "A", "modality": "wes",
            "truth_allele1": "A*01:01", "truth_allele2": "A*02:01",
        }]
        result = hb.build_arbitration_accuracy_table(discordance_rows, harmonized_rows)
        self.assertEqual(result[0]["is_correct"], "0")

    def test_build_arbitration_accuracy_table_abstain(self):
        discordance_rows = [{
            "sample": "S1", "gene": "A", "scope": "cross_modality",
            "tag": "dna_rna_discordance", "detail": "rna_pair_differs_from_dna",
            "arbitrated_pair": "",
            "arbitration_rule": "high_confidence_conflict",
            "arbitration_outcome": "abstain",
        }]
        harmonized_rows = [{
            "sample": "S1", "gene": "A", "modality": "wes",
            "truth_allele1": "A*01:01", "truth_allele2": "A*02:01",
        }]
        result = hb.build_arbitration_accuracy_table(discordance_rows, harmonized_rows)
        self.assertEqual(result[0]["is_correct"], "")

    def test_build_arbitration_accuracy_table_skips_non_discordance_rows(self):
        discordance_rows = [
            {"sample": "S1", "gene": "A", "scope": "rnaseq", "tag": "possible_expression_bias", "detail": ""},
        ]
        result = hb.build_arbitration_accuracy_table(discordance_rows, [])
        self.assertEqual(result, [])

    def test_summarize_arbitration_rules_counts(self):
        rows = [
            {"arbitration_rule": "dna_stronger_consensus", "arbitration_outcome": "dna_wins", "is_correct": "1"},
            {"arbitration_rule": "dna_stronger_consensus", "arbitration_outcome": "dna_wins", "is_correct": "0"},
            {"arbitration_rule": "high_confidence_conflict", "arbitration_outcome": "abstain", "is_correct": ""},
            {"arbitration_rule": "rna_stronger_consensus", "arbitration_outcome": "rna_wins", "is_correct": "1"},
        ]
        summary = {r["arbitration_rule"]: r for r in hb.summarize_arbitration_rules(rows)}
        self.assertEqual(summary["dna_stronger_consensus"]["n_resolved"], 2)
        self.assertEqual(summary["dna_stronger_consensus"]["n_correct"], 1)
        self.assertEqual(summary["high_confidence_conflict"]["n_abstain"], 1)
        self.assertEqual(summary["high_confidence_conflict"]["n_resolved"], 0)
        self.assertEqual(summary["rna_stronger_consensus"]["n_correct"], 1)

    def test_modality_mean_confidence_callable_only(self):
        rows = [
            {"is_callable": "1", "confidence_score": "0.9"},
            {"is_callable": "0", "confidence_score": "0.1"},  # should be excluded
            {"is_callable": "1", "confidence_score": "0.7"},
        ]
        result = hb._modality_mean_confidence(rows)
        self.assertAlmostEqual(result, 0.8, places=3)

    def test_modality_mean_read_support_none_when_empty(self):
        result = hb._modality_mean_read_support([])
        self.assertIsNone(result)


_ciwd_spec = importlib.util.spec_from_file_location("ciwd_module", str(REPO / "bin" / "ciwd.py"))
ciwd = importlib.util.module_from_spec(_ciwd_spec)
_ciwd_spec.loader.exec_module(ciwd)


class CiwdCatalogueTest(unittest.TestCase):
    def setUp(self):
        self.cat = ciwd.load_ciwd()

    def test_catalogue_loads(self):
        self.assertGreater(len(self.cat), 3000)  # ~3249 two-field alleles in CIWD 3.0.0

    def test_common_intermediate_wd_notciwd_categories(self):
        # verified against the published CIWD 3.0.0 P-group table
        self.assertEqual(self.cat.category("A*02:01"), "common")
        self.assertEqual(self.cat.category("B*07:02"), "common")
        self.assertEqual(self.cat.category("DRB1*15:01"), "common")
        self.assertEqual(self.cat.category("A*01:04N"), "well_documented")
        self.assertEqual(self.cat.category("A*02:25"), "not_ciwd")  # 'o' in source

    def test_absent_allele_is_unknown(self):
        self.assertEqual(self.cat.category("A*99:99"), "unknown")

    def test_expression_suffix_normalised(self):
        # a null-suffixed truth allele resolves to its two-field catalogue entry
        self.assertEqual(self.cat.category("A*01:04"), self.cat.category("A*01:04N"))
        self.assertEqual(self.cat.category("HLA-A*02:01:01"), "common")

    def test_population_group_lookup(self):
        self.assertEqual(self.cat.category("A*02:01", "EURO"), "common")

    def test_genotype_stratum_is_rarer_allele(self):
        self.assertEqual(self.cat.genotype_stratum("A*02:01", "A*01:01"), "common")
        self.assertEqual(self.cat.genotype_stratum("A*02:01", "A*02:25"), "not_ciwd")
        self.assertEqual(self.cat.genotype_stratum("A*02:01", "A*99:99"), "unknown")

    def test_is_implausible_flags_rare_and_novel(self):
        self.assertFalse(self.cat.is_implausible("A*02:01"))
        self.assertTrue(self.cat.is_implausible("A*02:25"))
        self.assertTrue(self.cat.is_implausible("A*99:99"))

    def test_missing_table_degrades_to_unknown(self):
        empty = ciwd.load_ciwd(path=REPO / "assets" / "does_not_exist.tsv")
        self.assertEqual(len(empty), 0)
        self.assertEqual(empty.category("A*02:01"), "unknown")


class CiwdBenchmarkOutputTest(unittest.TestCase):
    def test_stratified_and_plausibility_tables_written(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "out"
            config = FIXTURES / "benchmark_config.yaml"
            subprocess.run(["python3", str(SCRIPT), "--config", str(config), "--output-dir", str(outdir)], check=True, cwd=str(REPO))

            strat_path = outdir / "tables" / "summary_ciwd_stratified.tsv"
            self.assertTrue(strat_path.exists())
            with strat_path.open() as handle:
                strat = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(strat)
            self.assertEqual(set(strat[0].keys()),
                             {"method", "modality", "ciwd_stratum", "gene_rows", "callable_rate", "overall_correct_call_rate"})
            # consensus methods are stratified alongside single tools
            self.assertIn("MajorityVote", {r["method"] for r in strat})
            # fixture uses common alleles -> stratum resolves to 'common'
            self.assertIn("common", {r["ciwd_stratum"] for r in strat})

            plaus_path = outdir / "tables" / "summary_ciwd_plausibility.tsv"
            self.assertTrue(plaus_path.exists())

            # harmonized rows carry the additive CIWD annotation columns
            with (outdir / "tables" / "harmonized_benchmark_rows.tsv").open() as handle:
                harmonized = list(csv.DictReader(handle, delimiter="\t"))
            for col in ("truth_ciwd_stratum", "call_implausible"):
                self.assertIn(col, harmonized[0])


if __name__ == "__main__":
    unittest.main()
