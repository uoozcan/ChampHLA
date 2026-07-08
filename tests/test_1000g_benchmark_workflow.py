import csv
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
FIXTURES = REPO / "tests" / "fixtures"
BUILD_SCRIPT = REPO / "bin" / "build_1000g_benchmark_manifests.py"
PHASE_GATED_SCRIPT = REPO / "bin" / "build_1000g_phase_gated_inputs.py"
RUN_SCRIPT = REPO / "bin" / "run_1000g_benchmark.py"


class ThousandGenomesBenchmarkWorkflowTest(unittest.TestCase):
    def test_manifest_builder_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "manifests"
            subprocess.run(
                [
                    "python3",
                    str(BUILD_SCRIPT),
                    "--truth",
                    str(FIXTURES / "truth.tsv"),
                    "--sequencing",
                    str(FIXTURES / "sequencing_source.tsv"),
                    "--output-dir",
                    str(outdir),
                    "--acquisition-date",
                    "2014-07-25",
                    "--supported-loci",
                    "A,B,C",
                ],
                check=True,
                cwd=str(REPO),
            )
            truth_manifest = outdir / "truth_manifest.tsv"
            sequencing_manifest = outdir / "sequencing_manifest.tsv"
            cohort_manifest = outdir / "cohort_manifest.tsv"
            self.assertTrue(truth_manifest.exists())
            self.assertTrue(sequencing_manifest.exists())
            self.assertTrue(cohort_manifest.exists())
            with cohort_manifest.open("r", encoding="utf-8") as handle:
                cohort_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(cohort_rows), 2)
            include_map = {row["sample"]: row for row in cohort_rows}
            self.assertEqual(include_map["S1"]["include"], "1")
            self.assertEqual(include_map["S2"]["include"], "1")

    def test_phase_gated_input_builder_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "phase"
            wgs = Path(tmpdir) / "wgs_results"
            wes = Path(tmpdir) / "wes_results"
            rna = Path(tmpdir) / "rna_results"
            (wgs / "S1" / "spechla").mkdir(parents=True)
            (wes / "by_tool" / "optitype").mkdir(parents=True)
            (rna / "S1" / "arcashla").mkdir(parents=True)
            (wgs / "S1" / "spechla" / "S1_spechla.txt").write_text("A\tA*01:01\tA*02:01\n", encoding="utf-8")
            (wes / "by_tool" / "optitype" / "S1_optitype.txt").write_text("A\tA*01:01\tA*02:01\n", encoding="utf-8")
            (rna / "S1" / "arcashla" / "S1_arcashla.txt").write_text("A\tA*01:01\tA*02:01\n", encoding="utf-8")
            subprocess.run(
                [
                    "python3",
                    str(PHASE_GATED_SCRIPT),
                    "--truth-csv",
                    str(FIXTURES / "truth.tsv"),
                    "--output-dir",
                    str(outdir),
                    "--wgs-results",
                    str(wgs),
                    "--wes-results",
                    str(wes),
                    "--rnaseq-results",
                    str(rna),
                    "--samples",
                    "S1",
                    "--supported-loci",
                    "A,B,C",
                ],
                check=True,
                cwd=str(REPO),
            )
            self.assertTrue((outdir / "truth_long.tsv").exists())
            self.assertTrue((outdir / "sequencing_source.tsv").exists())
            self.assertTrue((outdir / "tool_availability_seed.tsv").exists())
            with (outdir / "tool_availability_seed.tsv").open("r", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["status"] == "not_available" for row in rows))

    def test_realdata_runner_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "benchmark"
            subprocess.run(
                [
                    "python3",
                    str(RUN_SCRIPT),
                    "--config",
                    str(FIXTURES / "benchmark_1000g_config.yaml"),
                    "--output-dir",
                    str(outdir),
                ],
                check=True,
                cwd=str(REPO),
            )
            tables_dir = outdir / "tables"
            self.assertTrue((tables_dir / "truth_manifest.tsv").exists())
            self.assertTrue((tables_dir / "sequencing_manifest.tsv").exists())
            self.assertTrue((tables_dir / "cohort_manifest.tsv").exists())
            self.assertTrue((tables_dir / "consensus_runtime_weights.json").exists())
            self.assertTrue((tables_dir / "tool_confidence_weights.tsv").exists())
            self.assertTrue((tables_dir / "summary_full_cohort.tsv").exists())
            self.assertTrue((tables_dir / "method_comparison.tsv").exists())
            self.assertTrue((tables_dir / "meta_method_comparison.tsv").exists())
            self.assertTrue((tables_dir / "meta_consensus_calls.tsv").exists())
            self.assertTrue((tables_dir / "tool_availability_by_sample.tsv").exists())
            self.assertTrue((tables_dir / "tool_availability_by_modality.tsv").exists())
            self.assertTrue((tables_dir / "tool_pairwise_agreement.tsv").exists())
            self.assertTrue((tables_dir / "consensus_decision_trace.tsv").exists())
            self.assertTrue((tables_dir / "summary_full_cohort_multiresolution.tsv").exists())
            metadata = json.loads((tables_dir / "benchmark_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["truth_acquisition_date"], "2014-07-25")
            self.assertEqual(metadata["final_tri_modal_cohort_size"], 2)
            self.assertEqual(metadata["supported_loci"], ["A", "B", "C"])
            self.assertEqual(metadata["per_modality_sample_counts"]["after_filtering"], {"wgs": 2, "wes": 2, "rnaseq": 2})
            self.assertEqual(metadata["tool_coverage_policy"], "phase_gated")
            self.assertIn("not_available_tools_by_modality", metadata)
            self.assertIn("available_tools_by_modality", metadata)
            with (tables_dir / "reference_metadata.tsv").open("r", encoding="utf-8") as handle:
                reference_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(reference_rows[0]["truth_acquisition_date"], "2014-07-25")
            self.assertEqual(reference_rows[0]["supported_loci"], "A,B,C")

    def test_realdata_runner_probabilistic_mode_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "benchmark"
            config_payload = yaml.safe_load((FIXTURES / "benchmark_1000g_config.yaml").read_text(encoding="utf-8"))
            config_payload.setdefault("benchmark", {})
            config_payload["benchmark"]["mode"] = "probabilistic_recalibrated"
            config_payload["benchmark"]["probabilistic_calibration"] = {"method": "platt", "cv_strategy": "loo"}
            config_path = tmpdir / "benchmark_prob.yaml"
            config_path.write_text(yaml.safe_dump(config_payload), encoding="utf-8")
            subprocess.run(
                [
                    "python3",
                    str(RUN_SCRIPT),
                    "--config",
                    str(config_path),
                    "--output-dir",
                    str(outdir),
                    "--benchmark-mode",
                    "probabilistic_recalibrated",
                ],
                check=True,
                cwd=str(REPO),
            )
            tables_dir = outdir / "tables"
            self.assertTrue((tables_dir / "cross_validation_weight_summary.tsv").exists())
            self.assertTrue((tables_dir / "confidence_calibration_summary.tsv").exists())
            self.assertTrue((tables_dir / "tool_population_diagnostics.tsv").exists())
            self.assertTrue((tables_dir / "tool_population_calibration.tsv").exists())
            self.assertTrue((tables_dir / "population_method_comparison.tsv").exists())
            self.assertTrue((tables_dir / "population_resolution_summary.tsv").exists())
            with (tables_dir / "tool_confidence_weights.tsv").open("r", encoding="utf-8") as handle:
                weight_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(all(row["benchmark_mode"] == "probabilistic_recalibrated" for row in weight_rows))
            with (tables_dir / "reference_metadata.tsv").open("r", encoding="utf-8") as handle:
                reference_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(reference_rows[0]["benchmark_mode"], "probabilistic_recalibrated")

    def test_realdata_runner_ensemble_ablation_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "benchmark"
            config_payload = yaml.safe_load((FIXTURES / "benchmark_1000g_config.yaml").read_text(encoding="utf-8"))
            config_payload.setdefault("benchmark", {})
            config_payload["benchmark"]["ensemble_ablation"] = {
                "enabled": True,
                "tool_subsets": [
                    {"name": "all_tools", "tools": ["SpecHLA", "OptiType", "ArcasHLA"]},
                    {"name": "optitype_only", "tools": ["OptiType"]},
                ],
            }
            config_path = tmpdir / "benchmark_ablation.yaml"
            config_path.write_text(yaml.safe_dump(config_payload), encoding="utf-8")
            subprocess.run(
                [
                    "python3",
                    str(RUN_SCRIPT),
                    "--config",
                    str(config_path),
                    "--output-dir",
                    str(outdir),
                ],
                check=True,
                cwd=str(REPO),
            )
            tables_dir = outdir / "tables"
            self.assertTrue((tables_dir / "ablation_method_comparison.tsv").exists())
            self.assertTrue((tables_dir / "ablation_method_per_gene.tsv").exists())
            self.assertTrue((tables_dir / "ablation_summary.json").exists())
            with (tables_dir / "ablation_method_comparison.tsv").open("r", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["ablation_name"] == "all_tools" and row["method"] == "WeightedConsensus" for row in rows))
            self.assertTrue(any(row["ablation_name"] == "optitype_only" and row["method"] == "MetaConsensus" for row in rows))

    def test_realdata_runner_locus_expert_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "benchmark"
            config_payload = yaml.safe_load((FIXTURES / "benchmark_1000g_config.yaml").read_text(encoding="utf-8"))
            config_payload.setdefault("benchmark", {})
            config_payload["benchmark"]["locus_expert_consensus"] = {
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
            config_path.write_text(yaml.safe_dump(config_payload), encoding="utf-8")
            subprocess.run(
                [
                    "python3",
                    str(RUN_SCRIPT),
                    "--config",
                    str(config_path),
                    "--output-dir",
                    str(outdir),
                ],
                check=True,
                cwd=str(REPO),
            )
            tables_dir = outdir / "tables"
            self.assertTrue((tables_dir / "locus_expert_consensus_calls.tsv").exists())
            self.assertTrue((tables_dir / "locus_expert_decision_trace.tsv").exists())
            self.assertTrue((tables_dir / "locus_expert_method_comparison.tsv").exists())
            self.assertTrue((tables_dir / "locus_expert_method_per_gene.tsv").exists())
            self.assertTrue((tables_dir / "locus_expert_summary.json").exists())
            with (tables_dir / "locus_expert_method_comparison.tsv").open("r", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["panel_name"] == "fixture_panel" and row["method"] == "LocusExpertConsensus" for row in rows))

    def test_realdata_runner_champion_challenger_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "benchmark"
            config_payload = yaml.safe_load((FIXTURES / "benchmark_1000g_config.yaml").read_text(encoding="utf-8"))
            config_payload.setdefault("benchmark", {})
            config_payload["benchmark"]["champion_challenger"] = {
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
            config_path.write_text(yaml.safe_dump(config_payload), encoding="utf-8")
            subprocess.run(
                [
                    "python3",
                    str(RUN_SCRIPT),
                    "--config",
                    str(config_path),
                    "--output-dir",
                    str(outdir),
                ],
                check=True,
                cwd=str(REPO),
            )
            tables_dir = outdir / "tables"
            self.assertTrue((tables_dir / "champion_challenger_calls.tsv").exists())
            self.assertTrue((tables_dir / "champion_challenger_trace.tsv").exists())
            self.assertTrue((tables_dir / "champion_challenger_method_comparison.tsv").exists())
            self.assertTrue((tables_dir / "champion_challenger_method_per_gene.tsv").exists())
            self.assertTrue((tables_dir / "champion_challenger_summary.json").exists())
            with (tables_dir / "champion_challenger_method_comparison.tsv").open("r", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["method"] == "ChampionChallenger" for row in rows))

    def test_realdata_runner_gated_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "benchmark"
            config_payload = yaml.safe_load((FIXTURES / "benchmark_1000g_config.yaml").read_text(encoding="utf-8"))
            config_payload.setdefault("benchmark", {})
            config_payload["benchmark"]["champion_challenger"] = {
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
            config_payload["benchmark"]["consensus_gating"] = {
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
            config_path.write_text(yaml.safe_dump(config_payload), encoding="utf-8")
            subprocess.run(
                [
                    "python3",
                    str(RUN_SCRIPT),
                    "--config",
                    str(config_path),
                    "--output-dir",
                    str(outdir),
                ],
                check=True,
                cwd=str(REPO),
            )
            tables_dir = outdir / "tables"
            self.assertTrue((tables_dir / "gated_consensus_calls.tsv").exists())
            self.assertTrue((tables_dir / "gated_consensus_trace.tsv").exists())
            self.assertTrue((tables_dir / "gated_method_comparison.tsv").exists())
            self.assertTrue((tables_dir / "gated_method_per_gene.tsv").exists())
            self.assertTrue((tables_dir / "gated_summary.json").exists())
            with (tables_dir / "gated_method_comparison.tsv").open("r", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["method"] == "GatedConsensus" for row in rows))

    def test_realdata_runner_threshold_sweep_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "benchmark"
            config_payload = yaml.safe_load((FIXTURES / "benchmark_1000g_config.yaml").read_text(encoding="utf-8"))
            config_payload["benchmark_analysis"] = {
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
            config_path = tmpdir / "benchmark_sweeps.yaml"
            config_path.write_text(yaml.safe_dump(config_payload), encoding="utf-8")
            subprocess.run(
                [
                    "python3",
                    str(RUN_SCRIPT),
                    "--config",
                    str(config_path),
                    "--output-dir",
                    str(outdir),
                ],
                check=True,
                cwd=str(REPO),
            )
            sweeps_dir = outdir / "sweeps"
            self.assertTrue((sweeps_dir / "wes_weighted_threshold_sweep.tsv").exists())
            self.assertTrue((sweeps_dir / "wes_champion_override_sweep.tsv").exists())
            self.assertTrue((sweeps_dir / "wes_sweep_summary.json").exists())

    def test_realdata_runner_threshold_sweep_outputs_rnaseq(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "benchmark"
            config_payload = yaml.safe_load((FIXTURES / "benchmark_1000g_config.yaml").read_text(encoding="utf-8"))
            config_payload["runs"] = [dict(run) for run in config_payload["runs"] if run.get("modality") == "rnaseq"]
            for run in config_payload["runs"]:
                if run.get("tool") == "SpecHLA":
                    run.pop("coverage_only", None)
            config_payload["benchmark_analysis"] = {
                "threshold_sweeps": {
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
            }
            config_path = tmpdir / "benchmark_sweeps_rna.yaml"
            config_path.write_text(yaml.safe_dump(config_payload), encoding="utf-8")
            subprocess.run(
                [
                    "python3",
                    str(RUN_SCRIPT),
                    "--config",
                    str(config_path),
                    "--output-dir",
                    str(outdir),
                ],
                check=True,
                cwd=str(REPO),
            )
            sweeps_dir = outdir / "sweeps"
            self.assertTrue((sweeps_dir / "rna_weighted_threshold_sweep.tsv").exists())
            self.assertTrue((sweeps_dir / "rna_champion_override_sweep.tsv").exists())
            self.assertTrue((sweeps_dir / "rna_sweep_summary.json").exists())


if __name__ == "__main__":
    unittest.main()
