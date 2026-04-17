import csv
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

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
            self.assertEqual(sorted({row["split"] for row in cohort_rows}), ["holdout", "training"])

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
            self.assertTrue((tables_dir / "tool_availability_by_sample.tsv").exists())
            self.assertTrue((tables_dir / "tool_availability_by_modality.tsv").exists())
            metadata = json.loads((tables_dir / "benchmark_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["truth_acquisition_date"], "2014-07-25")
            self.assertEqual(metadata["final_tri_modal_cohort_size"], 2)
            self.assertEqual(metadata["split_membership_summary"], {"holdout": 1, "training": 1})
            self.assertEqual(metadata["supported_loci"], ["A", "B", "C"])
            self.assertEqual(metadata["per_modality_sample_counts"]["after_filtering"], {"wgs": 2, "wes": 2, "rnaseq": 2})
            self.assertEqual(metadata["tool_coverage_policy"], "phase_gated")
            self.assertTrue(metadata["small_cohort_mode"])
            self.assertEqual(metadata["effective_training_samples"], ["S1"])
            self.assertEqual(metadata["effective_validation_samples"], [])
            self.assertEqual(metadata["effective_holdout_samples"], ["S2"])
            self.assertIn("not_available_tools_by_modality", metadata)
            self.assertIn("available_tools_by_modality", metadata)
            with (tables_dir / "reference_metadata.tsv").open("r", encoding="utf-8") as handle:
                reference_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(reference_rows[0]["truth_acquisition_date"], "2014-07-25")
            self.assertEqual(reference_rows[0]["supported_loci"], "A,B,C")


if __name__ == "__main__":
    unittest.main()
