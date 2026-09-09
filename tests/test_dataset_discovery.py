import json
import tempfile
import unittest
from pathlib import Path

from champhla_confirmation.dataset_discovery import (
    CROSSWALK_FIELDS,
    DATASET_FIELDS,
    PILOT_FIELDS,
    TRUTH_FIELDS,
    _validate_crosswalk,
    _validate_pilot,
    _validate_truth,
    audit_dataset_discovery_registry,
)
from champhla_confirmation.dataset_sources import (
    RELATED_GROUPS,
    locus_has_exact_two_field_truth,
    normalize_sample_name,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INVENTORY = PROJECT_ROOT / "external" / "dataset_discovery"


class DatasetDiscoveryTests(unittest.TestCase):
    def test_checked_in_inventory_is_valid_and_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = audit_dataset_discovery_registry(
                str(INVENTORY / "dataset_registry.tsv"),
                str(INVENTORY / "truth_registry.tsv"),
                str(INVENTORY / "sample_crosswalk.tsv"),
                str(INVENTORY / "pilot_manifest.tsv"),
                str(Path(tmp) / "audit.json"),
            )
        self.assertEqual(result["dataset_rows"], 11)
        self.assertEqual(result["truth_rows"], 10)
        self.assertEqual(len(result["complete_ready_subjects"]["rnaseq"]), 11)
        self.assertEqual(result["pilot_subject_target"], 12)
        self.assertFalse(result["acceptance"]["public_rna_pilot_ready"])
        self.assertTrue(result["acceptance"]["overlap_report_complete"])
        self.assertTrue(result["acceptance"]["truth_circularity_report_complete"])
        self.assertFalse(result["acceptance"]["relatedness_report_complete"])
        self.assertFalse(result["acceptance"]["public_wes_pilot_ready"])
        self.assertFalse(result["acceptance"]["public_matched_wes_rna_pilot_ready"])
        self.assertFalse(result["acceptance"]["matched_wgs_rna_candidate_accession_resolved"])
        self.assertFalse(result["acceptance"]["wes_minimum_closed_by_verified_cohort"])
        self.assertFalse(result["passed"])
        self.assertTrue(result["registry_valid"])
        self.assertEqual("executable", result["lane_states"]["NCI60_RNA"]["state"])
        self.assertEqual(11, result["lane_states"]["NCI60_RNA"]["complete_executable_subjects"])

    def test_pilot_rejects_truth_bearing_columns(self):
        row = {field: "" for field in PILOT_FIELDS}
        row.update({
            "cohort": "TEST", "subject": "S1", "modality": "rnaseq", "gene": "A",
            "source": "ENA", "source_accession": "SRR1", "input_type": "fastq",
            "read1_url": "https://example.test/R1", "read2_url": "https://example.test/R2",
            "access_class": "public", "independence_stratum": "donor_independent",
            "gate_status": "ready", "truth_allele1": "A*01:01",
        })
        with self.assertRaisesRegex(ValueError, "forbidden truth-bearing columns"):
            _validate_pilot([row])

    def test_ready_paired_fastq_requires_read_two(self):
        row = {field: "" for field in PILOT_FIELDS}
        row.update({
            "cohort": "TEST", "subject": "S1", "modality": "rnaseq", "gene": "A",
            "source": "ENA", "source_accession": "SRR1", "input_type": "fastq",
            "read1_url": "https://example.test/R1", "access_class": "public",
            "independence_stratum": "donor_independent", "gate_status": "ready",
        })
        with self.assertRaisesRegex(ValueError, "lacks read2_url"):
            _validate_pilot([row])

    def test_registry_contracts_include_required_audit_dimensions(self):
        self.assertTrue({"assay", "read_layout", "reference_build", "access_class",
                         "access_terms", "tool_compatibility", "sample_link_status",
                         "capacity_met"}.issubset(DATASET_FIELDS))
        self.assertTrue({"method", "resolution", "ambiguity_status", "rna_informed",
                         "circularity", "accuracy_eligible"}.issubset(TRUTH_FIELDS))

    def test_alias_ambiguity_null_and_derivative_rules(self):
        self.assertEqual(normalize_sample_name("MDA-MB_435"), "MDAMB435")
        self.assertEqual(RELATED_GROUPS["NCI-ADR-RES"][1], "derivative_or_duplicate")
        self.assertTrue(locus_has_exact_two_field_truth("0101N, 0101N"))
        self.assertFalse(locus_has_exact_two_field_truth("02"))
        self.assertFalse(locus_has_exact_two_field_truth("0101a"))
        self.assertFalse(locus_has_exact_two_field_truth("N.R."))

    def test_blocked_and_related_crosswalk_rows_require_evidence(self):
        row = {field: "" for field in CROSSWALK_FIELDS}
        row.update({
            "cohort": "TEST", "repository_sample": "S1", "donor": "D1",
            "modality": "wes", "sequencing_accession": "SRR1",
            "access_class": "public", "independence_stratum": "donor_independent",
            "gate_status": "blocked_truth_ambiguity", "existing_champhla_overlap": "none",
        })
        with self.assertRaisesRegex(ValueError, "lacks exclusion_reason"):
            _validate_crosswalk([row])
        row["gate_status"] = "excluded_related"
        row["exclusion_reason"] = "duplicate cell-line derivative"
        with self.assertRaisesRegex(ValueError, "lacks relationship evidence"):
            _validate_crosswalk([row])
        row["relatedness"] = "same donor derivative"
        with self.assertRaisesRegex(ValueError, "duplicate crosswalk row"):
            _validate_crosswalk([row, dict(row)])

    def test_circular_or_disputed_truth_cannot_be_accuracy_eligible(self):
        row = {field: "" for field in TRUTH_FIELDS}
        row.update({
            "truth_id": "T1", "cohort": "TEST", "method": "RNA reconciliation",
            "source_material": "RNA", "loci": "A,B,C", "resolution": "two-field",
            "ambiguity_status": "disputed", "typing_date": "2026",
            "reconciliation": "RNA-informed", "rna_informed": "yes",
            "truth_independence": "outcome-informed", "circularity": "partly_circular",
            "accuracy_eligible": "yes", "reported_subjects": "1", "usable_truth_count": "1",
            "access_class": "public", "source_accession": "TEST1",
            "evidence_url": "https://example.test", "audit_status": "conditional",
        })
        with self.assertRaisesRegex(ValueError, "cannot be accuracy eligible"):
            _validate_truth([row])

    def test_partial_loci_remain_explicit_in_truth_registry(self):
        rows = (INVENTORY / "truth_registry.tsv").read_text(encoding="utf-8")
        self.assertIn("A,B,C,DRB1,DQB1", rows)
        self.assertIn("partial class II by subject", rows)

    def test_dataset_roles_are_preassigned_and_not_pooled(self):
        roles = json.loads((PROJECT_ROOT / "configs" / "dataset_roles.json").read_text())
        self.assertTrue(roles["frozen_before_corrected_predictions"])
        self.assertIn("never pooled", roles["pooling_rule"])
        identifiers = [row["dataset_id"] for row in roles["datasets"]]
        self.assertEqual(len(identifiers), len(set(identifiers)))
        rejected = next(row for row in roles["datasets"]
                        if row["dataset_id"] == "REJECTED_OR_CIRCULAR_TRUTH")
        self.assertFalse(rejected["accuracy_eligible"])


if __name__ == "__main__":
    unittest.main()
