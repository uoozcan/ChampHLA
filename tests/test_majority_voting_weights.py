"""Weight lookup for the consensus vote.

The consensus step returned NO_CALL for every gene on a sample that five callers
had typed correctly, with `chosen_by=no_nonzero_weight`. Two reasons, both here:

  * conf/tool_weights_{wgs,wes,rna}.json are flat {tool: number} maps, and
    lookup_weight() understood only the nested runtime schema and a legacy
    raw_accuracy one -- everything else fell through to 0.0.
  * those files spell tools "OptiType" and "HLA-HD" while the rows in
    aggregated_calls.tsv say "optitype" and "hlahd", because that is what
    modules/aggregation.nf emits.

These pin both, so a weight file that ships with the pipeline keeps working with
the code that reads it.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

_spec = importlib.util.spec_from_file_location("mv", REPO / "bin" / "majority_voting.py")
mv = importlib.util.module_from_spec(_spec)
sys.modules["mv"] = mv
_spec.loader.exec_module(mv)


def row(tool, modality="wgs", gene="A"):
    return {"tool": tool, "modality": modality, "gene": gene}


class FlatWeightFileTests(unittest.TestCase):
    FLAT = {"ArcasHLA": 0.0809, "HLA-HD": 0.2686, "OptiType": 0.4974, "T1K": 0.3366}

    def test_a_flat_file_resolves(self):
        self.assertAlmostEqual(mv.lookup_weight(self.FLAT, row("OptiType"), False), 0.4974)

    def test_the_spelling_in_aggregated_calls_resolves(self):
        """aggregation.nf emits lower-case, separator-free tokens."""
        self.assertAlmostEqual(mv.lookup_weight(self.FLAT, row("optitype"), False), 0.4974)
        self.assertAlmostEqual(mv.lookup_weight(self.FLAT, row("hlahd"), False), 0.2686)
        self.assertAlmostEqual(mv.lookup_weight(self.FLAT, row("t1k"), False), 0.3366)

    def test_hla_hd_however_it_is_punctuated(self):
        for spelling in ("HLA-HD", "hlahd", "hla_hd", "Hla-Hd"):
            self.assertAlmostEqual(mv.lookup_weight(self.FLAT, row(spelling), False), 0.2686,
                                   msg=spelling)

    def test_a_tool_absent_from_the_file_scores_zero(self):
        """POLYSOLVER is absent from conf/tool_weights_wgs.json. Zero, not a crash."""
        self.assertEqual(mv.lookup_weight(self.FLAT, row("polysolver"), False), 0.0)

    def test_weights_are_clipped_to_the_unit_interval(self):
        self.assertEqual(mv.lookup_weight({"X": 5.0}, row("X"), False), 1.0)
        self.assertEqual(mv.lookup_weight({"X": -2.0}, row("X"), False), 0.0)


class ShippedWeightFileTests(unittest.TestCase):
    """The files in conf/ must work with this code, not merely parse."""

    def test_every_shipped_file_gives_a_nonzero_weight_to_its_tools(self):
        for modality, name in (("wgs", "tool_weights_wgs.json"),
                               ("wes", "tool_weights_wes.json"),
                               ("rnaseq", "tool_weights_rna.json")):
            path = REPO / "conf" / name
            weights = json.loads(path.read_text(encoding="utf-8"))
            for tool in weights:
                value = mv.lookup_weight(weights, row(tool, modality), False)
                self.assertGreater(value, 0.0, f"{name}: {tool} resolved to zero")

    def test_the_wes_panel_all_resolves(self):
        """The paper's WES panel; every member must carry weight."""
        weights = json.loads((REPO / "conf" / "tool_weights_wes.json").read_text(encoding="utf-8"))
        for tool in ("hlahd", "optitype", "polysolver", "spechla", "t1k"):
            self.assertGreater(mv.lookup_weight(weights, row(tool, "wes"), False), 0.0, tool)

    def test_the_rna_panel_all_resolves(self):
        weights = json.loads((REPO / "conf" / "tool_weights_rna.json").read_text(encoding="utf-8"))
        for tool in ("arcashla", "hlahd", "optitype", "t1k"):
            self.assertGreater(mv.lookup_weight(weights, row(tool, "rnaseq"), False), 0.0, tool)


class OtherSchemasStillWorkTests(unittest.TestCase):
    def test_the_nested_runtime_schema(self):
        nested = {"tool_weights": {"OptiType": {"wgs": {"final_weight": 0.9}}}}
        self.assertAlmostEqual(mv.lookup_weight(nested, row("OptiType"), False), 0.9)

    def test_the_nested_schema_also_tolerates_the_spelling(self):
        nested = {"tool_weights": {"HLA-HD": {"wgs": {"final_weight": 0.7}}}}
        self.assertAlmostEqual(mv.lookup_weight(nested, row("hlahd"), False), 0.7)

    def test_equal_mode_ignores_the_file_entirely(self):
        """Equal-weight plurality is the project's primary method."""
        self.assertEqual(mv.lookup_weight({}, row("anything"), False, equal_mode=True), 1.0)


if __name__ == "__main__":
    unittest.main()
