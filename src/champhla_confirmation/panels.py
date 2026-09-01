from __future__ import annotations

GENES = ("A", "B", "C")
MODALITIES = ("wgs", "wes", "rnaseq")
PANELS = {
    "wgs": ("HLA-HD", "Kourami", "OptiType", "SpecHLA", "T1K"),
    "wes": ("HLA-HD", "OptiType", "POLYSOLVER", "SpecHLA", "T1K"),
    "rnaseq": ("ArcasHLA", "HLA-HD", "OptiType", "T1K"),
}

METHOD_BASELINE = "SimpleTwoThirdsConsensus"
METHOD_PLURALITY = "SimplePluralityLex"
METHOD_RAW_CC = "ChampionChallenger"
METHOD_GUARDED_CC = "TwoThirdsGuardedCC"

