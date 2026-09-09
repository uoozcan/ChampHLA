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
METHOD_MV_FLOOR = "MVFlooredCC"

PRIMARY_METHOD = METHOD_PLURALITY
PLURALITY_METHOD_VERSION = "pair-plurality-v2"
METHOD_ALIASES = {
    "MajorityVote": METHOD_PLURALITY,
    "Majority Vote": METHOD_PLURALITY,
    "PairLevelPluralityConsensus": METHOD_PLURALITY,
    "Pair-level plurality consensus": METHOD_PLURALITY,
}


def canonical_method(value: str) -> str:
    """Return the stable serialized method identifier for accepted aliases."""
    token = value.strip()
    return METHOD_ALIASES.get(token, token)
