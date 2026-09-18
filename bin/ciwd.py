#!/usr/bin/env python3
"""CIWD 3.0.0 allele-commonness lookup for the champHLA benchmark.

CIWD is an *allele-classification* catalogue (common / intermediate / well-documented /
not-CIWD across seven population groups), **not** a ground-truth genotype dataset — it carries no
per-sample calls and no reads. We use it only to stratify benchmark concordance by allele commonness
and to flag biologically implausible calls. See assets/README_CIWD.md.

Standalone by design (no import of hla_benchmark, to avoid a circular import): it applies its own
two-field normalisation, matching hla_benchmark.normalize_allele semantics — strip `HLA-`, keep
`GENE*field1:field2`, drop the expression suffix (N/Q/L/S/A/C). The catalogue is indexed under both
the published allele (e.g. `A*01:04N`) and its two-field key (`A*01:04`), so callers can pass either.
"""
import re
from pathlib import Path

ASSET_DIR = Path(__file__).resolve().parent.parent / "assets"
DEFAULT_TABLE = ASSET_DIR / "ciwd_3.0.0.tsv"
GGROUP_TABLE = ASSET_DIR / "ciwd_3.0.0_ggroup.tsv"

POPULATION_GROUPS = ["AFA", "API", "EURO", "MENA", "HIS", "NAM", "UNK"]

# Source category codes -> canonical labels used throughout the benchmark.
CATEGORY_MAP = {
    "C": "common",
    "I": "intermediate",
    "WD": "well_documented",
    "O": "not_ciwd",
    "o": "not_ciwd",
    "": "not_ciwd",
}
# Rarity ordering (common = easiest ... unknown = novel/never observed).
CATEGORY_ORDER = ["common", "intermediate", "well_documented", "not_ciwd", "unknown"]
_RARITY = {name: i for i, name in enumerate(CATEGORY_ORDER)}

_EXPRESSION_SUFFIX = re.compile(r"[NQLSAC]$")


def normalize_two_field(raw):
    """`HLA-A*01:04:02N` / `A*01:04N` -> `A*01:04`; returns '' if unparseable."""
    if raw is None:
        return ""
    value = str(raw).strip().replace("HLA-", "").replace("HLA_", "").replace(" ", "")
    if "*" not in value:
        return ""
    gene, fields = value.split("*", 1)
    fields = _EXPRESSION_SUFFIX.sub("", fields)
    parts = [p for p in fields.split(":") if p]
    if not gene or not parts:
        return ""
    return "%s*%s" % (gene.upper(), ":".join(parts[:2]))


class CiwdCatalogue:
    """Loaded CIWD table; maps an allele to its commonness category per population group."""

    def __init__(self, rows):
        # index: normalized-2field allele -> {group_or_'total': category_label}
        self._index = {}
        self._raw_index = {}
        for row in rows:
            allele = row.get("allele", "").strip()
            if not allele:
                continue
            cats = {"total": CATEGORY_MAP.get(row.get("ciwd_total", "").strip(), "not_ciwd")}
            for group in POPULATION_GROUPS:
                cats[group.lower()] = CATEGORY_MAP.get(row.get("ciwd_%s" % group, "").strip(), "not_ciwd")
            cats["cwd2"] = row.get("cwd2_category", "").strip()
            self._raw_index[allele] = cats
            key = normalize_two_field(allele)
            if not key:
                continue
            # prefer an entry whose total is an actual CIWD category over a not_ciwd/null collision
            existing = self._index.get(key)
            if existing is None or (existing["total"] == "not_ciwd" and cats["total"] != "not_ciwd"):
                self._index[key] = cats

    def __len__(self):
        return len(self._index)

    def category(self, allele, group="total"):
        """Return one of CATEGORY_ORDER for `allele`; `unknown` if absent from the catalogue."""
        cats = self._raw_index.get(str(allele).strip()) or self._index.get(normalize_two_field(allele))
        if cats is None:
            return "unknown"
        return cats.get(group.lower() if group != "total" else "total", "unknown")

    def cwd2_category(self, allele):
        cats = self._raw_index.get(str(allele).strip()) or self._index.get(normalize_two_field(allele))
        return cats.get("cwd2", "") if cats else ""

    def genotype_stratum(self, allele1, allele2, group="total"):
        """Commonness of a genotype = the rarer of its two alleles (drives per-genotype stratifying)."""
        cats = [self.category(a, group) for a in (allele1, allele2) if a and str(a).strip()]
        if not cats:
            return "unknown"
        return max(cats, key=lambda c: _RARITY.get(c, len(CATEGORY_ORDER)))

    def is_implausible(self, allele, group="total"):
        """True for a call that is not-CIWD or absent — a candidate error / novel allele (QC flag)."""
        return self.category(allele, group) in ("not_ciwd", "unknown")


def load_ciwd(path=DEFAULT_TABLE):
    """Load a CIWD TSV (as written by assets/build_ciwd_tsv.py). Missing file -> empty catalogue."""
    path = Path(path)
    rows = []
    if path.exists():
        with path.open("r", encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return CiwdCatalogue(rows)


if __name__ == "__main__":
    import sys
    cat = load_ciwd()
    print("loaded %d two-field alleles from %s" % (len(cat), DEFAULT_TABLE))
    for a in sys.argv[1:]:
        print("%-14s total=%-16s EURO=%-16s cwd2=%s"
              % (a, cat.category(a), cat.category(a, "EURO"), cat.cwd2_category(a) or "-"))
