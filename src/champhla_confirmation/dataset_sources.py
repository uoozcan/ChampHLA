from __future__ import annotations

import re


NCI_ALIAS = {
    "MDA-MB435": "MDA-MB-435",
    "M 14": "M14",
    "MOLT 4": "MOLT-4",
    "LOX IMVI": "LOX-IMVI",
    "A549/ATCC": "A549",
    "NCI/ADR-RES": "NCI-ADR-RES",
}
NCI_RNA_ALIAS = {"786_0": "786-O"}
RELATED_GROUPS = {
    "M14": ("M14_MDA-MB-435", "canonical"),
    "MDA-MB-435": ("M14_MDA-MB-435", "derivative_or_duplicate"),
    "SNB-19": ("SNB-19_U251", "canonical"),
    "U251": ("SNB-19_U251", "derivative_or_duplicate"),
    "OVCAR-8": ("OVCAR-8_NCI-ADR-RES", "canonical"),
    "NCI-ADR-RES": ("OVCAR-8_NCI-ADR-RES", "derivative_or_duplicate"),
}


def normalize_sample_name(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper().replace("ATCC", ""))


def normalize_nci_name(value: str) -> str:
    cleaned = " ".join(value.strip().split())
    return NCI_ALIAS.get(cleaned, cleaned.replace(" ", "-") if cleaned in {"MOLT 4"} else cleaned)


def locus_has_exact_two_field_truth(value: str) -> bool:
    if not value or "N.R." in value or "NEW" in value.upper():
        return False
    for token in (part.strip() for part in value.split(",")):
        if token.lower().endswith("a"):
            return False
        if len(re.sub(r"[^0-9]", "", token)) < 4:
            return False
    return True


# Read-only compatibility names for the original discovery builder.
_norm_name = normalize_sample_name
_nci_name = normalize_nci_name
_locus_exact = locus_has_exact_two_field_truth
