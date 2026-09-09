from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree

from .io import read_tsv, sha256, write_json, write_tsv
from .registry import validate_registry


PLACEHOLDER_RE = re.compile(
    r"\[(?:AUTHOR|TODO|TBD|PLACEHOLDER)[^\]]*\]|\b(?:TODO|TBD)\b|author verification",
    re.I,
)
INVALID_WGS_RE = re.compile(r"(?:0\.399|0\.494|\+30\.66|79/411|205/411)")
RESULT_REF_RE = re.compile(r"\[RESULT:([A-Za-z0-9_.:-]+)\]")
OVERCLAIM_RE = re.compile(
    r"\b(?:best (?:tool|method) (?:for|across) all|universally best|all NGS types|"
    r"outperforms every (?:tool|caller|method)|equivalent to the best)\b",
    re.I,
)
NUMERIC_RESULT_RE = re.compile(
    r"(?:\b\d+\s*(?:/|of)\s*\d+\b|\b\d+(?:\.\d+)?\s*%|"
    r"\b0\.\d{3,}\b|[+−-]\d+(?:\.\d+)?\s*(?:percentage\s+)?points?\b)",
    re.I,
)


def extract_docx_text(path: str) -> str:
    with zipfile.ZipFile(path) as archive:
        root = ElementTree.fromstring(archive.read("word/document.xml"))
    paragraphs = []
    namespace = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    for paragraph in root.iter(namespace + "p"):
        text = "".join(node.text or "" for node in paragraph.iter(namespace + "t"))
        if text.strip():
            paragraphs.append(text)
    return "\n".join(paragraphs) + "\n"


def write_docx_text(source: str, output: str) -> None:
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(extract_docx_text(source), encoding="utf-8")


def audit_claims(manuscript: str, registry_path: str, claims_path: str, output: str,
                 root: str | None = None, supplement: str | None = None) -> dict:
    path = Path(manuscript)
    text = extract_docx_text(str(path)) if path.suffix.lower() == ".docx" else path.read_text(encoding="utf-8")
    supplement_path = Path(supplement) if supplement else None
    supplement_text = ""
    if supplement_path:
        supplement_text = (extract_docx_text(str(supplement_path))
                           if supplement_path.suffix.lower() == ".docx"
                           else supplement_path.read_text(encoding="utf-8"))
    audited_text = text + ("\n" + supplement_text if supplement_text else "")
    rows = read_tsv(registry_path)
    registry = {row["result_id"]: row for row in rows}
    claims = read_tsv(claims_path)
    failures = validate_registry(registry_path, root)
    warnings = []

    placeholders = sorted(set(PLACEHOLDER_RE.findall(audited_text)))
    if placeholders:
        failures.append(f"manuscript contains {len(placeholders)} placeholder marker(s)")
    if INVALID_WGS_RE.search(text):
        failures.append("manuscript contains a historical invalid-WGS performance number")
    if OVERCLAIM_RE.search(audited_text):
        failures.append("manuscript contains an unsupported universal or equivalence claim")
    main_refs = RESULT_REF_RE.findall(text)
    supplement_refs = RESULT_REF_RE.findall(supplement_text)
    diagnostic_invalid_refs = []
    for result_id in main_refs:
        if result_id not in registry:
            failures.append(f"unknown result reference: {result_id}")
        elif registry[result_id]["validity"] != "valid":
            failures.append(f"manuscript cites invalid result: {result_id}")
    for paragraph in re.split(r"\n\s*\n", supplement_text):
        for result_id in RESULT_REF_RE.findall(paragraph):
            if result_id not in registry:
                failures.append(f"unknown supplement result reference: {result_id}")
            elif registry[result_id]["validity"] != "valid":
                if not re.search(r"\b(?:invalid|diagnostic|excluded|defect)\b", paragraph, re.I):
                    failures.append(
                        f"supplement cites invalid result without diagnostic labeling: {result_id}"
                    )
                else:
                    diagnostic_invalid_refs.append(result_id)
    for label, content in (("main", text), ("supplement", supplement_text)):
        for paragraph in re.split(r"\n\s*\n", content):
            if NUMERIC_RESULT_RE.search(paragraph) and not RESULT_REF_RE.search(paragraph):
                failures.append(f"{label} numerical result lacks a registry reference")
    external_valid = any(
        row.get("evidence_role") == "independent_validation" and row["validity"] == "valid"
        for row in rows
    )
    if re.search(r"\b(?:externally validated|independent external validation (?:confirmed|demonstrated))\b", audited_text, re.I) and not external_valid:
        failures.append("external-validation language is unsupported by the result registry")
    abstract_match = re.search(
        r"(?ims)^## Abstract\s*$\n(.*?)(?=^##\s+|\Z)", text,
    )
    abstract = abstract_match.group(1) if abstract_match else ""
    for result_id in RESULT_REF_RE.findall(abstract):
        if result_id in registry and registry[result_id].get("abstract_allowed") != "1":
            failures.append(f"abstract cites a result not allowed in the abstract: {result_id}")

    route = "method" if "method_conditional" in str(path) else "benchmark" if "benchmark" in str(path) else "source_original"
    active_claims = [claim for claim in claims if claim.get("draft") in {route, "shared"}]
    claim_fields = set(claims[0]) if claims else set()
    required_claim_fields = {"claim_id", "draft", "section", "claim", "status", "evidence_result_id", "action"}
    if not claims:
        failures.append("claim audit is empty")
    elif required_claim_fields - claim_fields:
        failures.append(f"claim audit missing fields: {sorted(required_claim_fields - claim_fields)}")
    else:
        for claim in active_claims:
            result_id = claim["evidence_result_id"]
            if result_id and result_id not in registry:
                failures.append(f"claim {claim['claim_id']} references unknown result {result_id}")
            if claim["status"] in {"invalid", "remove"} and claim["draft"] in {"benchmark", "method"}:
                failures.append(f"active draft contains rejected claim record {claim['claim_id']}")
            if claim["status"] == "pending":
                warnings.append(f"pending claim: {claim['claim_id']}")

    result = {
        "schema_version": "manuscript-claim-audit-1",
        "manuscript": _portable_path(path, root),
        "manuscript_sha256": sha256(path),
        "word_count": len(text.split()),
        "supplement": _portable_path(supplement_path, root) if supplement_path else "",
        "supplement_sha256": sha256(supplement_path) if supplement_path else "",
        "supplement_word_count": len(supplement_text.split()),
        "draft_route": route,
        "result_references": sorted(set(RESULT_REF_RE.findall(audited_text))),
        "diagnostic_invalid_result_references": sorted(set(diagnostic_invalid_refs)),
        "failures": sorted(set(failures)),
        "warnings": sorted(set(warnings)),
        "submission_ready": not failures and not warnings,
    }
    write_json(output, result)
    return result


def _portable_path(path: Path, root: str | None) -> str:
    resolved = path.resolve()
    if root:
        try:
            return resolved.relative_to(Path(root).resolve()).as_posix()
        except ValueError:
            pass
    return resolved.as_posix()


def write_source_issue_register(source_docx: str, output_tsv: str, output_json: str) -> dict:
    text = extract_docx_text(source_docx)
    issues = []
    def add(category: str, severity: str, evidence: str, action: str) -> None:
        issues.append({
            "issue_id": f"MS-{len(issues)+1:03d}", "category": category,
            "severity": severity, "evidence": evidence, "required_action": action,
        })
    if len(text.split()) > 7000:
        add("structure", "high", f"Draft contains {len(text.split())} words.",
            "Replace with journal-length main text and move tutorials/exploratory work out of the main narrative.")
    if PLACEHOLDER_RE.search(text):
        add("placeholder", "high", "Author-confirmation or verification markers remain.",
            "Resolve every marker before submission.")
    if INVALID_WGS_RE.search(text):
        add("invalid_result", "critical", "Abstract/body contains performance from the invalid historical WGS slice.",
            "Remove from performance claims; retain only as a QC lesson.")
    if "MV-floor" in text:
        add("obsolete_headline", "critical", "MV-floor is presented as the reported operating method.",
            "Use the benchmark route by default; keep MV-floor as a comparator/ablation.")
    if "implemented but unbenchmarked long-read" in text or "longreads_hifi" in text:
        add("scope", "medium", "Unvalidated long-read support occupies the main narrative.",
            "Move long-read implementation to software documentation or limitations.")
    if "FIMM" in text and "survival" in text.lower():
        add("exploratory_scope", "high", "FIMM, LOH/HED, and survival analyses are embedded in the draft.",
            "Remove from the main manuscripts and retain only in a clearly exploratory archive.")
    add("claim_consistency", "critical", "Multiple source reports recommend incompatible headlines.",
        "Generate all result prose and tables from the authoritative result registry.")
    add("validation", "critical", "Qualifying three-modality confirmation is incomplete.",
        "Keep method manuscript conditional until frozen unseen-subject and HPRC gates pass.")
    add("comparison", "high", "The selected abstaining baseline can exaggerate fixed-denominator gain if callability is hidden.",
        "Report callability/called-only accuracy and always-call comparators beside it.")
    write_tsv(output_tsv, issues, ["issue_id", "category", "severity", "evidence", "required_action"])
    result = {
        "schema_version": "source-manuscript-issue-register-1",
        "source_sha256": sha256(source_docx), "word_count": len(text.split()),
        "issues": len(issues), "critical": sum(row["severity"] == "critical" for row in issues),
    }
    write_json(output_json, result)
    return result
