import hashlib
import json
import zipfile
from pathlib import Path

from champhla_recovery.submission import (
    audit_bibliography,
    audit_declarations,
    build_submission_package,
    render_docx,
    validate_docx_parity,
)


ROOT = Path(__file__).resolve().parents[1]


def _build(tmp_path: Path):
    return build_submission_package(
        ROOT,
        "manuscripts/benchmark/manuscript.md",
        "manuscripts/benchmark/supplementary.md",
        "manuscripts/submission/references.bib",
        "manuscripts/submission/declarations_template.md",
        tmp_path / "submission",
        tmp_path / "submission_manifest.tsv",
        tmp_path / "parity.json",
        tmp_path / "bibliography.json",
        tmp_path / "declarations.json",
    )


def test_build_submission_has_semantic_parity_and_embeds_figures(tmp_path: Path):
    result = _build(tmp_path)
    assert result["parity"]["passed"] is True
    assert result["parity"]["figure_order_valid"] is True
    docx = result["files"]["main_docx"]
    with zipfile.ZipFile(docx) as archive:
        media = [name for name in archive.namelist() if name.startswith("word/media/")]
        assert len(media) == 2
        assert "word/comments.xml" not in archive.namelist()


def test_docx_generation_is_byte_deterministic(tmp_path: Path):
    one = tmp_path / "one.docx"
    two = tmp_path / "two.docx"
    source = ROOT / "manuscripts/benchmark/manuscript.md"
    render_docx(source, one)
    render_docx(source, two)
    assert hashlib.sha256(one.read_bytes()).digest() == hashlib.sha256(two.read_bytes()).digest()


def test_parity_detects_document_text_tampering(tmp_path: Path):
    source = tmp_path / "source.md"
    source.write_text("# Title\n\nOriginal text.\n", encoding="utf-8")
    docx = tmp_path / "source.docx"
    render_docx(source, docx)
    replacement = tmp_path / "tampered.docx"
    with zipfile.ZipFile(docx) as old, zipfile.ZipFile(replacement, "w") as new:
        for info in old.infolist():
            data = old.read(info.filename)
            if info.filename == "word/document.xml":
                data = data.replace(b"Original text.", b"Altered text.")
            new.writestr(info, data)
    result = validate_docx_parity(source, replacement)
    assert result["text_parity"] is False


def test_bibliography_resolves_all_manuscript_citations():
    result = audit_bibliography(
        [ROOT / "manuscripts/benchmark/manuscript.md",
         ROOT / "manuscripts/benchmark/supplementary.md"],
        ROOT / "manuscripts/submission/references.bib",
    )
    assert result["passed"] is True
    assert result["missing_keys"] == []


def test_author_declarations_are_an_explicit_human_gate():
    result = audit_declarations(ROOT / "manuscripts/submission/declarations_template.md")
    assert result["passed"] is False
    assert result["unresolved_author_markers"] == 6


def test_submission_audit_json_is_machine_readable(tmp_path: Path):
    _build(tmp_path)
    audit = json.loads((tmp_path / "parity.json").read_text(encoding="utf-8"))
    assert audit["schema_version"] == "champhla-submission-parity-1"
