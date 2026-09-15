from __future__ import annotations

import csv
import hashlib
import re
import shutil
import struct
import zipfile
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET

from .io import write_json


W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
PIC = "http://schemas.openxmlformats.org/drawingml/2006/picture"
PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"
NS = {"w": W, "r": R, "wp": WP, "a": A, "pic": PIC}
for prefix, uri in NS.items():
    ET.register_namespace(prefix, uri)


@dataclass(frozen=True)
class Block:
    kind: str
    text: str = ""
    level: int = 0
    rows: tuple[tuple[str, ...], ...] = ()
    image: str = ""


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _plain(text: str) -> str:
    text = re.sub(r"`([^`]*)`", r"\1", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"\1", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"\1", text)
    return re.sub(r"\s+", " ", text).strip()


def _table_cells(line: str) -> tuple[str, ...]:
    return tuple(_plain(cell) for cell in line.strip().strip("|").split("|"))


def parse_markdown(path: str | Path) -> list[Block]:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    blocks: list[Block] = []
    paragraph: list[str] = []
    index = 0

    def flush() -> None:
        if paragraph:
            blocks.append(Block("paragraph", _plain(" ".join(paragraph))))
            paragraph.clear()

    while index < len(lines):
        raw = lines[index].rstrip()
        stripped = raw.strip()
        if not stripped:
            flush()
            index += 1
            continue
        image = re.fullmatch(r"!\[(.+)]\(([^)]+)\)", stripped)
        heading = re.match(r"^(#{1,6})\s+(.+)$", stripped)
        bullet = re.match(r"^[-*]\s+(.+)$", stripped)
        if image:
            flush()
            blocks.append(Block("image", _plain(image.group(1)), image=image.group(2)))
            index += 1
        elif heading:
            flush()
            blocks.append(Block("heading", _plain(heading.group(2)), level=len(heading.group(1))))
            index += 1
        elif bullet:
            flush()
            blocks.append(Block("bullet", _plain(bullet.group(1))))
            index += 1
        elif stripped.startswith("|") and stripped.endswith("|"):
            flush()
            table_lines: list[str] = []
            while index < len(lines):
                candidate = lines[index].strip()
                if not (candidate.startswith("|") and candidate.endswith("|")):
                    break
                table_lines.append(candidate)
                index += 1
            rows = [_table_cells(line) for line in table_lines]
            if len(rows) > 1 and all(re.fullmatch(r":?-{3,}:?", cell.replace(" ", ""))
                                     for cell in rows[1]):
                del rows[1]
            blocks.append(Block("table", rows=tuple(rows)))
        else:
            paragraph.append(stripped)
            index += 1
    flush()
    return blocks


def _qn(namespace: str, name: str) -> str:
    return f"{{{namespace}}}{name}"


def _paragraph(text: str, style: str | None = None) -> ET.Element:
    p = ET.Element(_qn(W, "p"))
    if style:
        ppr = ET.SubElement(p, _qn(W, "pPr"))
        ET.SubElement(ppr, _qn(W, "pStyle"), {_qn(W, "val"): style})
    run = ET.SubElement(p, _qn(W, "r"))
    node = ET.SubElement(run, _qn(W, "t"))
    if text.startswith(" ") or text.endswith(" "):
        node.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    node.text = text
    return p


def _table(rows: tuple[tuple[str, ...], ...]) -> ET.Element:
    table = ET.Element(_qn(W, "tbl"))
    props = ET.SubElement(table, _qn(W, "tblPr"))
    borders = ET.SubElement(props, _qn(W, "tblBorders"))
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        ET.SubElement(borders, _qn(W, edge), {
            _qn(W, "val"): "single", _qn(W, "sz"): "4", _qn(W, "color"): "B7B7B7",
        })
    for row_index, row in enumerate(rows):
        tr = ET.SubElement(table, _qn(W, "tr"))
        for value in row:
            tc = ET.SubElement(tr, _qn(W, "tc"))
            p = _paragraph(value)
            if row_index == 0:
                run_props = ET.SubElement(p.find(_qn(W, "r")), _qn(W, "rPr"))
                ET.SubElement(run_props, _qn(W, "b"))
            tc.append(p)
    return table


def _png_dimensions(path: Path) -> tuple[int, int]:
    data = path.read_bytes()[:24]
    if len(data) < 24 or data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"DOCX renderer accepts PNG figures only: {path}")
    return struct.unpack(">II", data[16:24])


def _drawing(rel_id: str, alt: str, image_id: int, width: int, height: int) -> ET.Element:
    max_width = 6.2 * 914400
    cx = min(max_width, width / 300 * 914400)
    cy = cx * height / width
    xml = f"""<w:p xmlns:w="{W}" xmlns:r="{R}" xmlns:wp="{WP}"
      xmlns:a="{A}" xmlns:pic="{PIC}"><w:r><w:drawing><wp:inline>
      <wp:extent cx="{int(cx)}" cy="{int(cy)}"/><wp:docPr id="{image_id}" name="Figure {image_id}" descr="{_xml_attr(alt)}"/>
      <a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">
      <pic:pic><pic:nvPicPr><pic:cNvPr id="{image_id}" name="image{image_id}.png"/><pic:cNvPicPr/></pic:nvPicPr>
      <pic:blipFill><a:blip r:embed="{rel_id}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>
      <pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{int(cx)}" cy="{int(cy)}"/></a:xfrm>
      <a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic>
      </a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>"""
    return ET.fromstring(xml)


def _xml_attr(text: str) -> str:
    return (text.replace("&", "&amp;").replace('"', "&quot;")
            .replace("<", "&lt;").replace(">", "&gt;"))


def _xml_bytes(element: ET.Element) -> bytes:
    return b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' + ET.tostring(
        element, encoding="utf-8",
    )


def _styles() -> bytes:
    root = ET.Element(_qn(W, "styles"))
    for style_id, name, size, bold in (
        ("Normal", "Normal", "22", False), ("Title", "Title", "32", True),
        ("Heading1", "heading 1", "28", True), ("Heading2", "heading 2", "25", True),
        ("Heading3", "heading 3", "23", True), ("ListBullet", "List Bullet", "22", False),
        ("Caption", "Caption", "20", False),
    ):
        style = ET.SubElement(root, _qn(W, "style"), {
            _qn(W, "type"): "paragraph", _qn(W, "styleId"): style_id,
        })
        ET.SubElement(style, _qn(W, "name"), {_qn(W, "val"): name})
        rpr = ET.SubElement(style, _qn(W, "rPr"))
        ET.SubElement(rpr, _qn(W, "sz"), {_qn(W, "val"): size})
        if bold:
            ET.SubElement(rpr, _qn(W, "b"))
    return _xml_bytes(root)


def render_docx(markdown_path: str | Path, output_path: str | Path) -> dict:
    source = Path(markdown_path).resolve()
    output = Path(output_path)
    blocks = parse_markdown(source)
    document = ET.Element(_qn(W, "document"))
    body = ET.SubElement(document, _qn(W, "body"))
    relationships: list[tuple[str, str]] = []
    media: list[tuple[str, bytes]] = []
    figure_order: list[str] = []
    for block in blocks:
        if block.kind == "heading":
            style = "Title" if block.level == 1 else f"Heading{min(block.level - 1, 3)}"
            body.append(_paragraph(block.text, style))
        elif block.kind == "paragraph":
            body.append(_paragraph(block.text))
        elif block.kind == "bullet":
            body.append(_paragraph(f"• {block.text}", "ListBullet"))
        elif block.kind == "table":
            body.append(_table(block.rows))
        elif block.kind == "image":
            image_path = (source.parent / block.image).resolve()
            if not image_path.is_file():
                raise FileNotFoundError(f"figure missing: {image_path}")
            width, height = _png_dimensions(image_path)
            image_id = len(media) + 1
            rel_id = f"rId{image_id}"
            media_name = f"image{image_id}.png"
            relationships.append((rel_id, f"media/{media_name}"))
            media.append((media_name, image_path.read_bytes()))
            body.append(_drawing(rel_id, block.text, image_id, width, height))
            body.append(_paragraph(block.text, "Caption"))
            figure_order.append(block.text)
    sect = ET.SubElement(body, _qn(W, "sectPr"))
    ET.SubElement(sect, _qn(W, "pgSz"), {_qn(W, "w"): "12240", _qn(W, "h"): "15840"})
    ET.SubElement(sect, _qn(W, "pgMar"), {
        _qn(W, "top"): "1080", _qn(W, "right"): "1080", _qn(W, "bottom"): "1080",
        _qn(W, "left"): "1080", _qn(W, "header"): "720", _qn(W, "footer"): "720",
    })

    rel_root = ET.Element("Relationships", xmlns=PKG_REL)
    for rel_id, target in relationships:
        ET.SubElement(rel_root, "Relationship", {
            "Id": rel_id,
            "Type": "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image",
            "Target": target,
        })
    package_rels = b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
</Relationships>'''
    content_types = b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/><Default Extension="png" ContentType="image/png"/>
<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
</Types>'''
    core = b'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title>ChampHLA journal-neutral submission</dc:title><dc:creator>ChampHLA authors</dc:creator><dcterms:created xsi:type="dcterms:W3CDTF">2000-01-01T00:00:00Z</dcterms:created><dcterms:modified xsi:type="dcterms:W3CDTF">2000-01-01T00:00:00Z</dcterms:modified></cp:coreProperties>'''
    parts = {
        "[Content_Types].xml": content_types,
        "_rels/.rels": package_rels,
        "docProps/core.xml": core,
        "word/document.xml": _xml_bytes(document),
        "word/styles.xml": _styles(),
        "word/_rels/document.xml.rels": _xml_bytes(rel_root),
    }
    parts.update({f"word/media/{name}": data for name, data in media})
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name in sorted(parts):
            info = zipfile.ZipInfo(name, (2000, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, parts[name])
    return {"blocks": blocks, "figure_order": figure_order, "sha256": _sha(output)}


def _block_signature(block: Block) -> tuple:
    if block.kind == "table":
        return ("table", block.rows)
    if block.kind == "image":
        return ("image", block.text)
    return (block.kind, block.text)


def _docx_signatures(path: Path) -> tuple[list[tuple], list[str], list[str]]:
    failures: list[str] = []
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        prohibited = [name for name in names if name.endswith("vbaProject.bin")
                      or "comments" in name or "people.xml" in name]
        if prohibited:
            failures.append(f"prohibited DOCX parts: {prohibited}")
        root = ET.fromstring(archive.read("word/document.xml"))
        signatures: list[tuple] = []
        figure_order: list[str] = []
        body = root.find("w:body", NS)
        for child in list(body):
            if child.tag == _qn(W, "sectPr"):
                continue
            if child.tag == _qn(W, "tbl"):
                rows = []
                for tr in child.findall("w:tr", NS):
                    rows.append(tuple(_plain("".join(tc.itertext()))
                                      for tc in tr.findall("w:tc", NS)))
                signatures.append(("table", tuple(rows)))
                continue
            style_node = child.find("w:pPr/w:pStyle", NS)
            style = style_node.get(_qn(W, "val")) if style_node is not None else ""
            text = _plain("".join(node.text or "" for node in child.findall(".//w:t", NS)))
            doc_pr = child.find(".//wp:docPr", NS)
            if doc_pr is not None:
                alt = doc_pr.get("descr", "")
                signatures.append(("image", alt))
                figure_order.append(alt)
            elif style == "Caption":
                continue
            elif style == "Title":
                signatures.append(("heading", text))
            elif style.startswith("Heading"):
                signatures.append(("heading", text))
            elif style == "ListBullet":
                signatures.append(("bullet", text[2:] if text.startswith("• ") else text))
            else:
                signatures.append(("paragraph", text))
        return signatures, figure_order, failures


def validate_docx_parity(markdown: Path, docx: Path) -> dict:
    expected_blocks = parse_markdown(markdown)
    expected = [_block_signature(block) for block in expected_blocks]
    actual, actual_figures, failures = _docx_signatures(docx)
    expected_figures = [block.text for block in expected_blocks if block.kind == "image"]
    if expected != actual:
        failures.append("Markdown/DOCX semantic block sequence differs")
    if expected_figures != actual_figures:
        failures.append("Markdown/DOCX figure order differs")
    return {
        "text_parity": expected == actual,
        "structure_valid": expected == actual and not any("prohibited" in item for item in failures),
        "figure_order_valid": expected_figures == actual_figures,
        "failures": failures,
    }


def audit_bibliography(markdown_paths: list[Path], bibliography: Path) -> dict:
    text = bibliography.read_text(encoding="utf-8")
    keys = re.findall(r"@\w+\s*\{\s*([^,\s]+)\s*,", text)
    citations: set[str] = set()
    for path in markdown_paths:
        citations.update(re.findall(r"(?<!\w)@([A-Za-z0-9_:.-]+)", path.read_text(encoding="utf-8")))
    duplicates = sorted({key for key in keys if keys.count(key) > 1})
    missing = sorted(citations - set(keys))
    uncited = sorted(set(keys) - citations)
    entries_without_doi = []
    for match in re.finditer(r"@\w+\s*\{\s*([^,\s]+)\s*,(.*?)(?=\n\})", text, re.S):
        if not re.search(r"\bdoi\s*=", match.group(2), re.I):
            entries_without_doi.append(match.group(1))
    passed = bool(keys) and not duplicates and not missing and not entries_without_doi
    return {
        "schema_version": "champhla-bibliography-audit-1", "passed": passed,
        "entry_count": len(keys), "citation_count": len(citations), "missing_keys": missing,
        "duplicate_keys": duplicates, "uncited_keys": uncited,
        "entries_without_doi": entries_without_doi,
    }


def audit_declarations(path: Path) -> dict:
    text = path.read_text(encoding="utf-8") if path.is_file() else ""
    required = ["Author list and affiliations", "Author contributions", "Competing interests",
                "Funding", "Ethics and consent", "Acknowledgements"]
    missing = [section for section in required if f"## {section}" not in text]
    markers = text.count("[AUTHOR INPUT REQUIRED]")
    return {
        "schema_version": "champhla-declarations-audit-1",
        "passed": not missing and markers == 0,
        "missing_sections": missing, "unresolved_author_markers": markers,
        "policy": "automation must not infer author identity, contributions, conflicts, funding, ethics, or acknowledgements",
    }


def build_submission_package(project_root: str | Path, main_markdown: str | Path,
                             supplement_markdown: str | Path, bibliography: str | Path,
                             declarations: str | Path, output_dir: str | Path,
                             manifest_output: str | Path, parity_output: str | Path,
                             bibliography_audit_output: str | Path,
                             declarations_audit_output: str | Path) -> dict:
    root = Path(project_root).resolve()
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    inputs = {
        "main": (root / main_markdown).resolve(),
        "supplement": (root / supplement_markdown).resolve(),
        "bibliography": (root / bibliography).resolve(),
        "declarations": (root / declarations).resolve(),
    }
    generated = {
        "main_markdown": output / "ChampHLA_manuscript.md",
        "main_docx": output / "ChampHLA_manuscript.docx",
        "supplement_markdown": output / "ChampHLA_supplement.md",
        "supplement_docx": output / "ChampHLA_supplement.docx",
        "bibliography": output / "references.bib",
        "declarations": output / "declarations_template.md",
    }
    for source_key, target_key in (("main", "main_markdown"), ("supplement", "supplement_markdown"),
                                   ("bibliography", "bibliography"), ("declarations", "declarations")):
        if inputs[source_key] != generated[target_key]:
            shutil.copyfile(inputs[source_key], generated[target_key])
    # Resolve figure paths against canonical source documents. The copied
    # Markdown is byte-identical, but its package directory intentionally has
    # a different relative path.
    render_docx(inputs["main"], generated["main_docx"])
    render_docx(inputs["supplement"], generated["supplement_docx"])
    main_parity = validate_docx_parity(generated["main_markdown"], generated["main_docx"])
    supplement_parity = validate_docx_parity(generated["supplement_markdown"], generated["supplement_docx"])
    failures = [f"main: {item}" for item in main_parity["failures"]]
    failures += [f"supplement: {item}" for item in supplement_parity["failures"]]
    parity = {
        "schema_version": "champhla-submission-parity-1",
        "passed": not failures,
        "text_parity": main_parity["text_parity"] and supplement_parity["text_parity"],
        "structure_valid": main_parity["structure_valid"] and supplement_parity["structure_valid"],
        "figure_order_valid": main_parity["figure_order_valid"] and supplement_parity["figure_order_valid"],
        "main": main_parity, "supplement": supplement_parity, "failures": failures,
    }
    bibliography_result = audit_bibliography(
        [generated["main_markdown"], generated["supplement_markdown"]], generated["bibliography"],
    )
    declarations_result = audit_declarations(generated["declarations"])
    write_json(parity_output, parity)
    write_json(bibliography_audit_output, bibliography_result)
    write_json(declarations_audit_output, declarations_result)
    fields = ["role", "path", "sha256", "bytes"]
    manifest_path = Path(manifest_output)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for role, path in generated.items():
            try:
                recorded_path = path.relative_to(root).as_posix()
            except ValueError:
                recorded_path = path.relative_to(output).as_posix()
            writer.writerow({"role": role, "path": recorded_path,
                             "sha256": _sha(path), "bytes": path.stat().st_size})
    return {"parity": parity, "bibliography": bibliography_result,
            "declarations": declarations_result, "files": generated}
