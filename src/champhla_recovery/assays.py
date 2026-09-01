from __future__ import annotations

import csv
import re
from pathlib import Path

from .io import sha256, write_json, write_tsv


SUBJECT_FROM_WES = re.compile(r"/data/[^/]+/((?:HG|NA)\d{5})/", re.I)
POP_FROM_WES = re.compile(r"/data/([^/]+)/", re.I)


def _https(url: str) -> str:
    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://"):]
    if url.startswith("ftp:/"):
        return "https://" + url[len("ftp:/"):]
    return url


def _comment_table(path: str, header_prefix: str) -> list[dict[str, str]]:
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith(header_prefix))
    header = lines[index].lstrip("#").split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[index + 1:] if line.strip()]


def build_assay_manifest(wgs_index: str, wes_index: str, rna_sdrf: str,
                         output: str, manifest: str) -> dict:
    rows = []
    for row in _comment_table(wgs_index, "#ENA_FILE_PATH"):
        url = _https(row["ENA_FILE_PATH"])
        rows.append({
            "sample_id": row["SAMPLE_NAME"].upper(), "modality": "wgs",
            "population": row.get("POPULATION", ""), "format": "CRAM",
            "genome_build": "GRCh38DH", "input_url": url,
            "index_url": url + ".crai", "input_md5": row.get("MD5SUM", ""),
            "index_md5": "", "source_record": row.get("RUN_ID", ""),
        })
    for row in _comment_table(wes_index, "#CRAM"):
        url = _https(row["CRAM"])
        subject_match = SUBJECT_FROM_WES.search(url)
        population_match = POP_FROM_WES.search(url)
        if not subject_match:
            continue
        rows.append({
            "sample_id": subject_match.group(1).upper(), "modality": "wes",
            "population": population_match.group(1) if population_match else "",
            "format": "CRAM", "genome_build": "GRCh38DH", "input_url": url,
            "index_url": _https(row.get("CRAI", "")), "input_md5": row.get("CRAM_MD5", ""),
            "index_md5": row.get("CRAI_MD5", ""), "source_record": Path(url).name,
        })
    with Path(rna_sdrf).open(newline="", encoding="utf-8", errors="replace") as handle:
        sdrf = list(csv.DictReader(handle, delimiter="\t"))
    seen = set()
    for row in sdrf:
        subject = row.get("Characteristics[individual]", "").upper()
        assay = row.get("Assay Name", "")
        if not subject or subject in seen:
            continue
        seen.add(subject)
        base = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/processed/"
        rows.append({
            "sample_id": subject, "modality": "rnaseq",
            "population": row.get("Factor Value[ancestry category]", ""),
            "format": "BAM", "genome_build": "GRCh37",
            "input_url": base + assay + ".bam", "index_url": base + assay + ".bam.bai",
            "input_md5": "", "index_md5": "", "source_record": row.get("Comment[ENA_RUN]", ""),
        })
    rows.sort(key=lambda row: (row["modality"], row["sample_id"]))
    write_tsv(output, rows, [
        "sample_id", "modality", "population", "format", "genome_build", "input_url",
        "index_url", "input_md5", "index_md5", "source_record",
    ])
    result = {
        "schema_version": "official-1000g-assay-manifest-1",
        "truth_blind": True,
        "counts": {modality: sum(row["modality"] == modality for row in rows)
                   for modality in ("wgs", "wes", "rnaseq")},
        "source_sha256": {"wgs": sha256(wgs_index), "wes": sha256(wes_index), "rnaseq": sha256(rna_sdrf)},
        "output_sha256": sha256(output),
    }
    write_json(manifest, result)
    return result

