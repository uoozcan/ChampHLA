#!/usr/bin/env python3
"""Build compact, truth-separated public metadata snapshots for dataset discovery.

The script consumes caller-supplied downloads from official APIs/pages. It does
not download sequence reads and never writes HLA allele values into the sample
crosswalk or pilot manifest.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import zipfile
from collections import defaultdict
from html.parser import HTMLParser
from pathlib import Path

from champhla_confirmation.io import read_tsv, sha256, write_json, write_tsv
from champhla_confirmation.dataset_sources import (
    NCI_RNA_ALIAS,
    RELATED_GROUPS,
    locus_has_exact_two_field_truth as _locus_exact,
    normalize_nci_name as _nci_name,
    normalize_sample_name as _norm_name,
)


class PmcTableParser(HTMLParser):
    def __init__(self, section_id: str):
        super().__init__()
        self.section_id = section_id
        self.in_section = False
        self.section_depth = 0
        self.in_row = False
        self.in_cell = False
        self.cell_parts: list[str] = []
        self.row: list[str] = []
        self.rows: list[list[str]] = []

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag == "section" and attributes.get("id") == self.section_id:
            self.in_section = True
            self.section_depth = 1
            return
        if self.in_section and tag == "section":
            self.section_depth += 1
        if not self.in_section:
            return
        if tag == "tr":
            self.in_row = True
            self.row = []
        elif tag in {"th", "td"} and self.in_row:
            self.in_cell = True
            self.cell_parts = []

    def handle_endtag(self, tag):
        if not self.in_section:
            return
        if tag in {"th", "td"} and self.in_cell:
            value = " ".join("".join(self.cell_parts).split())
            self.row.append(value)
            self.in_cell = False
        elif tag == "tr" and self.in_row:
            if self.row:
                self.rows.append(self.row)
            self.in_row = False
        elif tag == "section":
            self.section_depth -= 1
            if self.section_depth == 0:
                self.in_section = False

    def handle_data(self, data):
        if self.in_cell:
            self.cell_parts.append(data)


def _extract_nci_truth_metadata(html_path: str) -> list[dict[str, str]]:
    parser = PmcTableParser("T2")
    parser.feed(Path(html_path).read_text(encoding="utf-8", errors="replace"))
    header = parser.rows[0]
    rows = []
    for values in parser.rows[1:]:
        if len(values) != len(header):
            continue
        row = dict(zip(header, values))
        subject = _nci_name(row["Cell Line"])
        exact = all(_locus_exact(row[field]) for field in ("A locus", "B Locus", "Cw Locus"))
        related_group, relationship = RELATED_GROUPS.get(subject, ("", "unrelated_or_unresolved"))
        rows.append({
            "subject": subject,
            "truth_identifier": row["ID"],
            "tissue": row["Tissue"],
            "complete_abc_two_field": "1" if exact else "0",
            "ambiguity_status": "complete_two_field" if exact else "contains_one_field_or_unresolved",
            "related_group": related_group,
            "relationship": relationship,
            "truth_source": "PMC555742_Table_2",
        })
    return rows


def _ena_rows(path: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _nci_rna_snapshot(rows: list[dict[str, str]], truth_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    truth_by_norm = {_norm_name(row["subject"]): row for row in truth_rows}
    output = []
    for row in rows:
        sample_title = row["sample_title"]
        raw_name = sample_title.split(":", 1)[-1]
        truth = truth_by_norm.get(_norm_name(NCI_RNA_ALIAS.get(raw_name, raw_name)))
        if not truth:
            continue
        urls = [f"https://{value}" for value in row["fastq_ftp"].split(";") if value]
        output.append({
            "subject": truth["subject"], "sample_title": sample_title,
            "run_accession": row["run_accession"], "study_accession": row["study_accession"],
            "sample_accession": row["sample_accession"],
            "experiment_accession": row["experiment_accession"],
            "library_strategy": row["library_strategy"], "library_source": row["library_source"],
            "library_layout": row["library_layout"], "instrument_platform": row["instrument_platform"],
            "base_count": row["base_count"], "read_count": row["read_count"],
            "read1_url": urls[0] if urls else "", "read2_url": urls[1] if len(urls) > 1 else "",
        })
    return sorted(output, key=lambda row: row["subject"])


def _ccle_overlap_snapshot(rows: list[dict[str, str]], truth_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    truth_by_norm = {_norm_name(row["subject"]): row for row in truth_rows}
    output = []
    for row in rows:
        # CCLE titles append a tissue label that may itself contain underscores.
        sample_name = row["sample_title"].split("_", 1)[0]
        truth = truth_by_norm.get(_norm_name(sample_name))
        if not truth:
            continue
        urls = [f"https://{value}" for value in row["fastq_ftp"].split(";") if value]
        output.append({
            "subject": truth["subject"], "sample_title": row["sample_title"],
            "run_accession": row["run_accession"], "study_accession": row["study_accession"],
            "sample_accession": row["sample_accession"],
            "experiment_accession": row["experiment_accession"],
            "library_strategy": row["library_strategy"], "library_source": row["library_source"],
            "library_layout": row["library_layout"], "instrument_platform": row["instrument_platform"],
            "base_count": row["base_count"], "read_count": row["read_count"],
            "read1_url": urls[0] if urls else "", "read2_url": urls[1] if len(urls) > 1 else "",
        })
    return sorted(output, key=lambda row: row["subject"])


def _afgr_samples(zip_path: str) -> list[dict[str, str]]:
    member = "biorxiv-AFGR-SupplementaryTables/biorxiv-AFGR-SupplementaryTableS1.tsv"
    with zipfile.ZipFile(zip_path) as archive:
        text = archive.read(member).decode("utf-8-sig")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    return [row for row in rows if row["Population"] == "MKK" and row["Transcriptome"] == "yes"]


def _afgr_snapshot(zip_path: str, encode_json_path: str) -> list[dict[str, str]]:
    experiments = json.loads(Path(encode_json_path).read_text(encoding="utf-8"))["@graph"]
    rna_by_coriell = {}
    for experiment in experiments:
        if experiment.get("assay_term_name") != "RNA-seq":
            continue
        name = experiment.get("biosample_ontology", {}).get("term_name", "")
        rna_by_coriell[name] = experiment
    output = []
    for row in _afgr_samples(zip_path):
        coriell = row["LCL Coriell ID"]
        experiment = rna_by_coriell.get(coriell, {})
        file_accessions = [item.get("@id", "").strip("/").split("/")[-1]
                           for item in experiment.get("files", [])]
        output.append({
            "donor": row["1000G ID"], "coriell_id": coriell, "population": row["Population"],
            "sex": row["Gender"], "genome_reported": row["Genome"],
            "transcriptome_reported": row["Transcriptome"], "atac_seq_reported": row["ATAC-Seq"],
            "encode_experiment": experiment.get("accession", ""),
            "encode_biosample": (experiment.get("replicates", [{}])[0].get("library", {})
                                 .get("biosample", {}).get("accession", "")
                                 if experiment.get("replicates") else ""),
            "encode_file_accessions": ",".join(file_accessions),
            "encode_status": experiment.get("status", "missing"),
        })
    return sorted(output, key=lambda row: row["donor"])


def _dev_tokens(path: str) -> set[str]:
    result = set()
    for row in read_tsv(path):
        for field in ("subject", "sample", "sample_id", "aliases", "coriell_id"):
            result.update(_norm_name(value) for value in row.get(field, "").replace(";", ",").split(",")
                          if value.strip())
    return result


def _crosswalk(nci_truth, nci_rna, ccle_wes, afgr, development_path):
    truth_by_subject = {row["subject"]: row for row in nci_truth}
    dev = _dev_tokens(development_path)
    rows = []
    for modality, records, cohort in (("rnaseq", nci_rna, "NCI60"),
                                      ("wes", ccle_wes, "NCI60_CCLE")):
        for record in records:
            truth = truth_by_subject[record["subject"]]
            related = truth["relationship"]
            exclusion = ""
            gate = "ready" if truth["complete_abc_two_field"] == "1" else "blocked_truth_ambiguity"
            if gate == "blocked_truth_ambiguity":
                exclusion = "strict A/B/C truth contains a one-field, ambiguous, novel, or unresolved call"
            if related == "derivative_or_duplicate":
                gate, exclusion = "excluded_related", f"related_group={truth['related_group']}"
            overlap = "exact" if _norm_name(record["subject"]) in dev else "none"
            rows.append({
                "cohort": cohort, "repository_sample": record["sample_accession"],
                "donor": record["subject"], "aliases": record["sample_title"], "modality": modality,
                "sequencing_accession": record["run_accession"],
                "biosample_accession": record["sample_accession"],
                "experiment_accession": record["experiment_accession"],
                "file_accessions": record["run_accession"],
                "source_urls": ";".join((record["read1_url"], record["read2_url"])),
                "truth_identifier": truth["truth_identifier"], "sample_link_status": "verified",
                "existing_champhla_overlap": overlap, "relatedness": related,
                "independence_stratum": "donor_independent" if overlap == "none" else "new_library_overlap",
                "access_class": "public", "gate_status": gate, "exclusion_reason": exclusion,
            })
    for record in afgr:
        overlap = "exact" if (_norm_name(record["donor"]) in dev
                               or _norm_name(record["coriell_id"]) in dev) else "none"
        rows.append({
            "cohort": "AFGR_MKK", "repository_sample": record["encode_biosample"],
            "donor": record["donor"], "aliases": record["coriell_id"], "modality": "rnaseq",
            "sequencing_accession": record["encode_experiment"],
            "biosample_accession": record["encode_biosample"],
            "experiment_accession": record["encode_experiment"],
            "file_accessions": record["encode_file_accessions"], "source_urls": "",
            "truth_identifier": record["donor"], "sample_link_status": "unresolved_truth_file_mapping",
            "existing_champhla_overlap": overlap, "relatedness": "requires_coriell_family_audit",
            "independence_stratum": "donor_independent" if overlap == "none" else "new_library_overlap",
            "access_class": "mixed", "gate_status": "blocked_truth_crosswalk",
            "exclusion_reason": "EGA truth file not opened; donor-level mapping unverified",
        })
        if record["genome_reported"] == "yes":
            rows.append({
                "cohort": "AFGR_MKK", "repository_sample": record["donor"],
                "donor": record["donor"], "aliases": record["coriell_id"], "modality": "wgs",
                "sequencing_accession": "", "biosample_accession": "", "experiment_accession": "",
                "file_accessions": "", "source_urls": "", "truth_identifier": record["donor"],
                "sample_link_status": "literature_only", "existing_champhla_overlap": overlap,
                "relatedness": "requires_coriell_family_audit",
                "independence_stratum": "donor_independent" if overlap == "none" else "new_library_overlap",
                "access_class": "unresolved", "gate_status": "blocked_raw_accession",
                "exclusion_reason": "AFGR reports WGS but ENCODE AFGR collection has no WGS experiment",
            })
    return sorted(rows, key=lambda row: (row["cohort"], row["donor"], row["modality"]))


def _pilot(nci_truth, crosswalk):
    def tissue_group(value: str) -> str:
        return value.split(",", 1)[0].replace(" CA", "").strip()

    truth = {row["subject"]: row for row in nci_truth}
    by_key = {(row["cohort"], row["donor"], row["modality"]): row for row in crosswalk}
    complete = [row for row in nci_truth if row["complete_abc_two_field"] == "1"
                and row["relationship"] != "derivative_or_duplicate"]
    nci_subjects = sorted(row["subject"] for row in complete)
    represented_tissues = {tissue_group(row["tissue"]) for row in complete}
    dual_candidates = []
    for subject in sorted(truth):
        rna = by_key.get(("NCI60", subject, "rnaseq"))
        wes = by_key.get(("NCI60_CCLE", subject, "wes"))
        if rna and wes and subject not in nci_subjects and truth[subject]["relationship"] != "derivative_or_duplicate":
            dual_candidates.append(
                (subject, tissue_group(truth[subject]["tissue"]) not in represented_tissues)
            )
    if len(nci_subjects) < 12 and dual_candidates:
        # Select without allele values: prefer a new tissue among public dual-modality lines.
        nci_subjects.append(sorted(dual_candidates, key=lambda item: (-item[1], item[0]))[0][0])

    rows = []
    for subject in nci_subjects:
        for modality, cohort in (("rnaseq", "NCI60"), ("wes", "NCI60_CCLE")):
            source = by_key.get((cohort, subject, modality))
            if not source:
                continue
            urls = source["source_urls"].split(";")
            for gene in ("A", "B", "C"):
                rows.append({
                    "cohort": "NCI60_PUBLIC_PILOT", "subject": subject, "modality": modality,
                    "gene": gene, "source": "ENA/SRA",
                    "source_accession": source["sequencing_accession"], "input_type": "fastq",
                    "read1_url": urls[0], "read2_url": urls[1] if len(urls) > 1 else "",
                    "access_class": "public", "independence_stratum": source["independence_stratum"],
                    "gate_status": source["gate_status"],
                    "exclusion_reason": source["exclusion_reason"] or (
                        "strict A/B/C truth contains unresolved one-field alleles"
                        if source["gate_status"] == "blocked_truth_ambiguity" else ""
                    ),
                })

    mkk_rna = [row for row in crosswalk if row["cohort"] == "AFGR_MKK"
               and row["modality"] == "rnaseq" and row["existing_champhla_overlap"] == "none"]
    for source in sorted(mkk_rna, key=lambda row: row["donor"])[:12]:
        wgs = by_key.get(("AFGR_MKK", source["donor"], "wgs"))
        for modality, record in (("rnaseq", source), ("wgs", wgs)):
            if not record:
                continue
            for gene in ("A", "B", "C"):
                rows.append({
                    "cohort": "AFGR_MKK_GATED_PILOT", "subject": source["donor"],
                    "modality": modality, "gene": gene,
                    "source": "ENCODE" if modality == "rnaseq" else "AFGR_publication",
                    "source_accession": record["sequencing_accession"],
                    "input_type": "fastq" if modality == "rnaseq" else "unresolved",
                    "read1_url": "", "read2_url": "", "access_class": record["access_class"],
                    "independence_stratum": record["independence_stratum"],
                    "gate_status": record["gate_status"], "exclusion_reason": record["exclusion_reason"],
                })
    return sorted(rows, key=lambda row: (row["cohort"], row["subject"], row["modality"], row["gene"]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nci-hla-html", required=True)
    parser.add_argument("--nci-rna-ena", required=True)
    parser.add_argument("--ccle-wes-ena", required=True)
    parser.add_argument("--afgr-supplement-zip", required=True)
    parser.add_argument("--afgr-encode-json", required=True)
    parser.add_argument("--development", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    out = Path(args.output_dir)
    sources = out / "sources"

    nci_truth = _extract_nci_truth_metadata(args.nci_hla_html)
    nci_rna = _nci_rna_snapshot(_ena_rows(args.nci_rna_ena), nci_truth)
    ccle_wes = _ccle_overlap_snapshot(_ena_rows(args.ccle_wes_ena), nci_truth)
    afgr = _afgr_snapshot(args.afgr_supplement_zip, args.afgr_encode_json)
    crosswalk = _crosswalk(nci_truth, nci_rna, ccle_wes, afgr, args.development)
    pilot = _pilot(nci_truth, crosswalk)

    write_tsv(sources / "nci60_truth_eligibility.tsv", nci_truth)
    write_tsv(sources / "nci60_rna_ena.tsv", nci_rna)
    write_tsv(sources / "ccle_nci60_wes_ena.tsv", ccle_wes)
    write_tsv(sources / "afgr_mkk_encode.tsv", afgr)
    write_tsv(out / "sample_crosswalk.tsv", crosswalk)
    write_tsv(out / "pilot_manifest.tsv", pilot)
    write_json(sources / "source_manifest.json", {
        "schema_version": "dataset-discovery-source-snapshot-1",
        "generated_from_official_metadata": True,
        "raw_reads_downloaded": False,
        "input_sha256": {
            "nci_hla_html": sha256(args.nci_hla_html),
            "nci_rna_ena": sha256(args.nci_rna_ena),
            "ccle_wes_ena": sha256(args.ccle_wes_ena),
            "afgr_supplement_zip": sha256(args.afgr_supplement_zip),
            "afgr_encode_json": sha256(args.afgr_encode_json),
            "development": sha256(args.development),
        },
        "source_urls": {
            "nci_hla": "https://pmc.ncbi.nlm.nih.gov/articles/PMC555742/",
            "nci_rna": "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=study_accession%3D%22PRJNA433861%22",
            "ccle_wes": "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=study_accession%3D%22PRJNA523380%22%20AND%20library_strategy%3D%22WXS%22",
            "afgr_supplement": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10659267/supplementaryFiles",
            "afgr_encode": "https://www.encodeproject.org/search/?type=Experiment&searchTerm=AFGR&format=json&limit=all",
        },
        "snapshot_counts": {
            "nci_truth_subjects": len(nci_truth), "nci_rna_subjects": len(nci_rna),
            "ccle_nci60_wes_subjects": len(ccle_wes), "afgr_mkk_rna_subjects": len(afgr),
            "crosswalk_rows": len(crosswalk), "pilot_rows": len(pilot),
        },
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
