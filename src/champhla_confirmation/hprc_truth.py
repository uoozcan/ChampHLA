from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re

from .io import canonical_allele, read_tsv, sha256, write_json, write_tsv
from .panels import GENES


TRUTH_METHODS = {"HLA-ASM", "Immuannot"}
HEX64 = re.compile(r"^[0-9a-fA-F]{64}$")


def build_hprc_assembly_truth(calls_path: str | Path, truth_output: str | Path,
                              audit_output: str | Path) -> dict:
    rows = read_tsv(calls_path)
    required = {
        "sample_id", "haplotype", "gene", "method", "allele", "exon2_complete",
        "exon3_complete", "equally_supported_conflict", "source_sha256",
    }
    if not rows or required - set(rows[0]):
        raise ValueError(f"assembly truth calls missing columns: {sorted(required - set(rows[0]) if rows else required)}")
    indexed = defaultdict(list)
    for row in rows:
        gene = row["gene"].replace("HLA-", "").upper()
        if gene not in GENES:
            continue
        method = row["method"]
        if method not in TRUTH_METHODS:
            raise ValueError(f"unsupported assembly truth method: {method}")
        haplotype = row["haplotype"]
        if haplotype not in {"1", "2"}:
            raise ValueError(f"unsupported assembly haplotype: {haplotype}")
        if not HEX64.fullmatch(row.get("source_sha256", "")):
            raise ValueError(f"invalid assembly source SHA-256 for {row['sample_id']} {method}")
        indexed[(row["sample_id"], gene, method, haplotype)].append(row)
    samples = sorted({row["sample_id"] for row in rows})
    truth_rows, audit_rows = [], []
    for sample in samples:
        for gene in GENES:
            method_pairs = {}
            reasons = []
            source_hashes = set()
            for method in sorted(TRUTH_METHODS):
                alleles = []
                for haplotype in ("1", "2"):
                    matches = indexed.get((sample, gene, method, haplotype), [])
                    if len(matches) != 1:
                        reasons.append(f"{method}:haplotype_{haplotype}_record_count={len(matches)}")
                        continue
                    row = matches[0]
                    source_hashes.add(row["source_sha256"])
                    if row["exon2_complete"].lower() not in {"1", "true", "yes"}:
                        reasons.append(f"{method}:haplotype_{haplotype}:exon2_incomplete")
                    if row["exon3_complete"].lower() not in {"1", "true", "yes"}:
                        reasons.append(f"{method}:haplotype_{haplotype}:exon3_incomplete")
                    if row["equally_supported_conflict"].lower() in {"1", "true", "yes"}:
                        reasons.append(f"{method}:haplotype_{haplotype}:conflicting_allele")
                    try:
                        alleles.append(canonical_allele(row["allele"], gene))
                    except ValueError:
                        reasons.append(f"{method}:haplotype_{haplotype}:unresolved_allele")
                if len(alleles) == 2:
                    method_pairs[method] = tuple(sorted(alleles))
            resolved = not reasons and len(method_pairs) == 2 and len(set(method_pairs.values())) == 1
            if len(method_pairs) == 2 and len(set(method_pairs.values())) != 1:
                reasons.append("assembly_methods_disagree_at_two_field")
            pair = next(iter(method_pairs.values())) if resolved else ("", "")
            status = "resolved" if resolved else "unresolved"
            truth_rows.append({
                "cohort": "HPRC_RELEASE2_WGS", "subject": sample, "donor": sample,
                "independence_stratum": "donor_independent", "modality": "wgs",
                "gene": gene, "truth_allele1": pair[0], "truth_allele2": pair[1],
                "truth_status": status, "truth_source": "phased_assembly_dual_method",
                "source_sha256": ",".join(sorted(source_hashes)),
            })
            audit_rows.append({
                "sample_id": sample, "gene": gene, "truth_status": status,
                "HLA_ASM_pair": ",".join(method_pairs.get("HLA-ASM", ())),
                "Immuannot_pair": ",".join(method_pairs.get("Immuannot", ())),
                "reason": ";".join(reasons),
            })
    write_tsv(truth_output, truth_rows)
    write_tsv(Path(audit_output).with_suffix(".tsv"), audit_rows)
    summary = {
        "schema_version": "champhla-hprc-assembly-truth-audit-1",
        "truth_blind_to_short_read_predictions": True,
        "methods": sorted(TRUTH_METHODS),
        "subjects": len(samples), "loci": len(truth_rows),
        "resolved_loci": sum(row["truth_status"] == "resolved" for row in truth_rows),
        "unresolved_loci": sum(row["truth_status"] != "resolved" for row in truth_rows),
        "input_sha256": sha256(calls_path), "truth_sha256": sha256(truth_output),
        "audit_tsv_sha256": sha256(Path(audit_output).with_suffix(".tsv")),
        "passed_120_subject_roster": len(samples) == 120,
    }
    write_json(audit_output, summary)
    return summary
