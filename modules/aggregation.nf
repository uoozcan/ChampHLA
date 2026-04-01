/*
 * Aggregate per-tool HLA result files into a harmonized table for consensus.
 */

process AGGREGATE_RESULTS {
    tag "aggregate"
    label 'process_low'
    publishDir "${params.outdir}/consensus", mode: 'copy'

    input:
    path result_files

    output:
    path "aggregated_calls.tsv", emit: calls

    script:
    def modality = params.run_modality ?: (params.input_type == 'fastq' ? (params.seq_type == 'rna' ? 'rnaseq' : 'wes') : 'wgs')
    """
    printf '%s\n' ${result_files} > result_files.list
    python3 - << 'PYEOF'
import csv
import re
from pathlib import Path

with open("result_files.list", "r", encoding="utf-8") as fh:
    files = [Path(line.strip()) for line in fh if line.strip()]
modality = "${modality}"

# Expected naming from modules: <sample>_<tool>.txt
name_patterns = {
    "optitype": re.compile(r"^(?P<sample>.+)_optitype\.txt$"),
    "arcashla": re.compile(r"^(?P<sample>.+)_arcashla\.txt$"),
    "spechla": re.compile(r"^(?P<sample>.+)_spechla\.txt$"),
    "hlahd": re.compile(r"^(?P<sample>.+)_hlahd\.txt$"),
    "polysolver": re.compile(r"^(?P<sample>.+)_polysolver\.txt$"),
    "kourami": re.compile(r"^(?P<sample>.+)_kourami\.txt$"),
    "t1k": re.compile(r"^(?P<sample>.+)_t1k\.txt$"),
    "seq2hla": re.compile(r"^(?P<sample>.+)_seq2hla\.txt$"),
}

rows = []
for f in sorted(files):
    tool = None
    sample = None
    for t, pat in name_patterns.items():
        m = pat.match(f.name)
        if m:
            tool = t
            sample = m.group("sample")
            break
    if not tool:
        continue
    if not f.exists() or f.stat().st_size == 0:
        continue

    with f.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            gene = parts[0].strip().replace("HLA-", "")
            if gene.lower() == "gene":
                continue
            a1 = parts[1].strip()
            a2 = parts[2].strip()
            missing = {"", "-", "NA", "None", "none", "."}
            is_callable = int((a1 not in missing) and (a2 not in missing))
            rows.append({
                "sample": sample,
                "tool": tool,
                "modality": modality,
                "gene": gene,
                "allele1": a1,
                "allele2": a2,
                "is_callable": str(is_callable),
            })

out = Path("aggregated_calls.tsv")
with out.open("w", encoding="utf-8", newline="") as handle:
    fields = ["sample", "tool", "modality", "gene", "allele1", "allele2", "is_callable"]
    w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for row in rows:
        w.writerow(row)
PYEOF
    """
}
