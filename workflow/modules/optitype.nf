process OPTITYPE_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/optitype", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)
    val seq_type

    output:
    tuple val(sample_id), path("${sample_id}_optitype.txt"), emit: results
    path "versions.yml", emit: versions

    script:
    def type_flag = seq_type == 'rna' ? '--rna' : '--dna'
    """
    set -euo pipefail
    mkdir -p ${sample_id}
    trap 'rm -f R1.fastq R2.fastq' EXIT
    if [[ "${fastq1}" == *.gz ]]; then
        zcat ${fastq1} > R1.fastq
        zcat ${fastq2} > R2.fastq
    else
        ln -s \$(realpath ${fastq1}) R1.fastq
        ln -s \$(realpath ${fastq2}) R2.fastq
    fi
    test -s R1.fastq
    test -s R2.fastq
    OptiTypePipeline.py -i R1.fastq R2.fastq ${type_flag} -v -o ${sample_id} -p ${sample_id}
    result_tsv=\$(find ${sample_id} -type f -name '*_result.tsv' -size +0c -print -quit)
    test -n "\${result_tsv}"
    python3 - "\${result_tsv}" "${sample_id}_optitype.txt" <<'PY'
import csv, sys
with open(sys.argv[1], encoding="utf-8") as handle:
    row = next(csv.DictReader(handle, delimiter="\\t"))
with open(sys.argv[2], "w", encoding="utf-8") as out:
    out.write("Gene\\tAllele1\\tAllele2\\n")
    for gene in ("A", "B", "C"):
        a1 = row.get(gene + "1") or "-"
        a2 = row.get(gene + "2") or "-"
        out.write("HLA-{0}\\t{1}\\t{2}\\n".format(gene, a1, a2))
PY
    test -s ${sample_id}_optitype.txt
    printf '"%s":\\n    optitype: "container-pinned"\\n' "${task.process}" > versions.yml
    """
}
