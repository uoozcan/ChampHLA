process ARCASHLA_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/arcashla", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_arcashla.txt"), emit: results
    path("${sample_id}_arcashla.json"), emit: json_results
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    mkdir -p ${sample_id} arcas_patch
    if [[ "${fastq1}" == *.gz ]]; then
        ln -s \$(realpath ${fastq1}) ${sample_id}/${sample_id}.1.fq.gz
        ln -s \$(realpath ${fastq2}) ${sample_id}/${sample_id}.2.fq.gz
    else
        gzip -c ${fastq1} > ${sample_id}/${sample_id}.1.fq.gz
        gzip -c ${fastq2} > ${sample_id}/${sample_id}.2.fq.gz
    fi
    test -s ${sample_id}/${sample_id}.1.fq.gz
    test -s ${sample_id}/${sample_id}.2.fq.gz
    cp /home/arcasHLA-master/scripts/*.py arcas_patch/
    sed -i 's/count = counts\\[eq\\]/count = counts.get(eq, 0)/' arcas_patch/align.py
    export PYTHONPATH="\$(pwd)/arcas_patch:\${PYTHONPATH:-}"
    arcasHLA genotype ${sample_id}/${sample_id}.1.fq.gz ${sample_id}/${sample_id}.2.fq.gz \
        -o ${sample_id} -t ${task.cpus} --min_count 75 -v
    genotype_json=\$(find ${sample_id} -type f -name '*.genotype.json' -size +0c -print -quit)
    test -n "\${genotype_json}"
    cp "\${genotype_json}" ${sample_id}_arcashla.json
    python3 - "\${genotype_json}" "${sample_id}_arcashla.txt" <<'PY'
import json, sys
with open(sys.argv[1], encoding="utf-8") as handle:
    data = json.load(handle)
with open(sys.argv[2], "w", encoding="utf-8") as out:
    out.write("Gene\\tAllele1\\tAllele2\\n")
    for gene in ("A", "B", "C"):
        alleles = data.get(gene) or []
        a1 = alleles[0] if len(alleles) > 0 else "-"
        a2 = alleles[1] if len(alleles) > 1 else "-"
        out.write("HLA-{0}\\t{1}\\t{2}\\n".format(gene, a1, a2))
PY
    test -s ${sample_id}_arcashla.txt
    printf '"%s":\\n    arcashla: "container-pinned"\\n' "${task.process}" > versions.yml
    """
}
