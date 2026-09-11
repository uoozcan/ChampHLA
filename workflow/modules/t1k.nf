process T1K_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/t1k", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_t1k.txt"), emit: results
    path "versions.yml", emit: versions

    script:
    def reference_name = params.seq_type == 'rna' ? '_rna_seq.fa' : '_dna_seq.fa'
    def preset = params.seq_type == 'rna' ? 'hla' : 'hla-wgs'
    """
    set -euo pipefail
    mkdir -p t1k_out
    reference=${params.t1k_hlaidx}/${reference_name}
    test -s "\${reference}"
    test -s ${fastq1}
    test -s ${fastq2}
    run-t1k -1 ${fastq1} -2 ${fastq2} -f "\${reference}" --preset ${preset} \
        -t ${task.cpus} -o t1k_out/${sample_id} 2>t1k_out/t1k_stderr.log
    if [ -s t1k_out/${sample_id}_genotype.tsv ]; then
        genotype=t1k_out/${sample_id}_genotype.tsv
    elif [ -s t1k_out/${sample_id}.genotype.tsv ]; then
        genotype=t1k_out/${sample_id}.genotype.tsv
    else
        echo "T1K produced no genotype result" >&2
        exit 4
    fi
    python3 ${projectDir}/bin/parse_t1k_results.py "\${genotype}" ${sample_id}_t1k.txt
    test -s ${sample_id}_t1k.txt
    printf '"%s":\n    t1k: "1.0.9-r251"\n' "${task.process}" > versions.yml
    """
}
