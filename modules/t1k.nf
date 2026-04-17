/*
 * T1K Module
 * HLA typing from WGS/RNA-seq FASTQ data using the T1K assembler-based approach.
 * Covers Class I + II: A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1
 * Reference: https://github.com/mourisl/T1K
 */

process T1K_LONGREADS {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/t1k", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val platform

    output:
    tuple val(sample_id), path("${sample_id}_t1k.txt"), emit: results
    tuple val(sample_id), path("t1k_out/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    def preset = platform == 'hifi' ? 'HiFi' : 'ONT'
    """
    mkdir -p t1k_out

    HLA_FA=""
    if [ -n "${params.t1k_hlaidx}" ] && [ -f "${params.t1k_hlaidx}/hla_dna.fa" ]; then
        HLA_FA="${params.t1k_hlaidx}/hla_dna.fa"
    elif [ -n "${params.t1k_hlaidx}" ] && [ -f "${params.t1k_hlaidx}/_dna_seq.fa" ]; then
        HLA_FA="${params.t1k_hlaidx}/_dna_seq.fa"
    else
        HLA_FA=\$(find /usr/local/share -maxdepth 2 -type f -name 'hla_dna.fa' 2>/dev/null | head -1 || true)
    fi
    echo "[T1K] HLA reference: \${HLA_FA:-missing}"

    if [ -f "\${HLA_FA}" ]; then
        run-t1k \
            -u ${reads} \
            -f "\${HLA_FA}" \
            --preset ${preset} \
            -t ${task.cpus} \
            -o t1k_out/${sample_id} \
            2>t1k_out/t1k_stderr.log
    else
        echo "Missing T1K reference FASTA" > t1k_out/t1k_stderr.log
        exit 1
    fi

    GENOTYPE_FILE=\$(ls t1k_out/${sample_id}_genotype.tsv t1k_out/${sample_id}.genotype.tsv 2>/dev/null | head -1 || true)
    if [ -n "\${GENOTYPE_FILE}" ] && [ -f "\${GENOTYPE_FILE}" ]; then
        python3 ${projectDir}/bin/parse_t1k_results.py \
            "\${GENOTYPE_FILE}" \
            "${sample_id}_t1k.txt"
    else
        echo "# T1K long-read results for ${sample_id} (platform: ${preset})" > ${sample_id}_t1k.txt
        echo "# No results generated" >> ${sample_id}_t1k.txt
        cat t1k_out/t1k_stderr.log >> ${sample_id}_t1k.txt || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: \$(run-t1k --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo "1.0.9")
        platform: "${preset}"
    END_VERSIONS
    """
}

process T1K_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/t1k", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_t1k.txt"), emit: results
    tuple val(sample_id), path("t1k_out/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p t1k_out

    HLA_FA=""
    if [ -n "${params.t1k_hlaidx}" ] && [ -f "${params.t1k_hlaidx}/hla_dna.fa" ]; then
        HLA_FA="${params.t1k_hlaidx}/hla_dna.fa"
    elif [ -n "${params.t1k_hlaidx}" ] && [ -f "${params.t1k_hlaidx}/_dna_seq.fa" ]; then
        HLA_FA="${params.t1k_hlaidx}/_dna_seq.fa"
    else
        HLA_FA=\$(find /usr/local/share -maxdepth 2 -type f -name 'hla_dna.fa' 2>/dev/null | head -1 || true)
    fi
    echo "[T1K] HLA reference: \${HLA_FA:-missing}"

    F1="${fastq1}"
    F2="${fastq2}"
    T1K_PRESET=\$([ "${params.seq_type}" = "rna" ] && echo "hla" || echo "hla-wgs")

    if [ -f "\${HLA_FA}" ]; then
        run-t1k \
            -1 "\${F1}" \
            -2 "\${F2}" \
            -f "\${HLA_FA}" \
            --preset "\${T1K_PRESET}" \
            -t ${task.cpus} \
            -o t1k_out/${sample_id} \
            2>t1k_out/t1k_stderr.log
    else
        echo "Missing T1K reference FASTA" > t1k_out/t1k_stderr.log
        exit 1
    fi

    GENOTYPE_FILE=\$(ls t1k_out/${sample_id}_genotype.tsv t1k_out/${sample_id}.genotype.tsv 2>/dev/null | head -1 || true)
    if [ -n "\${GENOTYPE_FILE}" ] && [ -f "\${GENOTYPE_FILE}" ]; then
        python3 ${projectDir}/bin/parse_t1k_results.py \
            "\${GENOTYPE_FILE}" \
            "${sample_id}_t1k.txt"
    else
        echo "# T1K results for ${sample_id}" > ${sample_id}_t1k.txt
        echo "# No results generated" >> ${sample_id}_t1k.txt
        cat t1k_out/t1k_stderr.log >> ${sample_id}_t1k.txt || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: \$(run-t1k --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo "1.0.9")
    END_VERSIONS
    """
}
