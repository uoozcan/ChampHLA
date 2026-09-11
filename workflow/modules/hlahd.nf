process HLAHD {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/hlahd", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    val reference
    val hla_genes

    output:
    tuple val(sample_id), path("${sample_id}_hlahd.txt"), emit: results
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    test -x /app/hlahd.1.4.0/bin/hlahd.sh
    test -s ${params.hlahd_db}/HLA_gene.split.txt
    test -d ${params.hlahd_db}/dictionary
    test -d ${params.hlahd_db}/freq_data
    mkdir -p ${sample_id}
    samtools quickcheck -v ${bam}
    samtools view -H ${bam} | grep -Eq '^@SQ.*SN:(chr6|6)[[:space:]]'
    chr=\$(samtools view -H ${bam} | awk '/^@SQ.*SN:chr6[[:space:]]/{print "chr6"; exit} /^@SQ.*SN:6[[:space:]]/{print "6"; exit}')
    samtools view -b -h ${bam} "\${chr}:${params.hla_region_start}-${params.hla_region_end}" > hla_region.bam
    samtools view -b -f 4 ${bam} > unmapped.bam
    samtools merge -f merged.bam hla_region.bam unmapped.bam
    samtools sort -n -@ ${task.cpus} merged.bam -o sorted.bam
    samtools fastq -@ ${task.cpus} -1 R1.fastq -2 R2.fastq -0 /dev/null -s /dev/null sorted.bam
    test -s R1.fastq
    test -s R2.fastq
    read_len=\$(awk 'NR==2{print length(\$0); exit}' R1.fastq)
    test "\${read_len}" -gt 0
    min_tag=\$(( read_len < 100 ? read_len / 2 : 50 ))
    /app/hlahd.1.4.0/bin/hlahd.sh -t ${task.cpus} -m "\${min_tag}" -c 0.95 \
        -f ${params.hlahd_db}/freq_data R1.fastq R2.fastq \
        ${params.hlahd_db}/HLA_gene.split.txt ${params.hlahd_db}/dictionary \
        ${sample_id} ${sample_id}
    native=${sample_id}/${sample_id}/result/${sample_id}_final.result.txt
    test -s "\${native}"
    cp "\${native}" ${sample_id}_hlahd.txt
    printf '"%s":\n    hlahd: "1.4.0"\n' "${task.process}" > versions.yml
    rm -f hla_region.bam unmapped.bam merged.bam sorted.bam R1.fastq R2.fastq
    """
}

process HLAHD_FASTQ {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/hlahd", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)
    val hla_genes

    output:
    tuple val(sample_id), path("${sample_id}_hlahd.txt"), emit: results
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    test -x /app/hlahd.1.4.0/bin/hlahd.sh
    test -s ${params.hlahd_db}/HLA_gene.split.txt
    test -d ${params.hlahd_db}/dictionary
    test -d ${params.hlahd_db}/freq_data
    mkdir -p ${sample_id}
    if [[ "${fastq1}" == *.gz ]]; then
        zcat ${fastq1} > R1.fastq
        zcat ${fastq2} > R2.fastq
    else
        cp ${fastq1} R1.fastq
        cp ${fastq2} R2.fastq
    fi
    test -s R1.fastq
    test -s R2.fastq
    read_len=\$(awk 'NR==2{print length(\$0); exit}' R1.fastq)
    test "\${read_len}" -gt 0
    min_tag=\$(( read_len < 100 ? read_len / 2 : 50 ))
    /app/hlahd.1.4.0/bin/hlahd.sh -t ${task.cpus} -m "\${min_tag}" -c 0.95 \
        -f ${params.hlahd_db}/freq_data R1.fastq R2.fastq \
        ${params.hlahd_db}/HLA_gene.split.txt ${params.hlahd_db}/dictionary \
        ${sample_id} ${sample_id}
    native=${sample_id}/${sample_id}/result/${sample_id}_final.result.txt
    test -s "\${native}"
    cp "\${native}" ${sample_id}_hlahd.txt
    printf '"%s":\n    hlahd: "1.4.0"\n' "${task.process}" > versions.yml
    rm -f R1.fastq R2.fastq
    """
}
