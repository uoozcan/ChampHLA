process CRAM_TO_BAM {
    tag "$sample_id"
    label 'process_medium'
    input:
    tuple val(sample_id), path(cram)
    path reference_fasta
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    script:
    """
    set -euo pipefail
    samtools quickcheck -v ${cram}
    samtools view -@ ${task.cpus} -b -T ${reference_fasta} -o ${sample_id}.bam ${cram}
    samtools index -@ ${task.cpus} ${sample_id}.bam
    samtools quickcheck -v ${sample_id}.bam
    """
}

process BAM_TO_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    input:
    tuple val(sample_id), path(bam)
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: reads
    script:
    """
    set -euo pipefail
    samtools quickcheck -v ${bam}
    samtools sort -n -@ ${task.cpus} ${bam} -o ${sample_id}.namesort.bam
    samtools fastq -@ ${task.cpus} -1 ${sample_id}_R1.fastq.gz -2 ${sample_id}_R2.fastq.gz \
        -0 /dev/null -s /dev/null ${sample_id}.namesort.bam
    test -s ${sample_id}_R1.fastq.gz
    test -s ${sample_id}_R2.fastq.gz
    rm -f ${sample_id}.namesort.bam
    """
}

process EXTRACT_HLA_AND_CONVERT {
    tag "$sample_id"
    label 'process_medium'
    input:
    tuple val(sample_id), path(bam)
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: reads
    script:
    """
    set -euo pipefail
    samtools quickcheck -v ${bam}
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        samtools index -@ ${task.cpus} ${bam}
    fi
    samtools view -H ${bam} > hla_header.sam
    chr=\$(bash ${projectDir}/bin/detect_chr6_contig.sh hla_header.sam)
    samtools view -@ ${task.cpus} -b ${bam} "\${chr}:${params.hla_region_start}-${params.hla_region_end}" \
        > ${sample_id}.hla.bam
    samtools sort -n -@ ${task.cpus} ${sample_id}.hla.bam -o ${sample_id}.namesort.bam
    samtools fastq -@ ${task.cpus} -1 ${sample_id}_R1.fastq.gz -2 ${sample_id}_R2.fastq.gz \
        -0 /dev/null -s /dev/null ${sample_id}.namesort.bam
    test -s ${sample_id}_R1.fastq.gz
    test -s ${sample_id}_R2.fastq.gz
    rm -f ${sample_id}.hla.bam ${sample_id}.namesort.bam
    """
}
