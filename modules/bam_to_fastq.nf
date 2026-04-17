/*
 * Shared preprocessing utilities for BAM/CRAM inputs.
 */

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
    samtools view -@ ${task.cpus} -b -T ${reference_fasta} -o ${sample_id}.bam ${cram}
    samtools index -@ ${task.cpus} ${sample_id}.bam
    """
}

process EXTRACT_HLA_REGION {
    tag "$sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.hla.bam"), emit: bam

    script:
    """
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
      samtools index -@ ${task.cpus} ${bam}
    fi
    CHR=\$(samtools view -H ${bam} | awk '/^@SQ.*SN:chr6\t/{print "chr6"; exit} /^@SQ.*SN:6\t/{print "6"; exit}')
    if [ -z "\$CHR" ]; then
      CHR=6
    fi
    samtools view -@ ${task.cpus} -b ${bam} "\${CHR}:${params.hla_region_start}-${params.hla_region_end}" > ${sample_id}.hla.bam
    samtools index -@ ${task.cpus} ${sample_id}.hla.bam
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
    samtools sort -n -@ ${task.cpus} ${bam} -o ${sample_id}.namesort.bam
    samtools fastq -@ ${task.cpus} \
      -1 ${sample_id}_R1.fastq.gz \
      -2 ${sample_id}_R2.fastq.gz \
      -0 /dev/null -s /dev/null \
      ${sample_id}.namesort.bam
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
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
      samtools index -@ ${task.cpus} ${bam}
    fi
    CHR=\$(samtools view -H ${bam} | awk '/^@SQ.*SN:chr6\t/{print "chr6"; exit} /^@SQ.*SN:6\t/{print "6"; exit}')
    if [ -z "\$CHR" ]; then
      CHR=6
    fi
    samtools view -@ ${task.cpus} -b ${bam} "\${CHR}:${params.hla_region_start}-${params.hla_region_end}" > ${sample_id}.hla.bam
    samtools sort -n -@ ${task.cpus} ${sample_id}.hla.bam -o ${sample_id}.hla.namesort.bam
    samtools fastq -@ ${task.cpus} \
      -1 ${sample_id}_R1.fastq.gz \
      -2 ${sample_id}_R2.fastq.gz \
      -0 /dev/null -s /dev/null \
      ${sample_id}.hla.namesort.bam
    rm -f ${sample_id}.hla.bam ${sample_id}.hla.namesort.bam
    """
}
