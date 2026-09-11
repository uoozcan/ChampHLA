process KOURAMI {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/kourami", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_kourami.txt"), emit: results
    path("${sample_id}.kourami.result"), emit: raw
    path "versions.yml", emit: versions

    script:
    def kourami_dir = params.kourami_dir ?: '/usr/local/bin/kourami-0.9.6'
    def kourami_db = params.kourami_db ?: '/usr/local/bin/kourami-0.9.6/db'
    """
    set -euo pipefail
    test -d ${kourami_db}
    samtools quickcheck -v ${bam}
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        samtools index ${bam}
    fi
    samtools view -H ${bam} > hla_header.sam
    grep -E '^@SQ.*SN:(chr6|6)[[:space:]]' hla_header.sam > /dev/null
    chr=\$(awk '/^@SQ.*SN:chr6[[:space:]]/{print "chr6"; exit} /^@SQ.*SN:6[[:space:]]/{print "6"; exit}' hla_header.sam)
    samtools view -b ${bam} "\${chr}:${params.hla_region_start}-${params.hla_region_end}" | \
        samtools sort -n -@ ${task.cpus} | \
        samtools fastq -1 ${sample_id}._hla_1.fq.gz -2 ${sample_id}._hla_2.fq.gz -s /dev/null -
    test -s ${sample_id}._hla_1.fq.gz
    test -s ${sample_id}._hla_2.fq.gz
    bwa mem -t ${task.cpus} ${kourami_db}/All_FINAL_with_Decoy.fa.gz \
        ${sample_id}._hla_1.fq.gz ${sample_id}._hla_2.fq.gz -o ${sample_id}.panel.sam
    samtools sort -@ ${task.cpus} -o ${sample_id}.panel.bam ${sample_id}.panel.sam
    samtools index ${sample_id}.panel.bam
    jar=${kourami_dir}/target/Kourami.jar
    test -s "\${jar}"
    java -Xmx10g -jar "\${jar}" -d ${kourami_db} ${sample_id}.panel.bam -o ${sample_id}.kourami
    test -s ${sample_id}.kourami.result
    python3 ${projectDir}/bin/parse_kourami_results.py --input ${sample_id}.kourami.result \
        --sample ${sample_id} --output ${sample_id}_kourami.txt
    test -s ${sample_id}_kourami.txt
    printf '"%s":\\n    kourami: "0.9.6"\\n' "${task.process}" > versions.yml
    rm -f ${sample_id}.panel.sam ${sample_id}.panel.bam ${sample_id}.panel.bam.bai \
        ${sample_id}._hla_1.fq.gz ${sample_id}._hla_2.fq.gz
    """
}
