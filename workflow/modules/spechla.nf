process SPECHLA {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/spechla", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    val reference

    output:
    tuple val(sample_id), path("${sample_id}_spechla.txt"), emit: results
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    # Guarded with -r, not -x: the wrapper is invoked as `bash <script>` below, which
    # needs only read permission. The deployed wrapper is mode 660, so an -x guard
    # aborted the process under `set -e` with no output at all.
    test -r ${params.spechla_path}/script/whole/SpecHLA.sh
    mkdir -p ${sample_id}
    samtools quickcheck -v ${bam}
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        samtools index ${bam}
    fi
    samtools view -H ${bam} > hla_header.sam
    grep -E '^@SQ.*SN:(chr6|6)[[:space:]]' hla_header.sam > /dev/null
    chr=\$(awk '/^@SQ.*SN:chr6[[:space:]]/{print "chr6"; exit} /^@SQ.*SN:6[[:space:]]/{print "6"; exit}' hla_header.sam)
    samtools view -b ${bam} "\${chr}:${params.hla_region_start}-${params.hla_region_end}" > ${sample_id}/hla_region.bam
    samtools sort -n ${sample_id}/hla_region.bam -o ${sample_id}/namesort.bam
    samtools fastq -1 ${sample_id}/R1.fastq.gz -2 ${sample_id}/R2.fastq.gz \
        -0 /dev/null -s /dev/null ${sample_id}/namesort.bam
    test -s ${sample_id}/R1.fastq.gz
    test -s ${sample_id}/R2.fastq.gz
    cd ${sample_id}
    bash ${params.spechla_path}/script/whole/SpecHLA.sh -n ${sample_id} \
        -1 R1.fastq.gz -2 R2.fastq.gz -o . -j ${task.cpus} \
        -u ${params.spechla_exon_only}
    cd ..
    if [ -s "${sample_id}/hla.result.txt" ]; then
        cp ${sample_id}/hla.result.txt ${sample_id}_spechla.txt
    elif [ -s "${sample_id}/${sample_id}/hla.result.txt" ]; then
        cp ${sample_id}/${sample_id}/hla.result.txt ${sample_id}_spechla.txt
    else
        echo "SpecHLA produced no native result" >&2
        exit 4
    fi
    printf '"%s":\\n    spechla: "1.0.7-deployed-wrapper"\\n' "${task.process}" > versions.yml
    rm -f ${sample_id}/hla_region.bam ${sample_id}/namesort.bam \
        ${sample_id}/R1.fastq.gz ${sample_id}/R2.fastq.gz
    """
}
