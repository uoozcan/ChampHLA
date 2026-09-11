process POLYSOLVER {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/polysolver", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_polysolver.txt"), emit: results
    path("${sample_id}_polysolver_raw/winners.hla.nofreq.txt"), emit: raw
    path "versions.yml", emit: versions

    script:
    def build = params.reference_build == 'hg38' ? 'hg38' : 'hg19'
    """
    set -euo pipefail
    mkdir -p ${sample_id}_polysolver_raw picard_tmp
    ulimit -s unlimited
    export SAMTOOLS_DIR=/home/polysolver/binaries
    export _JAVA_OPTIONS="-Djava.io.tmpdir=\$(pwd)/picard_tmp"
    /home/polysolver/binaries/samtools sort -n ${bam} ${sample_id}_namesort
    /home/polysolver/binaries/samtools fixmate ${sample_id}_namesort.bam ${sample_id}_fixmate.bam
    /home/polysolver/binaries/samtools sort ${sample_id}_fixmate.bam ${sample_id}_fixed
    /home/polysolver/binaries/samtools index ${sample_id}_fixed.bam
    sed "s|TMP_DIR=/home/polysolver|TMP_DIR=\$(pwd)/picard_tmp|g" \
        /home/polysolver/scripts/shell_call_hla_type > patched_shell_call_hla_type
    chmod +x patched_shell_call_hla_type
    bash patched_shell_call_hla_type ${sample_id}_fixed.bam Unknown 0 ${build} STDFQ 0 \
        ${sample_id}_polysolver_raw
    native=${sample_id}_polysolver_raw/winners.hla.nofreq.txt
    test -s "\${native}"
    python3 ${projectDir}/bin/parse_polysolver_results.py --input "\${native}" \
        --sample ${sample_id} --output ${sample_id}_polysolver.txt
    test -s ${sample_id}_polysolver.txt
    printf '"%s":\n    polysolver: "v4"\n' "${task.process}" > versions.yml
    rm -f ${sample_id}_namesort.bam ${sample_id}_fixmate.bam \
        ${sample_id}_fixed.bam ${sample_id}_fixed.bam.bai
    """
}
