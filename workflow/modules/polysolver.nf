process POLYSOLVER {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/polysolver", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_polysolver.txt"), emit: results
    path("${sample_id}_polysolver_raw/winners.hla.nofreq.txt"), emit: raw
    path "polysolver_wrapper_transform.json", emit: wrapper_audit
    path "versions.yml", emit: versions

    script:
    def build = params.reference_build == 'hg38' ? 'hg38' : 'hg19'
    """
    set -euo pipefail
    mkdir -p ${sample_id}_polysolver_raw picard_tmp
    ulimit -s unlimited
    export SAMTOOLS_DIR=/home/polysolver/binaries
    export _JAVA_OPTIONS="-Djava.io.tmpdir=\$(pwd)/picard_tmp"

    # POLYSOLVER sorts its input twice. Handed the whole exome BAM that produced a 35 GB
    # work directory and exhausted the allocation, so it is given the HLA region plus
    # unmapped reads instead, exactly as HLAHD_BAM does. The unmapped reads matter:
    # POLYSOLVER uses them to recover HLA reads that failed to map to the reference.
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        /home/polysolver/binaries/samtools index ${bam}
    fi
    /home/polysolver/binaries/samtools view -H ${bam} > hla_header.sam
    chr=\$(bash ${projectDir}/bin/detect_chr6_contig.sh hla_header.sam)
    /home/polysolver/binaries/samtools view -b -h ${bam} "\${chr}:${params.hla_region_start}-${params.hla_region_end}" > hla_region.bam
    /home/polysolver/binaries/samtools view -b -f 4 ${bam} > unmapped.bam
    /home/polysolver/binaries/samtools merge -f polysolver_input.bam hla_region.bam unmapped.bam
    rm -f hla_region.bam unmapped.bam

    /home/polysolver/binaries/samtools sort -n polysolver_input.bam ${sample_id}_namesort
    /home/polysolver/binaries/samtools fixmate ${sample_id}_namesort.bam ${sample_id}_fixmate.bam
    /home/polysolver/binaries/samtools sort ${sample_id}_fixmate.bam ${sample_id}_fixed
    /home/polysolver/binaries/samtools index ${sample_id}_fixed.bam
    # The pinned wrapper hard-codes three hg38 intervals on contig `6`. Derive a wrapper
    # for the BAM's observed convention, while separately relocating its temporary files.
    # The transformer verifies the immutable source and exact substitution counts.
    python3 ${projectDir}/bin/patch_polysolver_wrapper.py \
        --source /home/polysolver/scripts/shell_call_hla_type \
        --output patched_shell_call_hla_type \
        --audit polysolver_wrapper_transform.json \
        --spec ${projectDir}/conf/polysolver_wrapper_patch.json \
        --contig "\${chr}" --tmp-dir "\$(pwd)/picard_tmp"
    test -s polysolver_wrapper_transform.json
    python3 ${projectDir}/bin/patch_polysolver_wrapper.py --verify-only \
        --source /home/polysolver/scripts/shell_call_hla_type \
        --output patched_shell_call_hla_type \
        --audit polysolver_wrapper_transform.json \
        --spec ${projectDir}/conf/polysolver_wrapper_patch.json \
        --contig "\${chr}" --tmp-dir "\$(pwd)/picard_tmp"
    bash patched_shell_call_hla_type ${sample_id}_fixed.bam Unknown 0 ${build} STDFQ 0 \
        ${sample_id}_polysolver_raw
    native=${sample_id}_polysolver_raw/winners.hla.nofreq.txt
    test -s "\${native}"
    python3 ${projectDir}/bin/parse_polysolver_results.py --input "\${native}" \
        --sample ${sample_id} --output ${sample_id}_polysolver.txt
    test -s ${sample_id}_polysolver.txt
    printf '"%s":\\n    polysolver: "v4"\\n' "${task.process}" > versions.yml
    rm -f ${sample_id}_namesort.bam ${sample_id}_fixmate.bam \
        ${sample_id}_fixed.bam ${sample_id}_fixed.bam.bai polysolver_input.bam
    """
}
