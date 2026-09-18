/*
 * POLYSOLVER Module
 * HLA Class I typing optimised for tumor/normal BAMs
 * Broad Institute tool widely used in cancer genomics (AML and somatic studies)
 * BAM input only — hg19/GRCh37 OR hg38 aligned coordinate-sorted BAM
 *
 * Container: docker://sachet/polysolver:v4  (pull: singularity pull polysolver.sif docker://sachet/polysolver:v4)
 *
 * Fixes applied (effective on native Linux / Puhti; novoalign does NOT work on WSL2):
 *   1. novoalign SIGSEGV: `ulimit -s unlimited` added; helps on some systems with small default stacks
 *      NOTE: novoalign (v2 and v3) crashes with SIGSEGV on WSL2 kernel 6.6 regardless of ulimit —
 *      WSL2-specific kernel incompatibility; works on native Linux (CentOS 7, Ubuntu 18.04+) and Puhti
 *   2. hg38 SAMTOOLS_DIR: env var unset in container → `export SAMTOOLS_DIR=/home/polysolver/binaries`
 *   3. hg38 Picard "Illegal mate state": fixmate pre-processing fixes inconsistent mate flags
 */

process POLYSOLVER {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/polysolver", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_polysolver.txt"), emit: results
    path("${sample_id}_polysolver_raw/winners.hla.nofreq.txt"), emit: raw, optional: true
    path "versions.yml", emit: versions

    script:
    def build = params.reference_build == 'hg38' ? 'hg38' : 'hg19'
    """
    mkdir -p ${sample_id}_polysolver_raw

    echo "[POLYSOLVER] Running on ${sample_id} (build: ${build})..."

    # Fix 1: novoalign stack overflow on modern kernels (default 8MB stack too small)
    ulimit -s unlimited

    # Fix 2: SAMTOOLS_DIR env var is unset in container for hg38 branch
    export SAMTOOLS_DIR=/home/polysolver/binaries

    # Fix 3: Picard SamToFastq "Illegal mate state" — pre-sort by name + fixmate + re-sort
    # Required when BAM has inconsistent mate flags (common in some pipelines)
    echo "[POLYSOLVER] Applying fixmate pre-processing for ${sample_id}..."
    /home/polysolver/binaries/samtools sort -n ${bam} ${sample_id}_namesort
    /home/polysolver/binaries/samtools fixmate ${sample_id}_namesort.bam ${sample_id}_fixmate.bam
    /home/polysolver/binaries/samtools sort ${sample_id}_fixmate.bam ${sample_id}_fixed
    /home/polysolver/binaries/samtools index ${sample_id}_fixed.bam
    POLYSOLVER_INPUT="${sample_id}_fixed.bam"

    # Fix 4: SortSam.jar has TMP_DIR=/home/polysolver hardcoded — patch script at runtime
    # /home/polysolver is read-only in Singularity; replace with writable CWD tmp
    mkdir -p picard_tmp
    PICARD_TMP=\$(pwd)/picard_tmp
    export _JAVA_OPTIONS="-Djava.io.tmpdir=\${PICARD_TMP}"
    # Fix 5: POLYSOLVER hardcodes its hg38 HLA regions without a chr prefix --
    # `samtools view $bam 6:29941260-29945884` -- but GRCh38_full_analysis_set BAMs
    # name the contig chr6. Every region query then matches nothing, POLYSOLVER
    # writes no winners file, and the run looks like a tool that found no alleles:
    #   [main_samview] region "6:29941260-29945884" specifies an unknown reference name
    # Detect the naming from the BAM header, the same way modules/spechla.nf does,
    # and rewrite the queries to match.
    CHR6=\$(/home/polysolver/binaries/samtools view -H ${bam} \
              | awk '/^@SQ.*SN:chr6\t/{print "chr6"; exit} /^@SQ.*SN:6\t/{print "6"; exit}')
    if [ -z "\${CHR6}" ]; then CHR6=6; fi
    echo "[POLYSOLVER] BAM names chromosome 6 as '\${CHR6}'"

    if [ "\${CHR6}" = "chr6" ]; then
        SED_CHR='s|\\\$bam 6:|\\\$bam chr6:|g'
    else
        SED_CHR='s|__panelhla_noop__|__panelhla_noop__|'
    fi

    sed -e "s|TMP_DIR=/home/polysolver|TMP_DIR=\${PICARD_TMP}|g" \
        -e "\${SED_CHR}" \
        /home/polysolver/scripts/shell_call_hla_type > patched_shell_call_hla_type
    chmod +x patched_shell_call_hla_type

    # POLYSOLVER args: BAM race includeFreq build format insertCalc outdir
    # race=Unknown (population-agnostic), includeFreq=0, insertCalc=0 (germline)
    # Record the status instead of discarding it with `|| true`. The process sets
    # errorStrategy 'ignore', so one bad sample still does not kill a run -- but a
    # crash is now distinguishable from a clean run that found nothing.
    set +e
    bash patched_shell_call_hla_type \
        \$POLYSOLVER_INPUT Unknown 0 ${build} STDFQ 0 ${sample_id}_polysolver_raw
    POLYSOLVER_RC=\$?
    set -e
    echo "[POLYSOLVER] shell_call_hla_type exit=\${POLYSOLVER_RC}"

    # Cleanup intermediate BAMs
    rm -f ${sample_id}_namesort.bam ${sample_id}_fixmate.bam ${sample_id}_fixed.bam ${sample_id}_fixed.bam.bai

    # Parse winners.hla.nofreq.txt → standard pipeline TSV
    if [ -f "${sample_id}_polysolver_raw/winners.hla.nofreq.txt" ]; then
        python3 ${projectDir}/bin/parse_polysolver_results.py \
            --input  ${sample_id}_polysolver_raw/winners.hla.nofreq.txt \
            --sample ${sample_id} \
            --output ${sample_id}_polysolver.txt
    else
        # Never manufacture an empty result. A placeholder file exits 0 and is
        # indistinguishable downstream from a caller that genuinely typed nothing.
        echo "POLYSOLVER produced no winners.hla.nofreq.txt for ${sample_id}" >&2
        echo "shell_call_hla_type exited \${POLYSOLVER_RC}; chromosome 6 named '\${CHR6}'." >&2
        echo "Contents of ${sample_id}_polysolver_raw:" >&2
        ls -la ${sample_id}_polysolver_raw >&2 || true
        exit 1
    fi

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    polysolver: "v4"
END_VERSIONS
    """
}
