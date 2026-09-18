/*
 * HLA-HD Module
 * Accurate HLA typing from WGS/WES data
 * Supports both BAM and FASTQ inputs
 */

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
    tuple val(sample_id), path("${sample_id}/${sample_id}/result/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    def ref_version = reference == 'hg19' ? '37' : '38'
    """
    # Create working directory
    mkdir -p ${sample_id}
    trap 'rm -f hla_region.bam unmapped.bam merged.bam sorted.bam R1.fastq R2.fastq core.*' EXIT

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index..."
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Extract HLA reads from BAM
    echo "[Step 1] Extracting HLA reads..."
    samtools view -b -h ${bam} chr6:${params.hla_region_start}-${params.hla_region_end} > hla_region.bam 2>/dev/null || \
    samtools view -b -h ${bam} 6:${params.hla_region_start}-${params.hla_region_end} > hla_region.bam

    # Add unmapped reads
    samtools view -b -f 4 ${bam} > unmapped.bam
    samtools merge -f merged.bam hla_region.bam unmapped.bam
    samtools sort -n -@ ${task.cpus} merged.bam -o sorted.bam

    # Convert to FASTQ
    samtools fastq -@ ${task.cpus} -1 R1.fastq -2 R2.fastq -0 /dev/null -s /dev/null sorted.bam

    # Run HLA-HD
    echo "[Step 2] Running HLA-HD..."
    HLAHD_BIN=\$(command -v hlahd.sh 2>/dev/null || true)
    [ -n "\$HLAHD_BIN" ] || HLAHD_BIN=/app/hlahd.1.4.0/bin/hlahd.sh
    [ -x "\$HLAHD_BIN" ] || { echo "HLA-HD executable not found" >&2; exit 127; }
    # Detect read length and set -m accordingly (stfr -L TSIZE filters reads shorter than TSIZE)
    READ_LEN=\$(awk 'NR==2{print length(\$0); exit}' R1.fastq 2>/dev/null || echo 100)
    MIN_TAG=\$(( READ_LEN < 100 ? READ_LEN / 2 : 50 ))
    echo "[Read length detected: \${READ_LEN}bp, using -m \${MIN_TAG}]"
    "\$HLAHD_BIN" -t ${task.cpus} -m \${MIN_TAG} -c 0.95 -f ${params.hlahd_db}/freq_data \
        R1.fastq R2.fastq ${params.hlahd_db}/HLA_gene.split.txt ${params.hlahd_db}/dictionary \
        ${sample_id} ${sample_id}

    # Parse results with per-gene read counts
    # Use cut to extract exactly 3 fields (avoids trailing fields for HLA-E etc.)
    echo "[Step 3] Parsing results with read counts...]"

    RESULT_DIR="${sample_id}/${sample_id}/result"
    if [ -f "\${RESULT_DIR}/${sample_id}_final.result.txt" ]; then
        while IFS='' read -r LINE; do
            GENE=\$(printf '%s\n' "\$LINE" | cut -f1)
            ALLELE1=\$(printf '%s\n' "\$LINE" | cut -f2)
            ALLELE2=\$(printf '%s\n' "\$LINE" | cut -f3)
            GENE_SHORT="\${GENE#HLA-}"
            READ_FILE="\${RESULT_DIR}/${sample_id}_\${GENE_SHORT}.read.txt"
            READS=0
            if [ -f "\${READ_FILE}" ]; then
                RAW=\$(awk 'NR==1{print \$2}' "\${READ_FILE}" 2>/dev/null)
                [[ "\${RAW}" =~ ^[0-9]+\$ ]] && READS=\${RAW}
            fi
            printf '%s\t%s\t%s\t%s\t%s\n' "\${GENE}" "\${ALLELE1}" "\${ALLELE2}" "\${READS}" "\${READS}"
        done < "\${RESULT_DIR}/${sample_id}_final.result.txt" > ${sample_id}_hlahd.txt
    else
        echo "# HLA-HD results for ${sample_id}" > ${sample_id}_hlahd.txt
        echo "# No results generated" >> ${sample_id}_hlahd.txt
    fi

    # Version info
    HLAHD_VER=\$([ -n "\$HLAHD_BIN" ] && "\$HLAHD_BIN" 2>&1 | grep -i version | head -1 || echo "1.4.0")
    SAM_VER=\$(samtools --version | head -1 | cut -d' ' -f2)
    printf '"${task.process}":\n    hlahd: %s\n    samtools: %s\n' "\$HLAHD_VER" "\$SAM_VER" > versions.yml
    """
}

/*
 * HLA-HD from paired FASTQ files
 */
process HLAHD_FASTQ {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/hlahd", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)
    val hla_genes

    output:
    tuple val(sample_id), path("${sample_id}_hlahd.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/${sample_id}/result/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create working directory
    mkdir -p ${sample_id}
    trap 'rm -f R1.fastq R2.fastq R1_pre.fastq R2_pre.fastq hla_pre.sam core.*' EXIT

    # Prepare FASTQ files (decompress if needed)
    if [[ "${fastq1}" == *.gz ]]; then
        zcat ${fastq1} > R1_pre.fastq
        zcat ${fastq2} > R2_pre.fastq
    else
        cp ${fastq1} R1_pre.fastq
        cp ${fastq2} R2_pre.fastq
    fi

    # Pre-filter: align against HLA-HD's own Bowtie2 index to extract HLA-enriched reads.
    # Without this, HLA-HD's internal FASTQ splitting fails silently on large RNA-seq inputs
    # (mapfile/*.fastq end up 0 bytes → "Unable to read file magic number" for all genes).
    echo "[Pre-filtering HLA reads using HLA-HD dictionary Bowtie2 index...]"
    HLA_BT2_IDX="${params.hlahd_db}/dictionary/all_exon_intron_N150.fasta"
    bowtie2 -p ${task.cpus} --no-unal -x "\${HLA_BT2_IDX}" \
        -1 R1_pre.fastq -2 R2_pre.fastq -S hla_pre.sam 2>/dev/null || true

    MAPPED_READS=\$(grep -c -v "^@" hla_pre.sam 2>/dev/null || echo 0)
    if [ -s hla_pre.sam ] && [ "\${MAPPED_READS}" -gt 100 ]; then
        echo "[Pre-filter: \${MAPPED_READS} HLA-mapped reads → using enriched subset]"
        samtools sort -n hla_pre.sam | samtools fastq -1 R1.fastq -2 R2.fastq -0 /dev/null -s /dev/null -
    else
        echo "[Pre-filter: too few mapped reads (\${MAPPED_READS}), using full input]"
        mv R1_pre.fastq R1.fastq
        mv R2_pre.fastq R2.fastq
    fi
    rm -f hla_pre.sam R1_pre.fastq R2_pre.fastq

    # Run HLA-HD
    echo "[Running HLA-HD from FASTQ...]"
    HLAHD_BIN=\$(command -v hlahd.sh 2>/dev/null || true)
    [ -n "\$HLAHD_BIN" ] || HLAHD_BIN=/app/hlahd.1.4.0/bin/hlahd.sh
    [ -x "\$HLAHD_BIN" ] || { echo "HLA-HD executable not found" >&2; exit 127; }
    # Detect read length and set -m accordingly (stfr -L TSIZE filters reads shorter than TSIZE)
    READ_LEN=\$(awk 'NR==2{print length(\$0); exit}' R1.fastq 2>/dev/null || echo 100)
    MIN_TAG=\$(( READ_LEN < 100 ? READ_LEN / 2 : 50 ))
    echo "[Read length detected: \${READ_LEN}bp, using -m \${MIN_TAG}]"
    "\$HLAHD_BIN" -t ${task.cpus} -m \${MIN_TAG} -c 0.95 -f ${params.hlahd_db}/freq_data \
        R1.fastq R2.fastq ${params.hlahd_db}/HLA_gene.split.txt ${params.hlahd_db}/dictionary \
        ${sample_id} ${sample_id}

    # Parse results with per-gene read counts
    echo "[Parsing results with read counts...]"
    RESULT_DIR="${sample_id}/${sample_id}/result"
    if [ -f "\${RESULT_DIR}/${sample_id}_final.result.txt" ]; then
        while IFS='' read -r LINE; do
            GENE=\$(printf '%s\n' "\$LINE" | cut -f1)
            ALLELE1=\$(printf '%s\n' "\$LINE" | cut -f2)
            ALLELE2=\$(printf '%s\n' "\$LINE" | cut -f3)
            GENE_SHORT="\${GENE#HLA-}"
            READ_FILE="\${RESULT_DIR}/${sample_id}_\${GENE_SHORT}.read.txt"
            READS=0
            if [ -f "\${READ_FILE}" ]; then
                RAW=\$(awk 'NR==1{print \$2}' "\${READ_FILE}" 2>/dev/null)
                [[ "\${RAW}" =~ ^[0-9]+\$ ]] && READS=\${RAW}
            fi
            printf '%s\t%s\t%s\t%s\t%s\n' "\${GENE}" "\${ALLELE1}" "\${ALLELE2}" "\${READS}" "\${READS}"
        done < "\${RESULT_DIR}/${sample_id}_final.result.txt" > ${sample_id}_hlahd.txt
    else
        echo "# HLA-HD results for ${sample_id}" > ${sample_id}_hlahd.txt
        echo "# No results generated" >> ${sample_id}_hlahd.txt
    fi

    # Cleanup
    rm -f R1.fastq R2.fastq

    # Version info
    HLAHD_VER=\$([ -n "\$HLAHD_BIN" ] && "\$HLAHD_BIN" 2>&1 | grep -i version | head -1 || echo "1.4.0")
    printf '"${task.process}":\n    hlahd: %s\n' "\$HLAHD_VER" > versions.yml
    """

    stub:
    """
    printf '# HLA-HD results for ${sample_id} (STUB)\nGene\tAllele1\tAllele2\tReads1\tReads2\nA\tA*02:01\tA*11:01\t50\t50\nB\tB*07:02\tB*08:01\t50\t50\nC\tC*03:04\tC*04:01\t50\t50\n' > ${sample_id}_hlahd.txt
    printf '"${task.process}":\n    hlahd: stub\n' > versions.yml
    """
}
