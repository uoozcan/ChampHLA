/*
 * Locityper Module
 * Locus-pangenome HLA genotyping (sequence-level), Class I/II, WGS + long-read.
 *
 * Locityper recruits reads to each locus and selects the haplotype pair that
 * best explains the alignment AND the read-depth profile against a locus
 * pangenome DB (see bin/build_locityper_db.sh). It contributes an INDEPENDENT,
 * depth-aware track distinct from the alignment/assembly tools already in the
 * pipeline, plus a genotype-quality (GQ) confidence and a novel-allele signal.
 *
 * BAM/CRAM input — coordinate-sorted, indexed. Intended for WGS and long-read
 * only (full intronic coverage); gated off for WES/RNA upstream in main.nf.
 *
 * Prerequisites (params):
 *   params.locityper_db    (required): path to the DB built by build_locityper_db.sh
 *   params.reference_fasta (required): reference the BAM was aligned to
 *   params.seq_type        : dna | longreads_hifi | longreads_ont (selects preset)
 *   params.locityper_loci  (optional): comma-separated locus allow-list
 */

process LOCITYPER {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/locityper", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_locityper.txt"), emit: results
    path("${sample_id}_locityper.confidence.tsv"), emit: confidence, optional: true
    path("${sample_id}.locityper_out"), emit: raw, optional: true
    path "versions.yml", emit: versions

    script:
    def db        = params.locityper_db
    def ref       = params.reference_fasta
    def seq_type  = (params.seq_type ?: 'dna').toLowerCase()
    def loci_arg  = params.locityper_loci ? "--loci ${params.locityper_loci}" : ''
    // Read technology preset: Locityper models long reads differently from short.
    def tech = (seq_type == 'longreads_hifi') ? 'hifi'
             : (seq_type == 'longreads_ont')  ? 'ont'
             : 'illumina'
    """
    echo "[Locityper] Running on ${sample_id} (tech=${tech})..."

    if [ -z "${db}" ] || [ ! -e "${db}" ]; then
        echo "[Locityper] ERROR: params.locityper_db not set or missing: '${db}'"
        printf "# Locityper results for ${sample_id}\n# WARNING: DB missing\nGene\tAllele1\tAllele2\n" > ${sample_id}_locityper.txt
        echo '"${task.process}": {locityper: "missing-db"}' > versions.yml
        exit 0
    fi

    [ -f "${bam}.bai" ] || [ -f "${bam}.crai" ] || samtools index ${bam}

    # Step 1: preprocess — estimate background read depth for the sample.
    locityper preproc \
        --input ${bam} \
        --reference ${ref} \
        --technology ${tech} \
        --threads ${task.cpus} \
        --output ${sample_id}.locityper_preproc || true

    # Step 2: genotype the target loci against the locus-pangenome DB.
    locityper genotype \
        --input ${bam} \
        --database ${db} \
        --preproc ${sample_id}.locityper_preproc \
        --reference ${ref} \
        --threads ${task.cpus} \
        --output ${sample_id}.locityper_out || true

    # Step 3: parse to standard pipeline TSV + confidence/novelty side-channel.
    if [ -d "${sample_id}.locityper_out" ]; then
        python3 ${projectDir}/bin/parse_locityper_results.py \
            --input  ${sample_id}.locityper_out \
            --sample ${sample_id} \
            --output ${sample_id}_locityper.txt \
            --confidence-output ${sample_id}_locityper.confidence.tsv \
            ${loci_arg}
    else
        printf "# Locityper results for ${sample_id}\n# WARNING: Locityper produced no output\nGene\tAllele1\tAllele2\n" > ${sample_id}_locityper.txt
    fi

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    locityper: "\$(locityper --version 2>/dev/null | head -n1 | sed 's/^[^0-9]*//' || echo unknown)"
	END_VERSIONS
    """

    stub:
    """
    printf '# Locityper results for ${sample_id} (STUB)\nGene\tAllele1\tAllele2\nA\tA*02:01\tA*11:01\nB\tB*07:02\tB*08:01\nC\tC*03:04\tC*04:01\n' > ${sample_id}_locityper.txt
    printf 'sample\tgene\tconfidence_score\tread_support\tgq\tweighted_dist\tmean_depth\tnovelty_score\n${sample_id}\tA\t0.9999\t30\t40\t3\t30\t0.0001\n${sample_id}\tB\t0.9999\t30\t38\t4\t30\t0.0002\n${sample_id}\tC\t0.9999\t30\t37\t4\t30\t0.0002\n' > ${sample_id}_locityper.confidence.tsv
    echo '"${task.process}": {locityper: "stub"}' > versions.yml
    """
}
