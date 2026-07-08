/*
 * Immuannot Module  (assembly annotation — interpretive layer)
 *
 * Annotates immune-gene structure in ASSEMBLED sequence (contigs) against
 * IPD-IMGT/HLA gene features: assigns allele identity, resolves exon/intron
 * completeness, and flags partial/hybrid genes. Two uses in CHAMPHLA:
 *   1. Gene-structure validation of assembly-based calls (Kourami/T1K/SpecHLA,
 *      long-read assemblies) — confirms/contradicts novelty signals.
 *   2. Produces the haplotype->IMGT allele crosswalk that backs an HPRC-derived
 *      Locityper DB (consumed by parse_locityper_results.py --crosswalk).
 *
 * Requires assembled contigs as input — NOT raw reads. Most ensemble members
 * emit allele calls rather than contigs, so this module is opt-in and wired to
 * an explicit contigs channel; it is the heaviest, lowest-priority integration
 * (see plan). Run on long-read assemblies for best results.
 *
 * Prerequisites (params):
 *   params.immuannot_refdata (required): Immuannot reference bundle directory
 *   params.immuannot_dir     (optional): Immuannot install dir (for Immuannot.sh)
 */

process IMMUANNOT {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/immuannot", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_immuannot.txt"), emit: results
    path("${sample_id}_immuannot.crosswalk.tsv"), emit: crosswalk, optional: true
    path("${sample_id}.immuannot.gtf.gz"), emit: gtf, optional: true
    path "versions.yml", emit: versions

    script:
    def refdata = params.immuannot_refdata
    def imdir   = params.immuannot_dir ? "${params.immuannot_dir}/" : ''
    """
    echo "[Immuannot] Annotating contigs for ${sample_id}..."

    if [ -z "${refdata}" ] || [ ! -e "${refdata}" ]; then
        echo "[Immuannot] ERROR: params.immuannot_refdata not set or missing: '${refdata}'"
        printf "# Immuannot results for ${sample_id}\n# WARNING: refdata missing\nGene\tAllele1\tAllele2\n" > ${sample_id}_immuannot.txt
        echo '"${task.process}": {immuannot: "missing-refdata"}' > versions.yml
        exit 0
    fi

    # Immuannot aligns contigs to IPD gene features (minimap2) and emits a GTF.
    bash ${imdir}Immuannot.sh \
        -c ${contigs} \
        -r ${refdata} \
        -o ${sample_id}.immuannot \
        -t ${task.cpus} || true

    if [ -f "${sample_id}.immuannot.gtf.gz" ]; then
        python3 ${projectDir}/bin/parse_immuannot_gtf.py \
            --input  ${sample_id}.immuannot.gtf.gz \
            --sample ${sample_id} \
            --output ${sample_id}_immuannot.txt \
            --crosswalk-output ${sample_id}_immuannot.crosswalk.tsv
    else
        printf "# Immuannot results for ${sample_id}\n# WARNING: Immuannot produced no GTF\nGene\tAllele1\tAllele2\n" > ${sample_id}_immuannot.txt
    fi

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    immuannot: "\$(${imdir}Immuannot.sh --version 2>/dev/null | head -n1 || echo unknown)"
	END_VERSIONS
    """

    stub:
    """
    printf '# Immuannot results for ${sample_id} (STUB)\nGene\tAllele1\tAllele2\nA\tA*02:01\tA*11:01\nB\tB*07:02\tB*08:01\nC\tC*03:04\tC*04:01\n' > ${sample_id}_immuannot.txt
    printf '# haplotype/contig\timgt_allele\ncontig1\tA*02:01\n' > ${sample_id}_immuannot.crosswalk.tsv
    echo '"${task.process}": {immuannot: "stub"}' > versions.yml
    """
}
