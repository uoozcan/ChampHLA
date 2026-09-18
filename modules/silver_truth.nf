/*
 * Silver-truth generation (optional).
 *
 * Collects the orthogonal genotyper outputs (Locityper + Immuannot) and builds
 * a silver-standard truth_long.tsv consumable by the benchmark harness. OFF by
 * default (params.generate_truth); the main typing/consensus scaffold is
 * unaffected when disabled.
 *
 * Agreement mode (default) needs BOTH sources; with only Locityper present the
 * generator abstains unless --allow-single-source is passed (truth_require_agreement=false).
 */

process GENERATE_SILVER_TRUTH {
    tag "silver_truth"
    label 'process_low'
    publishDir "${params.outdir}/silver_truth", mode: 'copy'

    input:
    path locityper_calls
    path locityper_conf
    path immuannot_calls

    output:
    path "truth_long.tsv",       emit: truth
    path "truth_provenance.tsv", emit: provenance

    script:
    def agreement = (params.truth_require_agreement == false) ? '--allow-single-source' : '--require-agreement'
    def minGq = params.truth_min_gq ?: 20.0
    def maxNov = params.truth_max_novelty ?: 0.05
    // Globs are passed directly: the shell expands matches, and on no-match the
    // literal pattern is handed through (generate_silver_truth.py's expand_paths
    // tolerates non-matching globs). Avoids \$(...) colliding with Groovy.
    """
    python3 ${projectDir}/bin/generate_silver_truth.py \\
        --locityper-calls *_locityper.txt \\
        --locityper-confidence *_locityper.confidence.tsv \\
        --immuannot-calls *_immuannot.txt \\
        --min-gq ${minGq} \\
        --max-novelty ${maxNov} \\
        ${agreement} \\
        --output truth_long.tsv \\
        --provenance-output truth_provenance.tsv
    """

    // The generator is pure host-python (no tool binary), so the stub runs the
    // REAL command — under -stub-run it still produces genuine truth from the
    // (stubbed) upstream Locityper/Immuannot outputs.
    stub:
    def agreement = (params.truth_require_agreement == false) ? '--allow-single-source' : '--require-agreement'
    def minGq = params.truth_min_gq ?: 20.0
    def maxNov = params.truth_max_novelty ?: 0.05
    """
    python3 ${projectDir}/bin/generate_silver_truth.py \\
        --locityper-calls *_locityper.txt \\
        --locityper-confidence *_locityper.confidence.tsv \\
        --immuannot-calls *_immuannot.txt \\
        --min-gq ${minGq} \\
        --max-novelty ${maxNov} \\
        ${agreement} \\
        --output truth_long.tsv \\
        --provenance-output truth_provenance.tsv
    """
}

workflow SILVER_TRUTH_WORKFLOW {
    take:
    ch_locityper          // tuple(sample, _locityper.txt)
    ch_locityper_conf     // path(_locityper.confidence.tsv)
    ch_immuannot          // tuple(sample, _immuannot.txt)

    main:
    ch_loc  = ch_locityper.map { sample_id, f -> f }.collect().ifEmpty([])
    ch_conf = ch_locityper_conf.collect().ifEmpty([])
    ch_imm  = ch_immuannot.map { sample_id, f -> f }.collect().ifEmpty([])

    GENERATE_SILVER_TRUTH(ch_loc, ch_conf, ch_imm)

    emit:
    truth = GENERATE_SILVER_TRUTH.out.truth
    provenance = GENERATE_SILVER_TRUTH.out.provenance
}
