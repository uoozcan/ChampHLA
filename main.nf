#!/usr/bin/env nextflow
/*
========================================================================================
    HLA TYPING MULTI-TOOL PIPELINE
========================================================================================
    Github : https://github.com/yourusername/hla-typing-pipeline
    Version: 2.0.0
========================================================================================
    Comprehensive HLA typing pipeline integrating multiple state-of-the-art tools:
    - OptiType: DNA/RNA HLA typing with high precision
    - ArcasHLA: Fast HLA typing from RNA-seq/DNA data
    - SpecHLA: Exome/genome HLA typing with variant calling
========================================================================================
*/

nextflow.enable.dsl = 2

// Import modules
include { OPTITYPE_FASTQ as OPTITYPE } from './modules/optitype'
include { ARCASHLA as ARCASHLA_BAM; ARCASHLA_FASTQ } from './modules/arcashla'
include { SPECHLA as SPECHLA_BAM; SPECHLA_FASTQ } from './modules/spechla'
include { HLAHD as HLAHD_BAM; HLAHD_FASTQ } from './modules/hlahd'
include { POLYSOLVER } from './modules/polysolver'
include { KOURAMI } from './modules/kourami'
include { T1K_FASTQ } from './modules/t1k'
include { SEQ2HLA } from './modules/seq2hla'
include { CRAM_TO_BAM; EXTRACT_HLA_REGION; BAM_TO_FASTQ; EXTRACT_HLA_AND_CONVERT } from './modules/bam_to_fastq'
include { AGGREGATE_RESULTS } from './modules/aggregation'
include { MAJORITY_VOTING_WORKFLOW } from './modules/majority_voting'

/*
========================================================================================
    PARAMETER VALIDATION
========================================================================================
*/

def helpMessage() {
    log.info"""
    ============================================================
    HLA TYPING MULTI-TOOL PIPELINE v2.0.0
    ============================================================
    
    Usage:
    nextflow run main.nf --input <path> [options]
    
    Required Arguments (choose one input method):
      --input              Path to directory of input files
      --input_type         Input type: 'bam', 'cram', or 'fastq' (required with --input)
      --input_samplesheet  Path to CSV samplesheet (alternative to --input)
                           BAM format:   sample_id,bam_path
                           FASTQ format: sample_id,fastq_1,fastq_2

    Pipeline Options:
      --tools              Comma-separated list of tools to run
                           Available: optitype,arcashla,spechla,hlahd,polysolver,kourami,t1k,seq2hla
                           Default: optitype,arcashla,hlahd
      --seq_type           Sequence type: dna, rna, longreads_hifi, longreads_ont [default: dna]
                           Automatically skips incompatible tools (e.g. polysolver/kourami for RNA)
      --outdir             Output directory [default: ./results]
      --reference_build    Reference genome build: hg19 or hg38 [default: hg38]
      --reference_fasta    Reference FASTA (required for CRAM input decoding)

    Tool-specific Options:
      --optitype_seq_type  OptiType sequence type: 'dna' or 'rna' [default: dna]
      --arcashla_genes     Genes for ArcasHLA [default: A,B,C,DQA1,DQB1,DRB1]
      --spechla_genes      Genes for SpecHLA [default: A,B,C,DQA1,DQB1,DRB1]
      --spechla_exon_only  SpecHLA exon-only mode (0 or 1) [default: 0]
                           Use 1 for exome/WES data!

    Consensus Options:
      --enable_majority_voting  Enable consensus calling [default: false]
      --weighting          Weighting strategy: 'equal' or 'calibrated' [default: calibrated]
                           'calibrated' auto-loads pre-computed weights for wgs/wes/rnaseq
                           Note: Consensus defaults to A/B/C focus; detailed per-tool outputs still contain all reported genes.
      --consensus_weight_file  Override: path to custom JSON weight file
      --mv_min_tools       Minimum tools for consensus [default: 2]
      --mv_resolution      Resolution level: 2 or 4 [default: 4]
      --mv_genes           Genes for voting [default: A,B,C]

    HPC Options:
      --slurm_account      SLURM account for job submission
      --slurm_partition    SLURM partition [default: small]
      --max_cpus           Maximum CPUs per job [default: 40]
      --max_memory         Maximum memory per job [default: 180.GB]
      --max_time           Maximum time per job [default: 24.h]

    Other Options:
      --save_intermediate  Save intermediate files [default: false]
      --extract_hla_region Extract HLA region before processing [default: false]
      --help               Show this help message

    Examples:
      # RNA-seq with samplesheet (auto-skips DNA-only tools)
      nextflow run main.nf \\
        --input_samplesheet samples_rna.csv \\
        --seq_type rna \\
        --tools optitype,arcashla,seq2hla \\
        --optitype_seq_type rna \\
        -profile singularity

      # Exome BAMs with SpecHLA (exon-only) and calibrated consensus
      nextflow run main.nf \\
        --input bam_files/ \\
        --input_type bam \\
        --seq_type wes \\
        --tools spechla,hlahd,optitype \\
        --spechla_exon_only 1 \\
        --enable_majority_voting \\
        -profile puhti,singularity

      # WGS CRAM files with full tool set
      nextflow run main.nf \\
        --input cram_files/ \\
        --input_type cram \\
        --reference_fasta /path/to/hg38.fa \\
        --tools optitype,arcashla,spechla,hlahd,kourami,t1k \\
        --enable_majority_voting \\
        -profile singularity -resume
    
    ============================================================
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.input && !params.input_samplesheet) {
    log.error "ERROR: Either --input (with --input_type) or --input_samplesheet is required!"
    helpMessage()
    exit 1
}

if (params.input && !params.input_type) {
    log.error "ERROR: --input_type is required when using --input (bam, cram, or fastq)!"
    helpMessage()
    exit 1
}

if (params.input_type == 'cram' && !params.reference_fasta) {
    log.error "ERROR: --reference_fasta is required for CRAM input!"
    exit 1
}

// Print parameter summary
log.info """
========================================================================================
    HLA TYPING MULTI-TOOL PIPELINE
========================================================================================
Input               : ${params.input_samplesheet ?: params.input}
Input type          : ${params.input_samplesheet ? 'samplesheet' : params.input_type}
Seq type            : ${params.seq_type ?: 'dna'}
Output dir          : ${params.outdir}
Tools               : ${params.tools}
Reference build     : ${params.reference_build}
Weighting           : ${params.weighting ?: 'calibrated'}
Max CPUs            : ${params.max_cpus}
Max memory          : ${params.max_memory}
Max time            : ${params.max_time}
========================================================================================
""".stripIndent()

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    // Parse tools to run
    def tools_list = params.tools.tokenize(',')*.trim()*.toLowerCase()

    // Filter tools by seq_type compatibility
    def seq_type = (params.seq_type ?: 'dna').toLowerCase()
    if (seq_type == 'rna') {
        def bam_only = ['polysolver', 'kourami']
        def removed = tools_list.findAll { it in bam_only }
        if (removed) log.warn "RNA-seq mode: skipping BAM-only tools: ${removed.join(', ')}"
        tools_list = tools_list.findAll { !(it in bam_only) }
    }
    if (seq_type in ['longreads_hifi', 'longreads_ont']) {
        def kept = ['t1k']
        def removed = tools_list.findAll { !(it in kept) }
        if (removed) log.warn "Long-read mode: skipping short-read tools: ${removed.join(', ')}"
        tools_list = tools_list.findAll { it in kept }
    }

    // -------------------------------------------------------------------------
    // Build input channels
    // Samplesheet takes priority over --input directory
    // -------------------------------------------------------------------------
    def eff_input_type = params.input_type ?: 'bam'
    def ch_input_bam   = Channel.empty()
    def ch_input_fastq = Channel.empty()

    if (params.input_samplesheet) {
        def ss_header = file(params.input_samplesheet).text.readLines()[0]?.toLowerCase() ?: ''
        if (ss_header.contains('bam_path') || ss_header.contains('cram_path')) {
            eff_input_type = ss_header.contains('cram_path') ? 'cram' : 'bam'
            ch_input_bam = Channel.fromPath(params.input_samplesheet)
                .splitCsv(header: true)
                .map { row ->
                    def f = file(row.bam_path ?: row.cram_path)
                    tuple(row.sample_id, f)
                }
        } else {
            eff_input_type = 'fastq'
            ch_input_fastq = Channel.fromPath(params.input_samplesheet)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample_id, file(row.fastq_1), file(row.fastq_2)) }
        }
    } else {
        if (params.input_type == 'fastq') {
            eff_input_type = 'fastq'
            ch_input_fastq = Channel
                .fromFilePairs("${params.input}/*_{R1,R2,1,2}*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: true)
                .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
        } else {
            eff_input_type = params.input_type
            ch_input_bam = Channel
                .fromPath("${params.input}/*.{bam,cram}", checkIfExists: true)
                .map { f -> tuple(f.baseName.replaceAll(/\.(bam|cram)$/, ''), f) }
        }
    }

    // CRAM → BAM conversion when needed
    if (eff_input_type == 'cram') {
        ch_ref = Channel.fromPath(params.reference_fasta, checkIfExists: true).first()
        CRAM_TO_BAM(ch_input_bam, ch_ref)
        ch_input_bam = CRAM_TO_BAM.out.bam
        eff_input_type = 'bam'
    }

    // Resolve modality once from effective input type + seq_type (override with --run_modality)
    def resolved_modality = params.run_modality
        ? params.run_modality.toString().toLowerCase()
        : ((seq_type == 'rna') ? 'rnaseq' : (eff_input_type == 'fastq' ? 'wes' : 'wgs'))
    log.info "Resolved modality   : ${resolved_modality}"

    // Initialize result channels
    ch_optitype_results  = Channel.empty()
    ch_arcashla_results  = Channel.empty()
    ch_spechla_results   = Channel.empty()
    ch_hlahd_results     = Channel.empty()
    ch_polysolver_results = Channel.empty()
    ch_kourami_results   = Channel.empty()
    ch_t1k_results       = Channel.empty()
    ch_seq2hla_results   = Channel.empty()

    if (eff_input_type == 'fastq') {
        // FASTQ path
        if ('optitype' in tools_list) {
            OPTITYPE(ch_input_fastq, params.optitype_seq_type)
            ch_optitype_results = OPTITYPE.out.results
        }
        if ('arcashla' in tools_list) {
            ARCASHLA_FASTQ(ch_input_fastq)
            ch_arcashla_results = ARCASHLA_FASTQ.out.results
        }
        if ('spechla' in tools_list) {
            SPECHLA_FASTQ(ch_input_fastq)
            ch_spechla_results = SPECHLA_FASTQ.out.results
        }
        if ('hlahd' in tools_list) {
            HLAHD_FASTQ(ch_input_fastq, params.hlahd_genes)
            ch_hlahd_results = HLAHD_FASTQ.out.results
        }
        if ('t1k' in tools_list) {
            T1K_FASTQ(ch_input_fastq)
            ch_t1k_results = T1K_FASTQ.out.results
        }
        if ('seq2hla' in tools_list) {
            if (seq_type == 'rna') {
                SEQ2HLA(ch_input_fastq)
                ch_seq2hla_results = SEQ2HLA.out.results
            } else {
                log.warn "seq2HLA is RNA-seq focused — set --seq_type rna to enable it."
            }
        }
        if ('polysolver' in tools_list || 'kourami' in tools_list) {
            log.warn "POLYSOLVER/Kourami require BAM input and are skipped for FASTQ."
        }

    } else {
        // BAM path
        if ('arcashla' in tools_list) {
            ARCASHLA_BAM(ch_input_bam, params.reference_build)
            ch_arcashla_results = ARCASHLA_BAM.out.results
        }
        if ('spechla' in tools_list) {
            SPECHLA_BAM(ch_input_bam, params.reference_build)
            ch_spechla_results = SPECHLA_BAM.out.results
        }
        if ('hlahd' in tools_list && resolved_modality != 'wes') {
            HLAHD_BAM(ch_input_bam, params.reference_build, params.hlahd_genes)
            ch_hlahd_results = HLAHD_BAM.out.results
        }
        if ('polysolver' in tools_list) {
            POLYSOLVER(ch_input_bam)
            ch_polysolver_results = POLYSOLVER.out.results
        }
        if ('kourami' in tools_list) {
            KOURAMI(ch_input_bam)
            ch_kourami_results = KOURAMI.out.results
        }

        // Tools that need FASTQ converted from BAM
        def need_fastq = ('optitype' in tools_list) || ('t1k' in tools_list) || (resolved_modality == 'wes' && 'hlahd' in tools_list)
        def use_region_fastq = params.extract_hla_region || (resolved_modality == 'wes')
        if (need_fastq) {
            if (use_region_fastq) {
                if (resolved_modality == 'wes' && !params.extract_hla_region) {
                    log.info "WES BAM input detected: enabling HLA-region extraction before FASTQ conversion to reduce storage footprint"
                }
                EXTRACT_HLA_AND_CONVERT(ch_input_bam)
                ch_fastq = EXTRACT_HLA_AND_CONVERT.out.reads
            } else {
                BAM_TO_FASTQ(ch_input_bam)
                ch_fastq = BAM_TO_FASTQ.out.reads
            }
            if ('optitype' in tools_list) {
                OPTITYPE(ch_fastq, params.optitype_seq_type)
                ch_optitype_results = OPTITYPE.out.results
            }
            if ('t1k' in tools_list) {
                T1K_FASTQ(ch_fastq)
                ch_t1k_results = T1K_FASTQ.out.results
            }
            if ('hlahd' in tools_list && resolved_modality == 'wes') {
                HLAHD_FASTQ(ch_fastq, params.hlahd_genes)
                ch_hlahd_results = HLAHD_FASTQ.out.results
            }
        }
        if ('seq2hla' in tools_list) {
            log.warn "seq2HLA is FASTQ/RNA-seq oriented — skipped for BAM input."
        }
    }

    // Consensus calling
    if (params.enable_majority_voting) {
        ch_modality = Channel.value(resolved_modality)
        MAJORITY_VOTING_WORKFLOW(
            ch_optitype_results,
            ch_arcashla_results,
            ch_spechla_results,
            ch_hlahd_results,
            ch_polysolver_results,
            ch_kourami_results,
            ch_t1k_results,
            ch_seq2hla_results,
            ch_modality
        )
    }
}

/*
========================================================================================
    COMPLETION MESSAGE
========================================================================================
*/

workflow.onComplete {
    log.info """
    ============================================================
    Pipeline completed!
    ============================================================
    Status      : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Results dir : ${params.outdir}
    ============================================================
    """.stripIndent()
    
    // Performance summary
    if (workflow.success) {
        def tools_used = params.tools.tokenize(',')*.trim()
        def performance_notes = []
        
        if ('spechla' in tools_used && params.input_type != 'fastq') {
            performance_notes << "- SpecHLA: Used BAM input with chromosome 6 extraction (OPTIMIZED)"
        }
        
        if ('arcashla' in tools_used && params.input_type != 'fastq') {
            performance_notes << "- ArcasHLA: Used BAM input directly (OPTIMIZED)"
        }
        
        if ('optitype' in tools_used && params.input_type != 'fastq') {
            if (params.extract_hla_region) {
                performance_notes << "- OptiType: Used HLA region extraction (SPACE EFFICIENT)"
            } else {
                performance_notes << "- OptiType: Full BAM-to-FASTQ conversion performed"
            }
        }
        
        if (performance_notes) {
            log.info """
    PERFORMANCE SUMMARY:
${performance_notes.join('\n    ')}
    ============================================================
            """.stripIndent()
        }
    }
}

workflow.onError {
    log.error """
    ============================================================
    Pipeline execution failed!
    ============================================================
    Error message: ${workflow.errorMessage}
    Error report : ${workflow.errorReport}
    ============================================================
    """.stripIndent()
}
