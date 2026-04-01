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
include { EXTRACT_HLA_REGION; BAM_TO_FASTQ; EXTRACT_HLA_AND_CONVERT } from './modules/bam_to_fastq'
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
    
    Required Arguments:
      --input              Path to input files (directory or file)
      --input_type         Input type: 'bam', 'cram', or 'fastq'
    
    Pipeline Options:
      --tools              Comma-separated list of tools to run
                           Available: optitype,arcashla,spechla,hlahd,polysolver,kourami,t1k,seq2hla
                           Default: optitype,arcashla,hlahd
      --outdir             Output directory [default: ./results]
      --reference_build    Reference genome build: hg19 or hg38 [default: hg38]
    
    Tool-specific Options:
      --optitype_seq_type  OptiType sequence type: 'dna' or 'rna' [default: dna]
      --arcashla_genes     Genes for ArcasHLA [default: A,B,C,DQA1,DQB1,DRB1]
      --spechla_genes      Genes for SpecHLA [default: A,B,C,DQA1,DQB1,DRB1]
      --spechla_exon_only  SpecHLA exon-only mode (0 or 1) [default: 0]
                           Use 1 for exome data!
    
    Majority Voting Options:
      --enable_majority_voting  Enable consensus calling [default: false]
      --mv_min_tools       Minimum tools for consensus [default: 2]
      --mv_resolution      Resolution level: 2 or 4 [default: 4]
      --mv_genes           Genes for voting [default: A,B,C,DQA1,DQB1,DRB1]
    
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
      # RNA-seq data with OptiType and ArcasHLA
      nextflow run main.nf \\
        --input samples/ \\
        --input_type fastq \\
        --tools optitype,arcashla \\
        --optitype_seq_type rna \\
        -profile singularity
      
      # Exome data with SpecHLA (exon-only mode)
      nextflow run main.nf \\
        --input bam_files/ \\
        --input_type bam \\
        --tools spechla \\
        --spechla_exon_only 1 \\
        --enable_majority_voting \\
        -profile puhti,singularity
      
      # DNA-seq with all tools and consensus calling
      nextflow run main.nf \\
        --input samples/ \\
        --input_type fastq \\
        --tools optitype,arcashla,spechla \\
        --enable_majority_voting \\
        --outdir results/ \\
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
if (!params.input) {
    log.error "ERROR: --input parameter is required!"
    helpMessage()
    exit 1
}

if (!params.input_type) {
    log.error "ERROR: --input_type parameter is required (bam, cram, or fastq)!"
    helpMessage()
    exit 1
}

// Print parameter summary
log.info """
========================================================================================
    HLA TYPING MULTI-TOOL PIPELINE
========================================================================================
Input               : ${params.input}
Input type          : ${params.input_type}
Output dir          : ${params.outdir}
Tools               : ${params.tools}
Reference build     : ${params.reference_build}
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
    tools_list = params.tools.tokenize(',')*.trim()*.toLowerCase()

    // Create input channel based on input type
    if (params.input_type == 'fastq') {
        // FASTQ files: expect paired-end reads
        ch_input = Channel
            .fromFilePairs("${params.input}/*_{R1,R2,1,2}*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: true)
            .map { sample_id, reads ->
                tuple(sample_id, reads[0], reads[1])
            }
    } else if (params.input_type == 'bam' || params.input_type == 'cram') {
        // BAM/CRAM files
        ch_input = Channel
            .fromPath("${params.input}/*.{bam,cram}", checkIfExists: true)
            .map { file ->
                def sample_id = file.baseName.replaceAll(/\.(bam|cram)$/, '')
                tuple(sample_id, file)
            }
    }

    // Initialize result channels
    ch_optitype_results = Channel.empty()
    ch_arcashla_results = Channel.empty()
    ch_spechla_results = Channel.empty()
    ch_hlahd_results = Channel.empty()
    ch_polysolver_results = Channel.empty()
    ch_kourami_results = Channel.empty()
    ch_t1k_results = Channel.empty()
    ch_seq2hla_results = Channel.empty()

    if (params.input_type == 'fastq') {
        // FASTQ input
        if ('optitype' in tools_list) {
            OPTITYPE(ch_input, params.optitype_seq_type)
            ch_optitype_results = OPTITYPE.out.results
        }

        if ('arcashla' in tools_list) {
            ARCASHLA_FASTQ(ch_input)
            ch_arcashla_results = ARCASHLA_FASTQ.out.results
        }

        if ('spechla' in tools_list) {
            SPECHLA_FASTQ(ch_input)
            ch_spechla_results = SPECHLA_FASTQ.out.results
        }

        if ('hlahd' in tools_list) {
            HLAHD_FASTQ(ch_input, params.hlahd_genes)
            ch_hlahd_results = HLAHD_FASTQ.out.results
        }

        if ('t1k' in tools_list) {
            T1K_FASTQ(ch_input)
            ch_t1k_results = T1K_FASTQ.out.results
        }

        if ('seq2hla' in tools_list) {
            if ((params.seq_type ?: '').toLowerCase() == 'rna') {
                SEQ2HLA(ch_input)
                ch_seq2hla_results = SEQ2HLA.out.results
            } else {
                log.warn "seq2HLA is RNA-seq focused. Set --seq_type rna to run seq2HLA with FASTQ input."
            }
        }

        if ('polysolver' in tools_list || 'kourami' in tools_list) {
            log.warn "POLYSOLVER/Kourami require BAM input and are skipped for FASTQ input."
        }

    } else {
        // BAM/CRAM input
        if ('arcashla' in tools_list) {
            ARCASHLA_BAM(ch_input, params.reference_build)
            ch_arcashla_results = ARCASHLA_BAM.out.results
        }

        if ('spechla' in tools_list) {
            SPECHLA_BAM(ch_input, params.reference_build)
            ch_spechla_results = SPECHLA_BAM.out.results
        }

        if ('hlahd' in tools_list) {
            HLAHD_BAM(ch_input, params.reference_build, params.hlahd_genes)
            ch_hlahd_results = HLAHD_BAM.out.results
        }

        if ('polysolver' in tools_list) {
            POLYSOLVER(ch_input)
            ch_polysolver_results = POLYSOLVER.out.results
        }

        if ('kourami' in tools_list) {
            KOURAMI(ch_input)
            ch_kourami_results = KOURAMI.out.results
        }

        // Tools that require FASTQ conversion from BAM/CRAM
        need_fastq = ('optitype' in tools_list) || ('t1k' in tools_list)
        if (need_fastq) {
            if (params.extract_hla_region) {
                EXTRACT_HLA_AND_CONVERT(ch_input)
                ch_fastq = EXTRACT_HLA_AND_CONVERT.out.reads
            } else {
                BAM_TO_FASTQ(ch_input)
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
        }

        if ('seq2hla' in tools_list) {
            log.warn "seq2HLA is FASTQ/RNA-seq oriented and is skipped for BAM/CRAM input."
        }
    }

    // Majority voting if enabled
    if (params.enable_majority_voting) {
        MAJORITY_VOTING_WORKFLOW(
            ch_optitype_results,
            ch_arcashla_results,
            ch_spechla_results,
            ch_hlahd_results,
            ch_polysolver_results,
            ch_kourami_results,
            ch_t1k_results,
            ch_seq2hla_results
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
