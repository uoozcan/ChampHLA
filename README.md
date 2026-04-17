# HLA Typing Multi-Tool Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
![Version](https://img.shields.io/badge/version-2.0.0-blue)

A comprehensive Nextflow pipeline for HLA typing using multiple state-of-the-art tools, optimized for high-performance computing environments like CSC Puhti.

## 🔬 Overview

This pipeline integrates multiple HLA typing tools to provide robust and accurate HLA genotyping from sequencing data:

- **OptiType**: High-precision HLA typing for DNA/RNA-seq data
- **ArcasHLA**: Fast HLA typing optimized for RNA-seq but works with DNA data
- **SpecHLA**: Exome/genome HLA typing with variant calling support

### Key Features

✨ **Multi-tool Integration**: Combines results from multiple HLA typing tools  
**Native Extended Toolset**: Includes HLA-HD, POLYSOLVER, Kourami, T1K, and Seq2HLA execution modules  
🎯 **Consensus Calling**: Optional majority voting for improved accuracy  
📊 **Flexible Input**: Supports BAM, CRAM, and FASTQ files  
🚀 **HPC Optimized**: Pre-configured for SLURM schedulers (CSC Puhti)  
🐳 **Container Ready**: Docker and Singularity support  
📈 **Exome Optimized**: Special handling for exome sequencing data  
💾 **Space Efficient**: Intelligent HLA region extraction to minimize storage  

## 📋 Requirements

- **Nextflow**: ≥23.04.0
- **Container Engine**: Docker or Singularity
- **Compute Resources**: 
  - Minimum: 4 CPUs, 16 GB RAM
  - Recommended: 8+ CPUs, 32+ GB RAM

### Pre-downloaded Containers (CSC Puhti)

If you're using CSC Puhti, the Singularity containers are already available at:
```
/scratch/project_2008084/hla_references/singularity_cache/containers/
├── optitype.sif
├── arcashla.sif
└── spechla.sif
```

## 🚀 Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/hla-typing-pipeline.git
cd hla-typing-pipeline
```

### 2. Test Installation

```bash
# Check Nextflow version
nextflow -version

# Show help message
nextflow run main.nf --help
```

### 3. Run the Pipeline

#### Example 1: RNA-seq Data (FASTQ)

```bash
nextflow run main.nf \
    --input samples/ \
    --input_type fastq \
    --tools optitype,arcashla \
    --optitype_seq_type rna \
    --outdir results/ \
    -profile docker
```

#### Example 2: Exome Data (BAM) with SpecHLA

```bash
nextflow run main.nf \
    --input bam_files/ \
    --input_type bam \
    --tools spechla \
    --spechla_exon_only 1 \
    --outdir results/ \
    -profile singularity
```

#### Example 3: DNA-seq with All Tools + Consensus

```bash
nextflow run main.nf \
    --input samples/ \
    --input_type fastq \
    --tools optitype,arcashla,spechla,hlahd,t1k \
    --enable_majority_voting \
    --outdir results/ \
    -profile docker -resume
```

### 4. On CSC Puhti (SLURM)

```bash
# Edit submit_slurm.sh with your parameters
nano submit_slurm.sh

# Submit to SLURM
sbatch submit_slurm.sh
```

## 📖 Usage

### Input Requirements

#### FASTQ Files
Paired-end FASTQ files with standard naming:
```
sample1_R1.fastq.gz, sample1_R2.fastq.gz
sample2_R1.fastq.gz, sample2_R2.fastq.gz
```

Acceptable patterns: `*_{R1,R2,1,2}*.{fastq,fq,fastq.gz,fq.gz}`

#### BAM/CRAM Files
```
sample1.bam
sample2.bam
```

Index files (`.bai`) should be in the same directory.

### Parameters

#### Required Parameters
| Parameter | Description |
|-----------|-------------|
| `--input` | Path to input directory or file |
| `--input_type` | Input type: `bam`, `cram`, or `fastq` |

#### Tool Selection
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--tools` | Comma-separated list of tools | `optitype,arcashla,hlahd` |

Available tools: `optitype`, `arcashla`, `spechla`, `hlahd`, `polysolver`, `kourami`, `t1k`, `seq2hla`

#### Tool-Specific Options

**OptiType**
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--optitype_seq_type` | Sequence type: `dna` or `rna` | `dna` |

**ArcasHLA**
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--arcashla_genes` | HLA genes to type | `A,B,C,DQA1,DQB1,DRB1` |

**SpecHLA**
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--spechla_genes` | HLA genes to type | `A,B,C,DQA1,DQB1,DRB1` |
| `--spechla_exon_only` | Exon-only mode (0 or 1) | `0` |

⚠️ **Important**: Set `--spechla_exon_only 1` for exome data!

#### Majority Voting
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--enable_majority_voting` | Enable consensus calling | `false` |
| `--mv_min_tools` | Minimum tools for consensus | `2` |
| `--mv_resolution` | Resolution level (2 or 4) | `2` |
| `--mv_genes` | Genes for voting (consensus focus) | `A,B,C` |

Consensus defaults to **HLA-A/B/C** for decision outputs. Detailed per-tool output files still include any additional loci each caller reports.

#### Resource Limits
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--max_cpus` | Maximum CPUs per job | `40` |
| `--max_memory` | Maximum memory per job | `180.GB` |
| `--max_time` | Maximum time per job | `24.h` |

#### HPC Options (CSC Puhti)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--slurm_account` | SLURM account | Required on Puhti |
| `--slurm_partition` | SLURM partition | `small` |

### Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Use Docker containers |
| `singularity` | Use Singularity containers |
| `puhti` | CSC Puhti HPC environment |
| `slurm` | Generic SLURM cluster |
| `test` | Small test dataset |

Profiles can be combined: `-profile puhti,singularity`

## 📁 Output Structure

```
results/
├── sample1/
│   ├── optitype/
│   │   ├── sample1_result.tsv
│   │   ├── sample1_coverage_plot.pdf
│   │   └── sample1.optitype.log
│   ├── arcashla/
│   │   ├── sample1.genotype.json
│   │   └── sample1.arcashla.log
│   ├── spechla/
│   │   ├── sample1_spechla.tsv
│   │   └── sample1_spechla.log
│   └── consensus/
│       ├── sample1.consensus.tsv
│       └── sample1.summary.txt
├── sample2/
│   └── ...
├── majority_voting/
│   ├── all_samples.consensus.tsv
│   └── consensus_summary.txt
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    └── execution_trace.txt
```

## 📊 1000 Genomes Benchmark Workflow

Scientific benchmarking now targets only real 1000 Genomes samples with public HLA ground truth and matched WGS, WES, and RNA-seq availability. The synthetic fixture in `tests/fixtures/` is retained only for CI and parser validation.

Build canonical manifests:

```bash
python3 bin/build_1000g_benchmark_manifests.py     --truth /path/to/1000g_hla_truth.tsv     --sequencing /path/to/1000g_sequencing_source.tsv     --output-dir results/1000g_manifests     --acquisition-date 2014-07-25     --supported-loci A,B,C,DRB1,DQB1
```

Run the split-aware real-data benchmark:

```bash
python3 bin/run_1000g_benchmark.py     --config conf/benchmark_1000g_config.example.yaml     --output-dir results/1000g_benchmark
```

Phase-gated real-data input generation from finished outputs:

```bash
python3 bin/build_1000g_phase_gated_inputs.py \
    --truth-csv /scratch/project_2008084/ozcanumu/hla_calibration/conf/ground_truth_data.csv \
    --wgs-results /scratch/project_2008084/hla_calibration/wgs/results \
    --wes-results /scratch/project_2008084/hla_calibration/wes_3sample/results \
    --rnaseq-results /scratch/project_2008084/hla_calibration/rna_3sample/results \
    --samples NA06985,NA06986,NA06994 \
    --supported-loci A,B,C,DRB1,DQB1 \
    --output-dir /scratch/project_2008084/pihla-publish/analysis/1000g_realdata/phase_gated_inputs
```

Then point `manifests.sequencing_source` to the generated `sequencing_source.tsv` and `truth.path` to `truth_long.tsv`.

The benchmark runner:
- assembles `truth_manifest.tsv`, `sequencing_manifest.tsv`, and `cohort_manifest.tsv`
- restricts analysis to strict tri-modal truth-backed samples
- learns confidence weights on training samples
- tunes abstention support on validation when available
- reports single-tool, majority-vote, and weighted-consensus results on holdout only

See `docs/BENCHMARK_FIGURES.md` for the real-data benchmark interface and outputs.

## 🔧 Advanced Usage

### Resume a Previous Run

```bash
nextflow run main.nf \
    --input samples/ \
    --input_type fastq \
    -resume
```

### Extract Only HLA Region (Space Efficient)

```bash
nextflow run main.nf \
    --input bam_files/ \
    --input_type bam \
    --tools optitype \
    --extract_hla_region \
    -profile singularity
```

### Custom HLA Region

```bash
nextflow run main.nf \
    --input bam_files/ \
    --input_type bam \
    --extract_hla_region \
    --hla_region_start 29000000 \
    --hla_region_end 33000000 \
    -profile singularity
```

### Save Intermediate Files

```bash
nextflow run main.nf \
    --input samples/ \
    --input_type bam \
    --save_intermediate \
    -profile docker
```

## 💡 Tips for CSC Puhti Users

### 1. Set Up Your Environment

```bash
# Load required modules
module load java/21
module load biopython-env/3.10.6
module load nextflow/25.10.0

# Set Singularity cache
export NXF_SINGULARITY_CACHEDIR=/scratch/project_2008084/hla_references/singularity_cache/containers
```

### 2. Edit SLURM Submission Script

The `submit_slurm.sh` script has all the parameters you need to modify:

```bash
# Open in editor
nano submit_slurm.sh

# Modify these sections:
INPUT_DIR="/scratch/project_2008084/your_samples"
INPUT_TYPE="fastq"
TOOLS="optitype,arcashla"
SLURM_ACCOUNT="project_2008084"  # Your project number
```

### 3. Submit Job

```bash
# Create logs directory
mkdir -p logs

# Submit to SLURM
sbatch submit_slurm.sh

# Check job status
squeue -u $USER

# Monitor log (replace JOBID)
tail -f logs/hla_pipeline_JOBID.out
```

### 4. For Exome Data

**Critical**: Always use `--spechla_exon_only 1` for exome data!

```bash
# In submit_slurm.sh, set:
SPECHLA_EXON_ONLY=1
TOOLS="spechla"  # or include with other tools
```

## 🔍 Troubleshooting

### Common Issues

#### 1. "Input parameter is required"
Make sure parameters are properly passed to Nextflow. Check your SLURM script has backslashes at line ends.

```bash
# ✅ CORRECT
nextflow run main.nf \
    --input samples/ \
    --input_type fastq

# ❌ WRONG (missing backslash causes parameters to run as separate commands)
nextflow run main.nf
    --input samples/
    --input_type fastq
```

#### 2. Container Not Found

```bash
# Set Singularity cache directory
export NXF_SINGULARITY_CACHEDIR=/scratch/project_2008084/hla_references/singularity_cache/containers
```

#### 3. SpecHLA No Results on Exome Data

Enable exon-only mode:
```bash
--spechla_exon_only 1
```

#### 4. Out of Memory

Reduce the number of parallel jobs:
```bash
# In nextflow.config or command line
--max_cpus 20
--max_memory 90.GB
```

#### 5. Resume Not Working

Clean work directory and restart:
```bash
rm -rf work/
nextflow run main.nf <your_params>
```

### Get Help

View the help message:
```bash
nextflow run main.nf --help
```

Check logs:
```bash
# Nextflow log
cat .nextflow.log

# SLURM logs
cat logs/hla_pipeline_*.out
cat logs/hla_pipeline_*.err

# Individual process logs
cat results/sample_name/tool_name/*.log
```

## 📊 Performance

### Typical Run Times (CSC Puhti)

| Data Type | Tools | Samples | Time | Resources |
|-----------|-------|---------|------|-----------|
| RNA-seq FASTQ | OptiType + ArcasHLA | 10 | ~30 min | 8 CPUs, 32 GB |
| Exome BAM | SpecHLA (exon-only) | 10 | ~2 hours | 8 CPUs, 32 GB |
| WGS FASTQ | All tools + voting | 5 | ~4 hours | 16 CPUs, 64 GB |

### Optimization Tips

1. **For Exome Data**: Always use `--spechla_exon_only 1`
2. **For BAM Input**: Use `--extract_hla_region` to save disk space
3. **For Large Cohorts**: Use `-resume` to recover from failures
4. **For Speed**: Request more CPUs with `--max_cpus`

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use this pipeline in your research, please cite:

```

}
```

And cite the individual tools used:

- **OptiType**: Szolek et al. (2014) Bioinformatics
- **ArcasHLA**: Orenbuch et al. (2020) Bioinformatics
- **SpecHLA**: Nariai et al. (2015) BMC Genomics

## 🙏 Acknowledgments

- CSC - IT Center for Science, Finland for providing HPC resources
- The developers of OptiType, ArcasHLA, and SpecHLA

## 📮 Contact

For questions or issues, please open an issue on GitHub or contact:
- **GitHub**: [@ozcanumu](https://github.com/yourusername)

---

**Version**: 2.0.0  
**Last Updated**: April 2026
