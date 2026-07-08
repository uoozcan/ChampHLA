# ChampHLA — Champion-Challenger HLA Typing Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
![Version](https://img.shields.io/badge/version-2.0.0-blue)

A reproducible Nextflow DSL2 pipeline integrating eight HLA typing tools with benchmark-calibrated Champion-Challenger ensemble consensus, plus two orthogonal tools (Locityper, Immuannot) for silver-standard truth generation. Supports WGS, WES, RNA-seq, and long-read sequencing.

## Get Running in 5 Minutes

```bash
# 1. Install Nextflow
curl -s https://get.nextflow.io | bash && mv nextflow ~/bin/

# 2. Clone
git clone https://github.com/uoozcan/ChampHLA.git && cd ChampHLA

# 3. Test
nextflow run main.nf -profile test,docker

# 4. Run on real data
# WES:
nextflow run main.nf --input bams/ --input_type bam --tools optitype,hlahd,t1k \
    --hlahd_db /path/to/db --t1k_hlaidx /path/to/idx --enable_majority_voting -profile docker
# RNA-seq:
nextflow run main.nf --input bams/ --input_type bam --seq_type rna \
    --tools optitype,arcashla,hlahd,t1k --optitype_seq_type rna --enable_majority_voting -profile docker
# HPC (SLURM):
nextflow run main.nf -params-file params.yaml -profile slurm,singularity -resume
```

See [docs/INSTALLATION.md](docs/INSTALLATION.md) for cross-platform installation (HPC, Linux, macOS, Windows) and [docs/USAGE.md](docs/USAGE.md) for detailed usage examples.

## Overview

This pipeline integrates multiple HLA typing tools to provide robust and accurate HLA genotyping from sequencing data:

- **OptiType**: High-precision HLA typing for DNA/RNA-seq data
- **ArcasHLA**: Fast HLA typing optimized for RNA-seq but works with DNA data
- **SpecHLA**: Exome/genome HLA typing with variant calling support

In addition to the eight allele callers, two **orthogonal tools** support silver-standard truth generation on cohorts that lack gold HLA truth (see [Silver-Standard Truth Generation](#-silver-standard-truth-generation-locityper--immuannot)):

- **Locityper**: depth-aware locus-pangenome genotyper (WGS + long-read); an independent, alignment-and-depth-based track with a genotype-quality/novelty confidence signal
- **Immuannot**: full-resolution HLA/KIR annotation of **assembled contigs** against IPD-IMGT/HLA gene features (validates gene structure and novelty)

### Key Features

✨ **Multi-tool Integration**: Combines results from multiple HLA typing tools  
**Native Extended Toolset**: Includes HLA-HD, POLYSOLVER, Kourami, T1K, and Seq2HLA execution modules  
🎯 **Consensus Calling**: Optional majority voting for improved accuracy  
📊 **Flexible Input**: Supports BAM, CRAM, and FASTQ files  
🚀 **HPC Optimized**: Pre-configured for SLURM schedulers (CSC Puhti & Roihu)  
🐳 **Container Ready**: Docker and Singularity support  
📈 **Exome Optimized**: Special handling for exome sequencing data  
💾 **Space Efficient**: Intelligent HLA region extraction to minimize storage  

## 📋 Requirements

- **Nextflow**: ≥23.04.0
- **Container Engine**: Docker or Singularity
- **Compute Resources**: 
  - Minimum: 4 CPUs, 16 GB RAM
  - Recommended: 8+ CPUs, 32+ GB RAM

### Pre-downloaded Containers (CSC Puhti / Roihu)

On CSC Puhti and Roihu, the Singularity/Apptainer containers are already available at:
```
/scratch/project_2008084/hla_references/singularity_cache/containers/
├── optitype.sif
├── arcashla.sif
└── spechla.sif
```

## 🚀 Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/uoozcan/ChampHLA.git
cd ChampHLA
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

### 4. On CSC HPC (SLURM: Puhti / Roihu)

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

Specialized / opt-in (see [Silver-Standard Truth Generation](#-silver-standard-truth-generation-locityper--immuannot)): `locityper` (WGS/long-read BAM only), `immuannot` (assembled contigs only). Both are primarily used as orthogonal truth sources rather than as ensemble members.

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

**Locityper** (depth-aware genotyper; **BAM/CRAM, WGS or long-read only** — gated off for WES/RNA)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--locityper_db` | Path to the locus-pangenome DB built by `bin/build_locityper_db.sh` (**required**) | `null` |
| `--reference_fasta` | Reference the BAM was aligned to (**required**) | `null` |
| `--locityper_loci` | Optional comma-separated locus allow-list (e.g. `A,B,C`) | `null` |
| `--seq_type` | Read-technology preset: `dna` (Illumina), `longreads_hifi`, or `longreads_ont` | `dna` |

**Immuannot** (annotates **assembled contigs**, not reads; opt-in)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--contigs_dir` | Directory of assembled contigs (`*.{fa,fasta,fa.gz,fasta.gz}`) to annotate (**required**) | `null` |
| `--immuannot_refdata` | Immuannot reference bundle directory (IPD gene features) (**required**) | `null` |
| `--immuannot_dir` | Optional Immuannot install dir (for `Immuannot.sh`) | `null` |

#### Majority Voting
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--enable_majority_voting` | Enable consensus calling | `false` |
| `--mv_min_tools` | Minimum tools for consensus | `2` |
| `--mv_resolution` | Resolution level (2 or 4) | `2` |
| `--mv_genes` | Genes for voting | `A,B,C,DQA1,DQB1,DRB1` |

Consensus defaults to **HLA-A/B/C** for decision outputs. Detailed per-tool output files still include any additional loci each caller reports.

#### Resource Limits
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--max_cpus` | Maximum CPUs per job | `40` |
| `--max_memory` | Maximum memory per job | `180.GB` |
| `--max_time` | Maximum time per job | `24.h` |

#### HPC Options (CSC Puhti / Roihu)
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
| `roihu` | CSC Roihu HPC environment |
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

## 🥈 Silver-Standard Truth Generation (Locityper + Immuannot)

Many WGS / long-read cohorts have **no gold HLA truth**, which blocks benchmarking the
ensemble. ChampHLA can build a **silver-standard** truth table from two methodologically
**orthogonal** tools — **Locityper** (read/depth-based genotyping) and **Immuannot**
(assembly annotation) — whose *agreement* is a defensible truth signal. This is a *silver*
standard: always validate it against gold truth before trusting it, and never benchmark a
truth-source tool against its own truth (an automatic guardrail enforces this).

**Prerequisite — build the Locityper DB** (once, pinned to an IPD-IMGT/HLA release):

```bash
bash bin/build_locityper_db.sh \
    --out /path/to/locityper_db \
    --ref /path/to/GRCh38.fa \
    --imgt-dir /path/to/IMGTHLA \      # clone of ANHIG/IMGTHLA (pinned tag)
    --imgt-version 3.59.0 \
    --loci A,B,C,DRB1,DQA1,DQB1,DPB1
```

**Example 1 — Locityper standalone** (WGS BAM/CRAM; long-read via `--seq_type longreads_hifi|longreads_ont`):

```bash
nextflow run main.nf \
    --input bam_files/ --input_type bam \
    --tools locityper \
    --locityper_db /path/to/locityper_db \
    --reference_fasta /path/to/GRCh38.fa \
    --seq_type dna \
    --outdir results/ -profile singularity
```

**Example 2 — Immuannot on assembled contigs** (annotation / novelty check):

```bash
nextflow run main.nf \
    --input bam_files/ --input_type bam \
    --tools immuannot \
    --contigs_dir /path/to/assembled_contigs \
    --immuannot_refdata /path/to/immuannot_refdata \
    --outdir results/ -profile singularity
```

**Example 3 — Generate silver truth** (run both truth sources + `--generate_truth`):

```bash
nextflow run main.nf \
    --input bam_files/ --input_type bam \
    --tools locityper,immuannot \
    --locityper_db /path/to/locityper_db \
    --reference_fasta /path/to/GRCh38.fa \
    --contigs_dir /path/to/assembled_contigs \
    --immuannot_refdata /path/to/immuannot_refdata \
    --generate_truth \
    --outdir results/ -profile singularity
```

This writes `results/silver_truth/truth_long.tsv` (+ `truth_provenance.tsv`), keeping only
calls where the two sources agree and pass the gating thresholds:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--generate_truth` | Build silver truth from Locityper + Immuannot | `false` |
| `--truth_min_gq` | Minimum Locityper genotype quality (GQ) to qualify | `20.0` |
| `--truth_max_novelty` | Maximum Locityper novelty score to qualify | `0.05` |
| `--truth_require_agreement` | Require Locityper + Immuannot to agree (`false` = allow a single high-confidence source) | `true` |

**Then**: validate against gold with `bin/validate_silver_truth.py`, and point a benchmark
config's `truth.path` at the generated `truth_long.tsv` (template:
`conf/benchmark_silver_truth_wgs.example.yaml`). The benchmark harness auto-drops any
Locityper/Immuannot runs from evaluation to prevent circularity.

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

## 💡 Tips for CSC HPC (Puhti & Roihu)

### 1. Set Up Your Environment

Load the workflow modules for your system, then run with the matching site profile.

```bash
# --- CSC Roihu ---
module load nextflow/25.10.4.11173 openjdk/25.0.1_8   # apptainer is available by default
module load python-data/3.12                          # for downstream figure/analysis scripts
# run with:  -profile roihu,singularity -params-file conf/roihu_params.yaml

# --- CSC Puhti ---
module load java/21
module load biopython-env/3.10.6
module load nextflow/25.10.0
# run with:  -profile puhti,singularity -params-file conf/puhti_params.yaml

# Set Singularity/Apptainer cache (same path on both systems)
export NXF_SINGULARITY_CACHEDIR=/scratch/project_2008084/hla_references/singularity_cache/containers
```

> Module versions change over time — check `module avail nextflow openjdk` (Roihu) or
> `module avail nextflow` (Puhti) for the currently installed versions.

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
*(indicative; Roihu compute nodes are larger — 384 CPU / 1.5 TB)*

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

## 🔬 Tool Versions

The following tool versions were used in the ChampHLA benchmark (n=131 WGS / 129 WES / 107 RNA-seq samples). Container tags are pinned where available; others were pulled at benchmark time.

| Tool | Version | Container / Source |
|------|---------|-------------------|
| OptiType | 1.3.5 | `fred2/optitype:latest` (Docker Hub) |
| ArcasHLA | 0.5.0 | `nmdpbioinformatics/arcas-hla:latest` (Docker Hub) |
| SpecHLA | v1.1 | No Docker image — install from [deepomicslab/SpecHLA](https://github.com/deepomicslab/SpecHLA) |
| HLA-HD | 1.4.0 | `quay.io/biocontainers/hla-hd:1.4.0--hdfd78af_0` |
| POLYSOLVER | v4 | `sachet/polysolver:v4` (Docker Hub) |
| Kourami | 0.9.6 | `zlskidmore/kourami:latest` (Docker Hub) |
| T1K | 1.0.9 | `quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0` |
| Seq2HLA | 2.3 | `fred2/seq2hla:2.3` (Docker Hub) |
| Locityper | 1.4.5.1 | `eichlerlab/locityper:1.4.5.1` |
| Immuannot | latest | `nmdpbioinformatics/immuannot` |

**SpecHLA installation note:** SpecHLA has no pre-built Docker image. Install natively following the instructions at https://github.com/deepomicslab/SpecHLA, then set `--spechla_path` to the installation directory. On CSC Puhti/Roihu a local installation is available at the path specified in `conf/puhti_params.yaml` / `conf/roihu_params.yaml` (pass `--use_local_spechla true --spechla_path <dir>`).

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use this pipeline in your research, please cite the software and the
individual tools it integrates.

### Pipeline

The ChampHLA manuscript is currently in preparation. A peer-reviewed reference
will be added here upon publication. In the meantime, please cite the software:

```bibtex
@software{champhla,
  author  = {Özcan, Umut Onur},
  title   = {ChampHLA: Champion-Challenger HLA Typing Pipeline},
  year    = {2025},
  url      = {https://github.com/uoozcan/ChampHLA},
  version = {2.0.0}
}
```

```bibtex
@article{champhla_manuscript,
  author  = {Özcan, Umut Onur and others},
  title   = {ChampHLA: benchmark-calibrated Champion-Challenger ensemble consensus for HLA typing},
  journal = {},
  year    = {},
  volume  = {},
  pages   = {},
  doi     = {},
  note    = {Manuscript in preparation}
}
```

### Integrated tools

Please cite the tools you enable via `--tools`, and the workflow engine:

- **OptiType** — Szolek A, Schubert B, Mohr C, Sturm M, Feldhahn M, Kohlbacher O. OptiType: precision HLA typing from next-generation sequencing data. *Bioinformatics.* 2014;30(23):3310–3316.
- **arcasHLA** — Orenbuch R, Filip I, Comito D, Shaman J, Pe'er I, Rabadan R. arcasHLA: high-resolution HLA typing from RNAseq. *Bioinformatics.* 2020;36(1):33–40.
- **SpecHLA** — Wang S, Wang M, Chen L, Pan G, Wang Y, Li SC. SpecHLA enables full-resolution HLA typing from sequencing data. *Cell Reports Methods.* 2023;3(9):100589.
- **HLA-HD** — Kawaguchi S, Higasa K, Shimizu M, Yamada R, Matsuda F. HLA-HD: An accurate HLA typing algorithm for next-generation sequencing data. *Human Mutation.* 2017;38(7):788–797.
- **POLYSOLVER** — Shukla SA, Rooney MS, Rajasagi M, et al. Comprehensive analysis of cancer-associated somatic mutations in class I HLA genes. *Nature Biotechnology.* 2015;33(11):1152–1158.
- **Kourami** — Lee H, Kingsford C. Kourami: graph-guided assembly for novel human leukocyte antigen allele discovery. *Genome Biology.* 2018;19(1):16.
- **T1K** — Song L, Bai G, Liu XS, Li B, Li H. Efficient and accurate KIR and HLA genotyping with massively parallel sequencing data. *Genome Research.* 2023;33(6):923–931.
- **seq2HLA** — Boegel S, Löwer M, Schäfer M, et al. HLA typing from RNA-Seq sequence reads. *Genome Medicine.* 2012;4(12):102.
- **Locityper** — Prodanov T, Plender EG, Seebohm G, Meuth SG, Eichler EE, Marschall T. Locityper enables targeted genotyping of complex polymorphic genes. *Nature Genetics.* 2025. doi:10.1038/s41588-025-02362-4.
- **Immuannot** — Zhou Y, Song L, Li H. Full-resolution HLA and KIR gene annotations for human genome assemblies. *Genome Research.* 2024;34(11):1931–1941.
- **Nextflow** (workflow engine) — Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. *Nature Biotechnology.* 2017;35(4):316–319.

## 🙏 Acknowledgments

- CSC - IT Center for Science, Finland for providing HPC resources
- The developers of the integrated tools: OptiType, arcasHLA, SpecHLA, HLA-HD, POLYSOLVER, Kourami, T1K, seq2HLA, Locityper, and Immuannot
- The Nextflow team for the workflow framework

## 📮 Contact

For questions or issues, please open an issue on GitHub or contact:
- **GitHub**: [@uoozcan](https://github.com/uoozcan)

---

**Version**: 2.0.0  
**Last Updated**: July 2026
