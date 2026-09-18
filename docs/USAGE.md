# ChampHLA Usage Guide

## Input Formats

ChampHLA accepts three input methods:

### Directory scan (`--input` + `--input_type`)

| Input type | Flag | Expected files |
|---|---|---|
| BAM | `--input_type bam` | `*.bam` with `.bai` index in the same directory |
| CRAM | `--input_type cram` | `*.cram` (requires `--reference_fasta`) |
| FASTQ | `--input_type fastq` | Paired-end: `*_{R1,R2,1,2}*.{fastq,fq,fastq.gz,fq.gz}` |

### Samplesheet (`--input_samplesheet`)

A CSV file mapping sample IDs to file paths. Header determines format:

**BAM samplesheet:**
```csv
sample_id,bam_path
NA12878,/data/NA12878.bam
NA12891,/data/NA12891.bam
```

**FASTQ samplesheet:**
```csv
sample_id,fastq_1,fastq_2
NA12878,/data/NA12878_R1.fastq.gz,/data/NA12878_R2.fastq.gz
```

---

## Basic Usage by Modality

### Whole-Exome Sequencing (WES)

```bash
nextflow run main.nf \
    --input bam_files/ \
    --input_type bam \
    --seq_type dna \
    --tools optitype,hlahd,t1k,spechla,polysolver,kourami,arcashla \
    --spechla_exon_only 1 \
    --extract_hla_region \
    --enable_majority_voting \
    --outdir results/ \
    -profile docker
```

Key WES flags:
- `--spechla_exon_only 1` — required for exome data when using SpecHLA
- `--extract_hla_region` — extracts chr6:28–34 Mb before FASTQ conversion (saves disk)

### Whole-Genome Sequencing (WGS)

```bash
nextflow run main.nf \
    --input bam_files/ \
    --input_type bam \
    --seq_type dna \
    --tools optitype,hlahd,t1k,spechla,kourami,arcashla \
    --enable_majority_voting \
    --outdir results/ \
    -profile docker
```

### RNA-seq

```bash
nextflow run main.nf \
    --input bam_files/ \
    --input_type bam \
    --seq_type rna \
    --tools optitype,arcashla,hlahd,t1k,spechla,seq2hla \
    --optitype_seq_type rna \
    --enable_majority_voting \
    --outdir results/ \
    -profile docker
```

Notes:
- `--seq_type rna` automatically skips POLYSOLVER and Kourami (DNA-only tools)
- `--optitype_seq_type rna` activates OptiType's RNA mode
- ArcasHLA performs best on RNA-seq; its WGS/WES accuracy is near-zero

### Long-Read Sequencing

```bash
nextflow run main.nf \
    --input fastq_files/ \
    --input_type fastq \
    --seq_type longreads_hifi \
    --tools t1k \
    --t1k_hlaidx /path/to/t1k_hla.fa \
    --outdir results/ \
    -profile docker
```

T1K is currently the only tool with validated long-read support. Use `longreads_hifi` for PacBio HiFi or `longreads_ont` for Oxford Nanopore.

---

## Consensus Voting

Enable ensemble consensus calling with `--enable_majority_voting`. ChampHLA supports two weighting strategies:

### Calibrated weights (default, recommended)

```bash
--enable_majority_voting \
--weighting calibrated
```

Pre-computed weights from the 1000 Genomes benchmark are applied automatically based on the detected modality. Weight files are bundled in `conf/tool_weights_{wgs,wes,rna}.json`.

### Equal weights

```bash
--enable_majority_voting \
--weighting equal
```

All tools contribute equally to the vote.

### Controlling consensus scope

```bash
--mv_genes A,B,C          # genes to include in consensus (default)
--mv_min_tools 2           # minimum tools required for a consensus call
--mv_resolution 2          # allele resolution: 2-field (default) or 4-field
```

---

## HPC Submission (SLURM)

### Single-sample SLURM job

```bash
#!/bin/bash
#SBATCH --job-name=champhla
#SBATCH --account=project_XXXXXXX
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out

module load java/17
module load nextflow

export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_cache

nextflow run main.nf \
    --input /path/to/sample.bam \
    --input_type bam \
    --tools optitype,hlahd,t1k \
    --hlahd_db /path/to/hlahd_db \
    --t1k_hlaidx /path/to/t1k_hla.fa \
    --outdir results/ \
    --max_cpus 16 \
    --max_memory 64.GB \
    -profile slurm,singularity \
    -resume
```

### Array job for multiple samples

```bash
#!/bin/bash
#SBATCH --job-name=champhla_array
#SBATCH --array=1-30
#SBATCH --account=project_XXXXXXX
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_list.txt)

nextflow run main.nf \
    --input /data/${SAMPLE}.bam \
    --input_type bam \
    --tools optitype,hlahd,t1k \
    --outdir results/${SAMPLE} \
    -profile slurm,singularity \
    -resume
```

### Resource tuning

Nextflow automatically retries failed jobs with increased resources. Override defaults:

```bash
--max_cpus 40          # cap per-job CPUs (default: 40)
--max_memory 180.GB    # cap per-job memory (default: 180 GB)
--max_time 24.h        # cap per-job time (default: 24 h)
```

---

## Output Structure

```
results/
├── sample1/
│   ├── optitype/
│   │   ├── sample1_result.tsv
│   │   └── sample1_coverage_plot.pdf
│   ├── arcashla/
│   │   └── sample1.genotype.json
│   ├── hlahd/
│   │   └── sample1_final.result.txt
│   ├── spechla/
│   │   └── sample1_spechla.txt
│   ├── t1k/
│   │   └── sample1_t1k.tsv
│   └── consensus/
│       ├── sample1.consensus.tsv
│       └── sample1.summary.txt
├── majority_voting/
│   ├── all_samples.consensus.tsv
│   └── consensus_summary.txt
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.html
```

Key output files:
- **Per-tool results** — each tool writes its native output plus a normalised TSV
- **`consensus/`** — per-sample consensus calls when `--enable_majority_voting` is set
- **`majority_voting/`** — cohort-wide aggregated consensus
- **`pipeline_info/`** — Nextflow execution metrics (timing, resource usage, DAG)

---

## Troubleshooting

### "Input parameter is required"

Ensure parameters are correctly passed. Multi-line commands need `\` at the end of each line:
```bash
# Correct
nextflow run main.nf \
    --input samples/ \
    --input_type fastq

# Wrong — each line runs as a separate command
nextflow run main.nf
    --input samples/
```

### Container not found (Singularity)

Set the Singularity cache directory and ensure SIF files exist:
```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/cache
ls $NXF_SINGULARITY_CACHEDIR/*.sif
```

### SpecHLA produces no results on exome data

Enable exon-only mode:
```bash
--spechla_exon_only 1
```

### Out of memory

Reduce parallelism or exclude memory-heavy tools:
```bash
--max_cpus 8
--max_memory 32.GB
# Or use a lighter tool panel:
--tools optitype,arcashla,hlahd
```

### Resume not working

Nextflow caches results by content hash in the `work/` directory. If the `work/` directory was deleted, caching is lost. Never delete `work/` during an active run.

```bash
# Resume a failed run
nextflow run main.nf <your_params> -resume
```

### HLA-HD database error

HLA-HD requires a bowtie2-indexed database. Provide it via:
```bash
--hlahd_db /path/to/hlahd_database_directory
```

### Viewing logs

```bash
# Nextflow log
cat .nextflow.log

# SLURM job logs
cat logs/champhla_*.out

# Per-tool logs within results
cat results/sample_name/tool_name/*.log
```
