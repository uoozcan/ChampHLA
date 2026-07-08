# ChampHLA Installation Guide

This guide covers installing ChampHLA across different environments: HPC clusters, local Linux/macOS, and Windows via WSL2.

## Prerequisites

| Requirement | Version | Notes |
|---|---|---|
| Nextflow | ≥ 23.04.0 | Java-based workflow engine |
| Java | ≥ 11 | Required by Nextflow |
| Docker **or** Singularity | Docker ≥ 20.10 / Singularity ≥ 3.7 | Container runtime |
| Git | any | To clone the repository |

## Quick Install

```bash
# 1. Install Nextflow (if not already installed)
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/   # or any directory on your PATH

# 2. Clone ChampHLA
git clone https://github.com/uoozcan/ChampHLA.git
cd ChampHLA

# 3. Verify with the test profile
nextflow run main.nf -profile test,docker --help
```

## Selecting Which Tools to Install

ChampHLA integrates eight HLA typing tools, but you **do not need to install all of them**. The `--tools` parameter controls which tools run, and you only need containers and databases for the tools you select.

| Tool | Container available | Extra database required | Modalities |
|---|---|---|---|
| OptiType | Yes | No | WGS, WES, RNA |
| ArcasHLA | Yes | No | RNA (primary), WGS/WES (low accuracy) |
| HLA-HD | Yes | Yes — bowtie2 HLA DB | WGS, WES, RNA |
| T1K | Yes | Yes — HLA index file | WGS, WES, RNA, long-reads |
| SpecHLA | Yes (new) / local install | No (bundled) | WGS, WES, RNA |
| Kourami | Yes | Optional — graph DB | WGS, WES |
| POLYSOLVER | Yes | No | WES |
| Seq2HLA | Yes | No | RNA |

A minimal high-accuracy panel for WES:
```bash
--tools optitype,hlahd,t1k
```

A minimal RNA-seq panel:
```bash
--tools optitype,arcashla,hlahd --optitype_seq_type rna
```

---

## HPC Cluster (SLURM + Singularity)

Most HPC clusters do not allow Docker. Use Singularity instead.

### 1. Install Nextflow

```bash
# Many HPC systems provide Nextflow as a module
module load java/17    # or java/11, java/21
module load nextflow   # if available

# Otherwise, install manually
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow ~/bin/
```

### 2. Clone the repository

```bash
git clone https://github.com/uoozcan/ChampHLA.git
cd ChampHLA
```

### 3. Pre-pull Singularity containers

Docker images are automatically converted to Singularity SIF files on first use, but this is slow and may fail on login nodes. Pre-pull to a shared directory:

```bash
# Create a shared cache directory
export NXF_SINGULARITY_CACHEDIR=/path/to/shared/singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

# Pull the containers you need
singularity pull $NXF_SINGULARITY_CACHEDIR/optitype.sif docker://fred2/optitype:latest
singularity pull $NXF_SINGULARITY_CACHEDIR/arcashla.sif docker://nmdpbioinformatics/arcas-hla:latest
singularity pull $NXF_SINGULARITY_CACHEDIR/hlahd.sif docker://quay.io/biocontainers/hla-hd:1.4.0--hdfd78af_0
singularity pull $NXF_SINGULARITY_CACHEDIR/t1k.sif docker://quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0
singularity pull $NXF_SINGULARITY_CACHEDIR/polysolver.sif docker://sachet/polysolver:v4
singularity pull $NXF_SINGULARITY_CACHEDIR/kourami.sif docker://zlskidmore/kourami:latest
singularity pull $NXF_SINGULARITY_CACHEDIR/seq2hla.sif docker://fred2/seq2hla:2.3
```

### 4. Set up tool-specific databases

**HLA-HD** requires a bowtie2-indexed HLA database:
```bash
# Download from the HLA-HD website or build from IMGT/HLA
--hlahd_db /path/to/hlahd_database
```

**T1K** requires an HLA index file:
```bash
# Bundled with T1K container at /usr/local/share/t1k/hlaidx
# Or provide your own:
--t1k_hlaidx /path/to/t1k_hla.fa
```

### 5. Configure and run

Create a params YAML file (recommended over command-line arguments for HPC):

```yaml
# params.yaml
input_samplesheet: '/path/to/samples.csv'
input_type: 'bam'
seq_type: 'dna'
tools: 'optitype,hlahd,t1k,spechla'
outdir: '/path/to/results'
singularity_cache_dir: '/path/to/shared/singularity_cache'
hlahd_db: '/path/to/hlahd_database'
t1k_hlaidx: '/path/to/t1k_hla.fa'
max_cpus: 16
max_memory: '64.GB'
max_time: '24.h'
slurm_account: 'your_project'
slurm_partition: 'small'
```

Submit via SLURM:
```bash
#!/bin/bash
#SBATCH --job-name=champhla
#SBATCH --account=your_project
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/champhla_%j.out

module load java/17
module load nextflow

export NXF_SINGULARITY_CACHEDIR=/path/to/shared/singularity_cache

nextflow run main.nf \
    -params-file params.yaml \
    -profile slurm,singularity \
    -resume
```

The Nextflow head process needs only ~4 CPUs and 16 GB; individual tool jobs are submitted as separate SLURM tasks with resources defined in `nextflow.config`.

### 6. CSC Puhti users

A pre-tuned `puhti` profile is included. See [INSTALL_PUHTI.md](INSTALL_PUHTI.md) for Puhti-specific instructions, including pre-cached containers and module loading.

```bash
module load java/21
module load nextflow/25.10.0
export NXF_SINGULARITY_CACHEDIR=/scratch/project_2008084/hla_references/singularity_cache/containers

nextflow run main.nf \
    -params-file conf/puhti_params.yaml \
    -profile puhti,singularity \
    -resume
```

---

## HPC Cluster (PBS/LSF)

ChampHLA's `slurm` profile can be adapted for PBS or LSF by changing the executor in `nextflow.config`:

```groovy
// For PBS/Torque:
process {
    executor = 'pbs'
    queue = 'batch'
}

// For LSF:
process {
    executor = 'lsf'
    queue = 'normal'
}
```

See the [Nextflow executor documentation](https://www.nextflow.io/docs/latest/executor.html) for full configuration details.

---

## Local Linux / macOS (Docker)

### 1. Install Docker

**Linux:**
```bash
curl -fsSL https://get.docker.com | sh
sudo usermod -aG docker $USER
# Log out and back in for group changes to take effect
```

**macOS:**
Install [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop/) and allocate at least 32 GB RAM in Settings → Resources.

### 2. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/   # or ~/bin/
```

### 3. Clone and run

```bash
git clone https://github.com/uoozcan/ChampHLA.git
cd ChampHLA

# Test run
nextflow run main.nf -profile test,docker

# Real data — WES example
nextflow run main.nf \
    --input /path/to/bam_files/ \
    --input_type bam \
    --seq_type dna \
    --tools optitype,hlahd,t1k \
    --hlahd_db /path/to/hlahd_db \
    --t1k_hlaidx /path/to/t1k_hla.fa \
    --extract_hla_region \
    --enable_majority_voting \
    --outdir results/ \
    -profile docker
```

### Resource requirements

| Tool panel | Peak RAM | CPUs | Typical time (1 WES sample) |
|---|---|---|---|
| OptiType only | 2 GB | 2 | 1 min |
| OptiType + ArcasHLA + HLA-HD | 8 GB | 4 | 5 min |
| Full 7-tool panel (excl. SpecHLA) | 16 GB | 8 | 15 min |
| Full 8-tool panel | 32 GB | 8 | 20 min |

---

## Windows (WSL2 + Docker Desktop)

ChampHLA runs on Windows through the Windows Subsystem for Linux (WSL2).

### 1. Enable WSL2

Open PowerShell as Administrator:
```powershell
wsl --install -d Ubuntu-22.04
```

Restart your computer when prompted.

### 2. Install Docker Desktop

Download and install [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop/). In Settings:
- Enable **WSL 2 based engine**
- Under Resources → WSL Integration, enable your Ubuntu distribution
- Allocate at least 16 GB RAM (32 GB recommended)

### 3. Install Nextflow inside WSL2

Open your Ubuntu terminal:
```bash
sudo apt update && sudo apt install -y default-jdk
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### 4. Clone and run

```bash
git clone https://github.com/uoozcan/ChampHLA.git
cd ChampHLA

nextflow run main.nf \
    --input /mnt/c/Users/yourname/data/bam_files/ \
    --input_type bam \
    --tools optitype,hlahd,t1k \
    --hlahd_db /path/to/hlahd_db \
    --t1k_hlaidx /path/to/t1k_hla.fa \
    --outdir results/ \
    -profile docker
```

Note: Store input data on the Linux filesystem (e.g., `/home/user/data/`) rather than `/mnt/c/` for better I/O performance.

---

## SpecHLA Installation

SpecHLA is the only tool in ChampHLA that historically lacked a pre-built container. ChampHLA now provides two installation paths:

### Option A: Docker container (recommended)

A Dockerfile is provided at `containers/spechla/Dockerfile`. Build locally:

```bash
cd containers/spechla
docker build -t champhla/spechla:1.1 .
```

For Singularity:
```bash
singularity build spechla.sif docker-daemon://champhla/spechla:1.1
# Move to your cache directory
mv spechla.sif $NXF_SINGULARITY_CACHEDIR/
```

### Option B: Local installation (HPC fallback)

If the container does not work in your environment (e.g., nested Singularity issues), install SpecHLA natively:

```bash
git clone https://github.com/deepomicslab/SpecHLA.git /opt/SpecHLA
cd /opt/SpecHLA
bash install.sh
```

Then configure ChampHLA to use the local installation:
```bash
nextflow run main.nf \
    --tools spechla,optitype,hlahd \
    --use_local_spechla true \
    --spechla_path /opt/SpecHLA \
    -profile singularity
```

### Running without SpecHLA

SpecHLA is optional. Simply omit it from `--tools`:
```bash
--tools optitype,hlahd,t1k,arcashla
```

The remaining seven tools cover all three modalities and provide strong ensemble accuracy without SpecHLA.

---

## Verifying Your Installation

```bash
# Show help and verify Nextflow can parse the pipeline
nextflow run main.nf --help

# Run the built-in test (uses small synthetic fixtures)
nextflow run main.nf -profile test,docker

# Check that containers are accessible
docker images | grep -E "optitype|arcas|hla-hd|t1k"
# or for Singularity:
ls $NXF_SINGULARITY_CACHEDIR/*.sif
```

If the test profile completes without errors, your installation is ready.
