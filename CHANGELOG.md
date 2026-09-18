# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.1.0] - 2026-04-19

### Added
- `--weight-alpha` / `--weight-beta` CLI flags on `hla_benchmark.py` and `run_1000g_benchmark.py`
  to override the default 0.7/0.3 reliability-confidence split at runtime
- `--figures-dir` and `--tables-dir` are now required args in `generate_html_report.py`;
  hard-coded scratch paths and username removed for portability
- Weight formula sensitivity analysis: three variants (1.0/0.0, 0.5/0.5, 0.0/1.0) evaluated
  via `slurm_weight_sensitivity.sh`; results in `analysis/weight_sensitivity/`
- Wilson 95% CIs on all `overall_correct_call_rate` summary columns in benchmark output tables
- `analysis/weight_sensitivity/sensitivity_comparison.tsv`: cross-modality summary of
  WeightedConsensus performance across all weight variants

### Changed
- `bin/generate_html_report.py`: `--figures-dir` and `--tables-dir` are now required (no defaults)

## [3.0.0] - 2026-04-18

### Added
- 8-tool HLA typing ensemble: OptiType, ArcasHLA, SpecHLA, HLA-HD, Kourami, T1K, POLYSOLVER, Seq2HLA
- `read_confidence_v2` confidence weighting: `final_weight = 0.7 × base_reliability + 0.3 × effective_confidence_score`
- Guardrail system per tool × modality: `applied`, `poor_calibration`, `no_confidence` statuses with shrinkage factors
- Tri-modal benchmark (WGS n=99, WES n=51, RNA-seq n=50) with 60/20/20 training/validation/holdout split
- 1000 Genomes NYGC 30× CRAM integration: 49-sample wave-2 WGS batch downloaded from EBI and typed
- Publication-quality figure generation (`bin/generate_figures_v2.py`): 8 figures PNG 300dpi + PDF
- HTML benchmark report v4 (`analysis/generate_report_v4.py`): standalone 3MB report embedding all figures
- IMGT/HLA v3.59.0 pinned for reproducibility
- Abstention tradeoff analysis (`abstention_tradeoff.tsv`): callable rate vs accuracy at configurable thresholds
- Discordance taxonomy (`discordance_summary.tsv`): DNA/RNA discordance, expression bias, technical conflict tags
- SLURM array scripts with `%10` throttle for disk-safe batch processing

### Key benchmark results (holdout)
- WES (n=12): WeightedConsensus 91.7%, MajorityVote 94.4%
- RNA-seq (n=7): WeightedConsensus/MajorityVote 100%, T1K 90.5%
- WGS (n=20): WeightedConsensus 26.7%, callable rate 56.7% at threshold 0.45
- WGS best single tool: OptiType 43.3%

### Changed
- Benchmark scope: `partial_realdata_benchmark` (A/B/C loci; DRB1/DQB1 excluded pending truth source)
- Truth source: Gourraud et al. 2014 1000 Genomes (HLA-A, -B, -C only)
- Benchmark split expanded from 30/10/10 to 60/20/20 to accommodate 100-sample tri-modal cohort

---

## [2.0.0] - 2025-11-21

### Added
- Complete multi-tool HLA typing pipeline with Nextflow DSL2
- Support for three HLA typing tools:
  - OptiType (DNA/RNA-seq)
  - ArcasHLA (RNA-seq optimized, works with DNA)
  - SpecHLA (Exome/WGS with variant calling)
- Majority voting / consensus calling module
- Flexible input support (BAM, CRAM, FASTQ)
- HLA region extraction for space-efficient processing
- Chromosome 6 extraction optimization for SpecHLA
- CSC Puhti HPC profile with SLURM support
- Docker and Singularity container support
- Comprehensive error handling and logging
- Automatic BAM index creation
- Chromosome naming convention detection (chr6 vs 6)
- Exon-only mode for SpecHLA (critical for exome data)
- SLURM submission script with proper parameter handling
- Comprehensive documentation (README, QUICKSTART guide)
- Pipeline execution reports (timeline, trace, DAG)

### Features
- **Smart Input Handling**: Automatically detects and processes different input formats
- **Optimized for Exome Data**: Special exon-only mode for SpecHLA
- **Space Efficient**: Optional HLA region extraction to reduce disk usage
- **Resume Capability**: Nextflow's built-in resume functionality
- **Parallel Processing**: Efficiently processes multiple samples simultaneously
- **Comprehensive Logging**: Detailed logs for each tool and process
- **Resource Management**: Configurable CPU, memory, and time limits
- **Container Ready**: Pre-configured for Docker and Singularity

### Fixed
- SLURM script parameter handling (parameters now properly passed to Nextflow)
- Container paths for CSC Puhti pre-downloaded SIF files
- SpecHLA chromosome naming issues (chr6 vs 6)
- BAM indexing in various scenarios
- Error handling for missing or empty files

### Configuration
- CSC Puhti profile with pre-downloaded containers
- SLURM executor configuration
- Resource labels for different process types
- Automatic retry on failure
- Comprehensive parameter validation

### Documentation
- Complete README with usage examples
- Quick Start guide for new users
- Detailed parameter descriptions
- Troubleshooting section
- CSC Puhti specific instructions
- Performance benchmarks

### Known Issues
- SpecHLA requires exon-only mode (`--spechla_exon_only 1`) for exome data
- Some tools may produce empty results for low-coverage samples
- Chromosome 6 naming must be detected correctly from BAM headers

### Compatibility
- Nextflow: ≥23.04.0
- Docker/Singularity required
- Tested on CSC Puhti (SLURM)
- Compatible with hg19 and hg38 reference builds

## [1.0.0] - 2025-11-XX

### Initial Development
- Basic pipeline structure
- Single tool support
- Initial testing on 1000 Genomes data

---

## Planned for Future Releases

### [2.1.0] - Planned
- [ ] Additional HLA typing tools (HLA-HD, HLA*LA, xHLA)
- [ ] Advanced consensus algorithms
- [ ] Quality score integration
- [ ] Phasing information
- [ ] Interactive HTML reports
- [ ] Database storage support

### [2.2.0] - Planned
- [ ] Cloud deployment (AWS, Google Cloud)
- [ ] GUI interface
- [ ] Real-time monitoring dashboard
- [ ] Automated quality control
- [ ] Population frequency analysis

### [3.0.0] - Planned
- [ ] Machine learning-based consensus
- [ ] Novel allele detection
- [ ] Integration with clinical databases
- [ ] Multi-sample joint calling
