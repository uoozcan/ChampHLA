# Nextflow Modules

Native execution modules in this directory:

- `optitype.nf`
- `arcashla.nf`
- `spechla.nf`
- `hlahd.nf`
- `polysolver.nf`
- `kourami.nf`
- `t1k.nf`
- `seq2hla.nf`
- `bam_to_fastq.nf`
- `aggregation.nf`
- `majority_voting.nf`

The consensus flow is:
1. Each tool emits `<sample>_<tool>.txt`.
2. `AGGREGATE_RESULTS` builds `aggregated_calls.tsv`.
3. `MAJORITY_VOTING` runs `bin/majority_voting.py` with runtime weights.
