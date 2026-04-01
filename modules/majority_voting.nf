/*
 * Consensus workflow wrapper around aggregated calls and weighted voting.
 */

include { AGGREGATE_RESULTS } from './aggregation'

process MAJORITY_VOTING {
    tag "consensus"
    label 'process_low'
    publishDir "${params.outdir}/consensus", mode: 'copy'

    input:
    path calls_tsv

    output:
    path "consensus_calls.tsv", emit: consensus
    path "runtime_weights.json", emit: weights

    script:
    def modality = params.run_modality ?: (params.input_type == 'fastq' ? (params.seq_type == 'rna' ? 'rnaseq' : 'wes') : 'wgs')
    def minWeight = params.consensus_min_weight ?: 0.0
    def useGeneSpecific = params.consensus_use_gene_specific_weights ? '--use-gene-specific' : ''
    def weightFile = params.consensus_weight_file ?: ''
    """
    if [ -n "${weightFile}" ] && [ -f "${weightFile}" ]; then
      cp "${weightFile}" runtime_weights.json
    else
      python3 - << 'PYEOF'
import csv
import json
from collections import defaultdict

modality = "${modality}"
weights = {"tool_weights": defaultdict(dict), "gene_weights": defaultdict(lambda: defaultdict(dict))}

with open("${calls_tsv}", "r", encoding="utf-8") as handle:
    r = csv.DictReader(handle, delimiter="\t")
    seen = defaultdict(set)
    for row in r:
        tool = row.get("tool", "").strip()
        gene = row.get("gene", "").strip()
        if not tool or not gene:
            continue
        seen[tool].add(gene)

for tool, genes in seen.items():
    weights["tool_weights"][tool][modality] = {"final_weight": 1.0}
    for gene in genes:
        weights["gene_weights"][tool][modality][gene] = {"final_weight": 1.0}

# convert defaultdicts
weights["tool_weights"] = dict(weights["tool_weights"])
weights["gene_weights"] = {k: dict(v) for k, v in weights["gene_weights"].items()}

with open("runtime_weights.json", "w", encoding="utf-8") as out:
    json.dump(weights, out, indent=2, sort_keys=True)
PYEOF
    fi

    python3 ${projectDir}/bin/majority_voting.py \
      --calls ${calls_tsv} \
      --weights runtime_weights.json \
      --output consensus_calls.tsv \
      --min-weight ${minWeight} \
      ${useGeneSpecific}
    """
}

workflow MAJORITY_VOTING_WORKFLOW {
    take:
    ch_optitype
    ch_arcashla
    ch_spechla
    ch_hlahd
    ch_polysolver
    ch_kourami
    ch_t1k
    ch_seq2hla

    main:
    ch_all_results = ch_optitype.map { sample_id, f -> f }
        .mix(ch_arcashla.map { sample_id, f -> f })
        .mix(ch_spechla.map { sample_id, f -> f })
        .mix(ch_hlahd.map { sample_id, f -> f })
        .mix(ch_polysolver.map { sample_id, f -> f })
        .mix(ch_kourami.map { sample_id, f -> f })
        .mix(ch_t1k.map { sample_id, f -> f })
        .mix(ch_seq2hla.map { sample_id, f -> f })

    AGGREGATE_RESULTS(ch_all_results.collect())
    MAJORITY_VOTING(AGGREGATE_RESULTS.out.calls)

    emit:
    calls = AGGREGATE_RESULTS.out.calls
    consensus = MAJORITY_VOTING.out.consensus
    weights = MAJORITY_VOTING.out.weights
}
