#!/usr/bin/env bash
#
# build_locityper_db.sh — Build a versioned Locityper locus-pangenome database
# for the classical HLA loci, pinned to a specific IPD-IMGT/HLA release.
#
# The database is the haplotype reference Locityper genotypes against. We build
# it from IPD-IMGT/HLA genomic (`*_gen`) sequences so that haplotype names carry
# the IMGT allele id directly (making the haplotype->allele crosswalk trivial).
# Optionally, HPRC pangenome haplotypes can be added for full structural
# diversity; in that case a crosswalk TSV must be supplied to the parser.
#
# A manifest.json recording the IMGT release, reference build, loci and tool
# versions is written alongside the DB so the genotyping run is reproducible.
#
# Usage:
#   build_locityper_db.sh \
#       --out      assets/locityper_db \
#       --ref      /path/to/GRCh38.fa \
#       --imgt-dir /path/to/IMGTHLA \         # clone of ANHIG/IMGTHLA (pinned tag)
#       --imgt-version 3.59.0 \
#       [--loci A,B,C,DRB1,DQA1,DQB1,DPB1] \
#       [--kmer 25] [--threads 8]
#
# Requires: locityper, jellyfish (k-mer counting), samtools, minimap2 on PATH.

set -euo pipefail

OUT=""
REF=""
IMGT_DIR=""
IMGT_VERSION=""
LOCI="A,B,C,DRB1,DQA1,DQB1,DPB1"
KMER=25
THREADS=8

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out)          OUT="$2"; shift 2 ;;
        --ref)          REF="$2"; shift 2 ;;
        --imgt-dir)     IMGT_DIR="$2"; shift 2 ;;
        --imgt-version) IMGT_VERSION="$2"; shift 2 ;;
        --loci)         LOCI="$2"; shift 2 ;;
        --kmer)         KMER="$2"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

for req in OUT REF IMGT_DIR IMGT_VERSION; do
    if [[ -z "${!req}" ]]; then
        echo "ERROR: --${req,,} is required" >&2
        exit 1
    fi
done

for tool in locityper jellyfish samtools minimap2; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found on PATH" >&2; exit 1; }
done

mkdir -p "$OUT"
echo "[build_locityper_db] Output: $OUT  IMGT: $IMGT_VERSION  loci: $LOCI"

# Step 1: build the jellyfish k-mer counts of the reference genome used by
# Locityper for read recruitment / background depth modelling.
JF="${OUT}/genome.jf"
if [[ ! -f "$JF" ]]; then
    echo "[build_locityper_db] Counting ${KMER}-mers in reference..."
    jellyfish count --canonical --lower-count 2 --out-counter-len 2 \
        --mer-len "$KMER" --threads "$THREADS" --size 3G \
        --output "$JF" "$REF"
fi

# Step 2: initialise the Locityper DB against the reference + k-mer counts.
echo "[build_locityper_db] Initialising Locityper DB..."
locityper add --database "$OUT" --reference "$REF" --jellyfish "$JF" || true

# Step 3: add each classical HLA locus from its IPD-IMGT/HLA genomic alignment.
IFS=',' read -r -a LOCUS_ARR <<< "$LOCI"
for locus in "${LOCUS_ARR[@]}"; do
    gen_fasta="${IMGT_DIR}/fasta/${locus}_gen.fasta"
    if [[ ! -f "$gen_fasta" ]]; then
        echo "[build_locityper_db] WARNING: missing ${gen_fasta}; skipping ${locus}" >&2
        continue
    fi
    echo "[build_locityper_db] Adding locus HLA-${locus} from ${gen_fasta}..."
    # `locityper add` aligns the locus haplotype sequences into the DB. Haplotype
    # names are taken from the FASTA headers (IMGT allele ids => trivial crosswalk).
    locityper add --database "$OUT" --reference "$REF" \
        --locus "HLA-${locus}" --sequences "$gen_fasta" \
        --threads "$THREADS" || \
        echo "[build_locityper_db] WARNING: failed to add ${locus}" >&2
done

# Step 4: reproducibility manifest.
MANIFEST="${OUT}/manifest.json"
cat > "$MANIFEST" <<EOF
{
  "imgt_hla_version": "${IMGT_VERSION}",
  "reference": "$(basename "$REF")",
  "loci": "${LOCI}",
  "kmer": ${KMER},
  "built_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "locityper_version": "$(locityper --version 2>/dev/null | head -n1 || echo unknown)",
  "source": "IPD-IMGT/HLA genomic (*_gen) sequences",
  "crosswalk_required": false
}
EOF

echo "[build_locityper_db] Done. Manifest: ${MANIFEST}"
