#!/bin/bash
#SBATCH --job-name=exon_usage_pverr
#SBATCH --output=exon_usage.%j.out
#SBATCH --error=exon_usage.%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --account=pi_hputnam_uri_edu

# Script created in /project/pi_hputnam_uri_edu/dbecks/Becker_RNASeq/scripts
# Change to scratch workspace for fast I/O
cd /scratch4/workspace/danielle_becker_uri_edu-becker_exon_usage

# Paths
WS=/scratch4/workspace/danielle_becker_uri_edu-becker_exon_usage
OUT_DIR=${WS}/exon_usage
GTF=${WS}/refs/Pocillopora_verrucosa_modified_HIv1.gtf
BAM_DIR=/project/pi_hputnam_uri_edu/dbecks/Becker_RNASeq/data/mapped

mkdir -p "${OUT_DIR}"

# Load modules
module load python/3.8.18
module load samtools/1.19.2

echo "Starting exon usage pipeline at $(date)"

# 1) Flatten GTF → exon bins (DEXSeq style)
FLAT_GFF=${OUT_DIR}/Pverr_flattened_exon_bins.gff
if [ ! -f "${FLAT_GFF}" ]; then
    echo "Flattening GTF..."
    dexseq_prepare_annotation.py "${GTF}" "${FLAT_GFF}"
    echo "Flattened GFF created: $(wc -l ${FLAT_GFF}) lines"
fi

# 2) Index BAMs (safe to run multiple times)
echo "Indexing BAMs..."
for BAM in "${BAM_DIR}"/*.bam; do
    if [ ! -f "${BAM}.bai" ]; then
        echo "Indexing ${BAM##*/}"
        samtools index "${BAM}"
    fi
done

# 3) Count reads per exon bin for each sample
echo "Counting reads per exon bin (this takes time)..."
N_BAMS=$(ls "${BAM_DIR}"/*.bam | wc -l)
echo "Processing ${N_BAMS} BAM files..."

for BAM in "${BAM_DIR}"/*.bam; do
    SAMPLE=$(basename "${BAM}" .sam.sorted.bam)
    COUNTS_FILE=${OUT_DIR}/${SAMPLE}.exon_counts.txt

    if [ -f "${COUNTS_FILE}" ]; then
        echo "✓ ${SAMPLE}: counts exist ($(wc -l < "${COUNTS_FILE}") lines)"
        continue
    fi

    echo "Counting for ${SAMPLE}..."
    dexseq_count.py \
        -p yes \                    # paired-end (change to 'no' if single-end)
        -s no \                     # unstranded (change to 'yes'/'reverse' if needed)
        "${FLAT_GFF}" \
        "${BAM}" \
        "${COUNTS_FILE}"

    if [ $? -eq 0 ]; then
        echo "✓ ${SAMPLE}: $(wc -l < "${COUNTS_FILE}") exon bins counted"
    else
        echo "✗ ${SAMPLE}: dexseq_count.py failed"
    fi
done

echo "Pipeline complete at $(date)"
echo "Exon counts saved to: ${OUT_DIR}"
ls -lh ${OUT_DIR}/*.exon_counts.txt | head -5
EOF

chmod +x run_exon_usage_pverr.sh
