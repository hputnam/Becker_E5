#!/bin/bash
#SBATCH --job-name=exon_usage_pverr
#SBATCH --output=exon_usage.%j.out
#SBATCH --error=exon_usage.%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --account=pi_hputnam_uri_edu
#SBATCH --constraint=avx512

cd /scratch4/workspace/danielle_becker_uri_edu-becker_exon_usage

WS=/scratch4/workspace/danielle_becker_uri_edu-becker_exon_usage
OUT_DIR=${WS}/exon_usage
GTF=${WS}/refs/Pocillopora_verrucosa_modified_HIv1.gtf
BAM_DIR=/project/pi_hputnam_uri_edu/dbecks/Becker_RNASeq/data/mapped

mkdir -p "${OUT_DIR}"

# Load modules (Rsubread/featureCounts comes with R-bundle-Bioconductor)
module load uri/main R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1
module load samtools/1.19.2

echo "Starting exon usage pipeline at $(date)"

# 1) Index BAMs
echo "Indexing BAMs..."
for BAM in "${BAM_DIR}"/*.bam; do
    if [ ! -f "${BAM}.bai" ]; then
        echo "Indexing ${BAM##*/}"
        samtools index "${BAM}"
    fi
done

# 2) featureCounts via Rscript (DEXSeq compatible output)
BAM_LIST=$(echo "${BAM_DIR}"/*.bam | tr ' ' '\n' | paste -s -d ' ')
echo "Processing BAMs: ${BAM_LIST}"

Rscript -e "
library(Rsubread)

# List of BAM files
bams <- c(${BAM_LIST//./\"\\\"\n\"\\\"} )

# Run featureCounts
fc <- featureCounts(
    files = bams,
    annot.ext = '${GTF}',
    isGTFAnnotationFile = TRUE,
    GTF.featureType = 'exon',
    GTF.attrType = 'gene_id',  # or 'gene_name' if your GTF uses that
    isPairedEnd = TRUE,
    strandSpecific = 0,
    nthreads = 4,
    countMultiMappingReads = FALSE
)

# Save counts matrix (DEXSeq ready)
counts <- fc\$counts
write.table(counts, file='${OUT_DIR}/Pverr_exon_counts.txt', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)

cat('featureCounts complete. Output shape:', dim(counts), '\n')
cat('Saved to: ${OUT_DIR}/Pverr_exon_counts.txt\n')
"

echo "Pipeline complete at $(date)"
echo "Exon counts saved to: ${OUT_DIR}/Pverr_exon_counts.txt"
head -20 "${OUT_DIR}/Pverr_exon_counts.txt"
ls -lh "${OUT_DIR}"/*.txt
