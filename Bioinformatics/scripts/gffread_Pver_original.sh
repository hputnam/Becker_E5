#!/bin/bash
#SBATCH --time 06:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --export=NONE
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

cd /work/pi_hputnam_uri_edu/refs/REFS/Pverr

# load modules needed
module load gffread/0.12.7

# Combined code, want to see if these are any different
gffread -E /work/pi_hputnam_uri_edu/refs/REFS/Pverr/Pver_genome_assembly_v1.0.gff3 -T -o /work/pi_hputnam_uri_edu/refs/REFS/Pverr/Pocillopora_verrucosa_original_HIv1.gtf

echo "Check how many fields are in each row of the file; currently there are rows with two different lenghts: 10 and 12"
awk '{print NF}' Pocillopora_verrucosa_original_HIv1.gtf | sort -u

# Use awk to add "gene_id = TRANSCRIPT_ID" to each of the rows that only have a transcript id listed (the non-transcript lines)

awk 'BEGIN {FS=OFS="\t"} {if ($9 ~ /transcript_id/ && $9 !~ /gene_id/) {match($9, /transcript_id "([^"]+)";/, a); $9 = $9 " gene_id \"" a[1] "\";"}; print}' Pocillopora_verrucosa_original_HIv1.gtf > Pocillopora_verrucosa_original_HIv1_modified.gtf

echo "Check how many fields are in each row of the file; Now all the rows should be the same length and only one value should be printed, 12"
awk '{print NF}' Pocillopora_verrucosa_original_HIv1_modified.gtf | sort -u

# remove the non-modified file
rm Pocillopora_verrucosa_original_HIv1.gtf

# rename the modified file
mv Pocillopora_verrucosa_original_HIv1_modified.gtf Pocillopora_verrucosa_original_HIv1.gtf
