# Run on URI HPC andromeda in a .sh script for BP as tthe processing time was too long on desktop in RStudio
# sbatch goMWU_job.sh
# Submitted batch job 333229 on 20240726 at 14:07

#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=danielle_becker@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

cd /data/putnamlab/dbecks/Becker_E5/goMWU
module load R/4.2.2-foss-2022b
Rscript GO_MWU_BP.R
