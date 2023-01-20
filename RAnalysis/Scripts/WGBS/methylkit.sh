#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=danielle_becker@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/
#SBATCH --error="script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script" #once your job is completed, any final job report comments will be put in this file

echo "Loading modules" $(date)
module load msPIPE/20220530-foss-2021b
Rscript methyl.R

file.list=list("/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/1.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/2.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/3.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/4.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/5.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/6.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/7.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/8.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/9.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/10.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/11.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/12.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/13.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/14.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/15.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/16.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/17.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/18.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/19.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/20.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/21.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/22.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/23.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/24.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/25.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/26.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/27.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/28.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/29.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/30.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/31.CpG_report.txt",
                  "/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/32.CpG_report.txt")

                  myobj = methRead(file.list,
                              sample.id=list("enriched1", "enriched2", "enriched3", "enriched4", "enriched5", "enriched6", "enriched7", "enriched8", "enriched9", "enriched10", "enriched11", "enriched12", "enriched13", "enriched14", "enrich  ed15", "enriched16", "control17", "control18", "control19", "control20", "control21", "control22", "control23", "control24", "control25", "control26", "control27", "control28", "control29", "control30", "control31", "control32"),
                              assembly="hg18",treatment=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), pipeline='bismarkCytosineReport', dbtype="tabix")
