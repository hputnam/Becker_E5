#!/bin/bash
#SBATCH --job-name="KofamScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=100GB
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/KEGG/

echo "Loading modules" $(date)
module load kofam_scan/1.3.0-foss-2019b
module load libyaml/0.1.5
module unload HMMER/3.3.1-foss-2019b
module load HMMER/3.3.2-gompi-2019b
module list

#echo "Starting analysis... downloading KO database" $(date)
#wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz #download KO database
#wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
#gunzip ko_list.gz
#tar xf profiles.tar.gz

echo "Beginning mapping" $(date)
/opt/software/kofam_scan/1.3.0-foss-2019b/exec_annotation -o Pver_KO_annot.tsv -k ./ko_list -p ./profiles/eukaryote.hal -E 0.00001 -f detail-tsv --report-unannotated /data/putnamlab/REFS/Pverr/Pver_proteins_names_v1.0.faa

echo "Analysis complete!" $(date)
