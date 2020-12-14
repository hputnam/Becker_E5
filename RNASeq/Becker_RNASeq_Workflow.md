# Check for required software on Bluewave
### Update software to latest version
- fastqc  
- MultiQC  
- fastp  
- HiSat2  
- Samtools  
- StringTie  
- gffcompare  
- Python  

# 1) Obtain Reference Genome
[Buitrago-López et al 2020](https://academic.oup.com/gbe/article/12/10/1911/5898631)

```
cd /data/putnamlab/REFS/

mkdir Pverr
```

### Genome scaffolds
```
wget http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.fasta.gz
```
### Gene Models (CDS)
```
wget
http://pver.reefgenomics.org/download/Pver_genes_names_v1.0.fna.gz
```
### Gene Models (Proteins)
```
wget
http://pver.reefgenomics.org/download/Pver_proteins_names_v1.0.faa.gz
```
### Gene Models (GFF)
```
wget
http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.gff3.gz
```

### Check to ensure data transfer of genome files
reef genomics md5 checksums

Core files	MD5 hash
Genome scaffolds			fb4d03ba2a9016fabb284d10e513f873
Gene models (CDS)		019359071e3ab319cd10e8f08715fc71
Gene models (proteins)	438f1d59b060144961d6a499de016f55
Gene models (GFF3)		614efffa87f6e8098b78490a5804c857

Miscellaneous files	MD5 hash
Full transcripts	76b5d8d405798d5ca7f6f8cc4b740eb2


On Bluewaves
019359071e3ab319cd10e8f08715fc71  Pver_genes_names_v1.0.fna.gz


### Functional annotation
[functional annotatino xlsx file link](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/gbe/12/10/10.1093_gbe_evaa184/1/evaa184_supplementary_data.zip?Expires=1608305338&Signature=q7Y3SsLOTtlk5ZF6UMMLGy~HRtUFQqRbQt8ZiasnG3EeO11vmHcNgToGSowqYxQK1vibkmPEzMWDeS6u8qG~D20t7G31abz9zbpFrdW9T0cisHAQwY5g~lyK-WRFd-EDYW1eHFI4x~vU0G0xopva7kx1KlXdWxyZW86Fr7CDckFFvav78SAvZtmcvL8WuY4tWmEf33LK4ruuX7ZndqT8k~Kzag57phDdN1qleKWmeAf2wI-Wn8B4w-gV7UU4WQV1Ybs1wwdmexfPxH-DYEuSm-3T4sFq52FW1eRa8WD0V9XDUyysgajGh3sXRxHy-hEUUdnrzlVEk9~Doo9l9IIaUA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)


# make folder structure
```
mkdir data
cd data

```

# Downloaded files

## path where we stored the RAW fastq.gz files
```/data/putnamlab/KITT/hputnam/20201209_Becker_RNASeq_combo/combo```

# 2) Check file integrity 

a) Count all files to make sure all downloaded

```
ls -1 | wc -l
```

b) Verify data transfer integrity with md5sum

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/check_transfer.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw

md5sum /data/putnamlab/KITT/hputnam/20201209_Becker_RNASeq_combo/combo/*.gz > URIcheckmd5.md5

```

```
sbatch ../../scripts/check_transfer.sh 
```

### Checksum from Genewiz
```
/data/putnamlab/KITT/hputnam/20201209_Becker_RNASeq_combo/combo/md5sum_list.txt
```


### Checksum from files on Bluewaves



c) Count number of reads per file 

check for code after @ in fastq.gz files(e.g.,@GWNJ).

```
zcat *.gz | echo $((`wc -l`/4)) > rawread.counts.txt

```



# 3) Run FastQC

a) Make folders for raw FastQC results and scripts

b) Write script for checking quality with FastQC and submit as job on bluewaves

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/fastqc_raw.sh
```

```  
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/

module load all/FastQC/0.11.9-Java-11

for file in /data/putnamlab/KITT/hputnam/20201209_Becker_RNASeq_combo/combo/*.gz
do
fastqc $file --outdir /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw/qc
done
```

```
sbatch /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/fastqc_raw.sh
```


c) Make sure all files were processed

```
/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw/qc ls -1 | wc -l 
```

## Combined QC output into 1 file with MultiQC

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw/qc

```

c) Copy MultiQC files to local computer

```
scp -r hputnam@bluewaves.uri.edu:/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw/qc/*.html /Users/hputnam/MyProjects/Becker_E5/qc

```

The read count is lower than the reads quoted by Genewiz. They are re-sequencing to make up the difference. 


# 4) Trim and clean reads 

a) Make trimmed reads folder in all other results folders 

```

mkdir data/trimmed
cd trimmed

```

c) Write script for Trimming and run on bluewaves

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/trim.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw

module load fastp/0.19.7-foss-2018b

for file in "C17" "C18" "C19" "C20" "C21" "C22" "C23" "C24" "C25" "C26" "C27" "C28" "C29" "C30" "C31" "C32" \
"E1" "E2" "E3" "E4" "E5" "E6" "E7" "E8" "E9" "E10" "E11" "E12" "E13" "E14" "E15" "E16" \
do
fastp --in1 /data/putnamlab/KITT/hputnam/20201124_Becker_RNASeq/${file}_R1_001.fastq.gz --in2 /data/putnamlab/KITT/hputnam/20201124_Becker_RNASeq/${file}_R2_001.fastq.gz --out1 /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/trimmed/${file}_R1_001.clean.fastq.gz --out2 /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/trimmed/${file}_R2_001.clean.fastq.gz --qualified_quality_phred 20 --unqualified_percent_limit 10 --length_required 100 detect_adapter_for_pe --cut_right cut_right_window_size 5 cut_right_mean_quality 20
done
```
```
sbatch /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/trim.sh
```


# 5) Check quality of trimmed files 

a) Check number of files 

```
/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/trimmed/ ls -1 | wc -l
```

b) Check number of reads

```
zgrep -c "@GWNJ" /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/trimmed/*.gz > trimmed_seq_counts


c) Run FastQC on trimmed data

mkdir trimmed_qc

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/fastqc_trimmed.sh
```

```  
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw

module load all/FastQC/0.11.9-Java-11

for file in /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/trimmed/*.gz
do
fastqc $file --outdir /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/trimmed/trimmed_qc
done
```

```
sbatch /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/fastqc_trimmed.sh
```

fastqc E16_R*_001.clean.fastq.gz --outdir /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/trimmed/trimmed_qc


d) Run MultiQC on trimmed data 

scp -r hputnam@bluewaves.uri.edu:/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/trimmed/trimmed_qc/*.html /Users/hputnam/MyProjects/Becker_E5/trimmed_qc

/Users/hputnam/MultiQC/scripts/multiqc .

# 6) Align reads 

a) Generate genome build

### Need to unzip genome files before running

```
gunzip Pver_genome_assembly_v1.0.fasta.gz  
gunzip Pver_genome_assembly_v1.0.gff3.gz
```

### HiSat2 Align reads to refernece genome
[HiSat2](https://daehwankimlab.github.io/hisat2/main/)
[HiSat2 Github](https://github.com/DaehwanKimLab/hisat2)

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/Hisat2_genome_build.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw

module load HISAT2/2.1.0-foss-2018b

hisat2-build -f /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0.fasta ./Pver_ref

```
sbatch /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/Hisat2_genome_build.sh
```

b) Align reads to genome

mkdir mapped

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/Hisat2_align2.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

module load HISAT2/2.1.0-foss-2018b

for i in "C17" "C18" "C19" "C20" "C21" "C22" "C23" "C24" "C25" "C26" "C27" "C28" "C29" "C30" "C31" "C32" "E1" "E2" "E3" "E4" "E5" "E6" "E7" "E8" "E9" "E10" "E11" "E12" "E13" "E14" "E15" "E16"
do
hisat2 -p 48 --rna-strandness RF --dta -q -x /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw/Pver_ref -1 /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/trimmed/${i}_R1_001.clean.fastq.gz \
-2 /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/trimmed/${i}_R2_001.clean.fastq.gz -S /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/${i}.sam
done
```

```
sbatch /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/Hisat2_align2.sh
```

## Sort and convert sam to bam

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/SAMtoBAM.sh
```

/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/*.sam

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

module load SAMtools/1.9-foss-2018b

for i in "C17" "C18" "C19" "C20" "C21" "C22" "C23" "C24" "C25" "C26" "C27" "C28" "C29" "C30" "C31" "C32" "E1" "E2" "E3" "E4" "E5" "E6" "E7" "E8" "E9" "E10" "E11" "E12" "E13" "E14" "E15" "E16"
do
samtools view -b -S /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/${i}.sam | samtools sort -o /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/${i}.sorted.bam -O bam
done
```

```
sbatch /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/SAMtoBAM.sh
```

### Remove Sam files to save space
```
rm /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/*.sam

```


# 7) Perform gene counts with stringTie

```
mkdir counts
cd counts

```


b) Assemble and estimate reads 

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/StringTie_Assemble.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

module load StringTie/1.3.5-foss-2018b


for i in "C17" "C18" "C19" "C20" "C21" "C22" "C23" "C24" "C25" "C26" "C27" "C28" "C29" "C30" "C31" "C32" "E1" "E2" "E3" "E4" "E5" "E6" "E7" "E8" "E9" "E10" "E11" "E12" "E13" "E14" "E15" "E16"
do
stringtie /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/${i}.sorted.bam -p 48 -e -G /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0.gff3 -o /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/${i}.gtf 
done
```

```
sbatch /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/StringTie_Assemble.sh
```



c) Merge stringTie gtf results 

```
ls *gtf > mergelist.txt

```
 nano samplelist.txt
 
 C17	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C17.gtf
 C18	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C18.gtf
 C19	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C19.gtf
 C20	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C20.gtf
 C21	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C21.gtf
 C22	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C22.gtf
 C23	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C23.gtf
 C24	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C24.gtf
 C25	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C25.gtf
 C26	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C26.gtf
 C27	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C27.gtf
 C28	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C28.gtf
 C29	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C29.gtf
 C30	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C30.gtf
 C31	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C31.gtf
 C32	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/C32.gtf
 E1		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E1.gtf
 E2		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E2.gtf
 E3		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E3.gtf
 E4		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E4.gtf
 E5		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E5.gtf
 E6		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E6.gtf
 E7		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E7.gtf
 E8		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E8.gtf
 E9		/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E9.gtf
 E10	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E10.gtf  
 E11	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E11.gtf     
 E12	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E12.gtf    
 E13	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E13.gtf    
 E14	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E14.gtf    
 E15	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E15.gtf  
 E16	/data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/E16.gtf    
   
f) Create gene matrix

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/GTFtoCounts.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

module load StringTie/1.3.5-foss-2018b
module load Python/2.7.15-foss-2018b


python prepDE.py  -i samplelist.txt -g Poc_gene_count_matrix.csv
```

```
sbatch /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/GTFtoCounts.sh
```


g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/GTF_merge/gene_count_ofav_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/
```

