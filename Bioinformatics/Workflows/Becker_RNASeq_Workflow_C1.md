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

# 1) Obtain Reference Genome (add how sym_c1 genome obtained here)

[Liu et al 2018](https://www.nature.com/articles/s42003-018-0098-3)
[Cladocopium goreaui - Clade C1](http://symbs.reefgenomics.org/download/)

```
cd /data/putnamlab/REFS/

mkdir Sym_C1
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

```
Core files	MD5 hash

Genome scaffolds    fb4d03ba2a9016fabb284d10e513f873
Gene models (CDS)   019359071e3ab319cd10e8f08715fc71
Gene models (proteins)    438f1d59b060144961d6a499de016f55
Gene models (GFF3)    614efffa87f6e8098b78490a5804c857

Miscellaneous files	MD5 hash
Full transcripts	            76b5d8d405798d5ca7f6f8cc4b740eb2

On Bluewaves
Pver_genes_names_v1.0.fna.gz  019359071e3ab319cd10e8f08715fc71 
```

### Functional annotation
[functional annotation xlsx file link](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/gbe/12/10/10.1093_gbe_evaa184/1/evaa184_supplementary_data.zip?Expires=1608305338&Signature=q7Y3SsLOTtlk5ZF6UMMLGy~HRtUFQqRbQt8ZiasnG3EeO11vmHcNgToGSowqYxQK1vibkmPEzMWDeS6u8qG~D20t7G31abz9zbpFrdW9T0cisHAQwY5g~lyK-WRFd-EDYW1eHFI4x~vU0G0xopva7kx1KlXdWxyZW86Fr7CDckFFvav78SAvZtmcvL8WuY4tWmEf33LK4ruuX7ZndqT8k~Kzag57phDdN1qleKWmeAf2wI-Wn8B4w-gV7UU4WQV1Ybs1wwdmexfPxH-DYEuSm-3T4sFq52FW1eRa8WD0V9XDUyysgajGh3sXRxHy-hEUUdnrzlVEk9~Doo9l9IIaUA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)


# make folder structure
```
#add the Sym_C1 folder from putnamlab REFS directory to sub folder, refs in Becker_E5 folder

ln -s ../../../../../REFS/Sym_C1/ ./refs/

```

# Downloaded files

## path where we stored the RAW fastq.gz files

```
# these data had already been downloaded to the Becker_E5 - RNA_seq folder from mapping reads to host genome using this command: 

ln -s ../../../../../KITT/hputnam/20201209_Becker_RNASeq_combo/combo/*.fastq.gz ./raw/

/data/putnamlab/KITT/hputnam/20201209_Becker_RNASeq_combo/combo

```

# 2) Check file integrity 

a) Count all files to make sure all downloaded

```
ls -1 | wc -l
```

b) Verify data transfer integrity with md5sum

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/check_transfer.sh
```

```
#!/bin/bash ###creating slurm script
#SBATCH -t 24:00:00 ###give script 24 hours to run
#SBATCH --nodes=1 --ntasks-per-node=1 ###on server, different nodes you can use for processing power, node just do one task 
#SBATCH --export=NONE 
#SBATCH --mem=100GB ###in server allocate 100GB amount of memory
#SBATCH --account=putnamlab ###primary account
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/raw ###path


md5sum /dbecks/Becker_E5/Becker_RNASeq/data/*.gz > URIcheckmd5.md5

md5sum /data/putnamlab/KITT/hputnam/20201209_Becker_RNASeq_combo/combo/*.gz > URIcheckmd5.md5

```

```
sbatch check_transfer.sh 
```
###Submitted batch job 1816196 20201230


### Checksum from Genewiz

```
317dd03e9704a73347c5ccabb86d1e18  C17_R1_001.fastq.gz
b6368b2789cfacc8e8adbadf45d6c533  C17_R2_001.fastq.gz
adb1af7c62b309f7100117f50ebe846d  C18_R1_001.fastq.gz
776842109ee56cd4619d5dec08fe005e  C18_R2_001.fastq.gz
a6d48537ec697d2040798a302bfa0aa1  C19_R1_001.fastq.gz
9b46b63800d83940b53196dccd4eda53  C19_R2_001.fastq.gz
e2a63dca87aa53b4756d51edb02b377f  C20_R1_001.fastq.gz
b2752d8b60c97f754cc24ededc1bf053  C20_R2_001.fastq.gz
f659265a858c029bfb4d2764da76701a  C21_R1_001.fastq.gz
712ffc2eb87aadd42d50f9c10c7c9b81  C21_R2_001.fastq.gz
3fe2c3f4c39c06adca53440c037d48fe  C22_R1_001.fastq.gz
ff775c4de3ad8a585f64367cc78e467d  C22_R2_001.fastq.gz
2f4f97aa4dd101b5a14ab0a4284f0816  C23_R1_001.fastq.gz
aabac96314290241fda7e7413d716b2c  C23_R2_001.fastq.gz
94fa7191998ac59d2e89851b5fb15431  C24_R1_001.fastq.gz
921198dae25a3c4515be0d534281a685  C24_R2_001.fastq.gz
feb863a003dfff84ca3868ec5b674e02  C25_R1_001.fastq.gz
3ae4179cd04cdf04bd5873ba5d6a0534  C25_R2_001.fastq.gz
225bafe67a42557b7a34ac201cf88273  C26_R1_001.fastq.gz
15aee2823967a782a859c039610b0b23  C26_R2_001.fastq.gz
82ec398d754d52fef140a325c80ee289  C27_R1_001.fastq.gz
8aebf98f2af591cf8a1149f263f6d86d  C27_R2_001.fastq.gz
1bf51f4961df4b7da8f284cd7202f6d3  C28_R1_001.fastq.gz
1b4cc1e1be421387abfd7da8d490235b  C28_R2_001.fastq.gz
f02e752d98e9a0168286becf52203c0f  C29_R1_001.fastq.gz
fd280e79ecaa4a51a1670293c742d5fd  C29_R2_001.fastq.gz
27153ad98aa8735c3474bd72dfe041cd  C30_R1_001.fastq.gz
8d977f9ebebea4c26c0cfc9bf7033a9b  C30_R2_001.fastq.gz
b6a2cf9d02329c3c73c4ec02c1441660  C31_R1_001.fastq.gz
52cdaf65c3cb976c168478678131066f  C31_R2_001.fastq.gz
cdfc0509137d2b8624a236b2ded04a47  C32_R1_001.fastq.gz
9148b0e6bbc9b7520c7e1c6caae13fc8  C32_R2_001.fastq.gz
8854b4f07a80e5c0e49a136bb5c76e75  E10_R1_001.fastq.gz
1cccd8da94d4a8a79f1681edc9349166  E10_R2_001.fastq.gz
e7277e7ac17b4643e748e1a3e00dfb04  E11_R1_001.fastq.gz
3ff35b0852ea83593dd438dcbb2deddd  E11_R2_001.fastq.gz
50ab9cee762e076b5183de2c8fc66f02  E12_R1_001.fastq.gz
e188fc0f1f0ecd1feb581dd01574667c  E12_R2_001.fastq.gz
af969bcc9b91dfd59c48f8cb87fb6317  E13_R1_001.fastq.gz
7759be7b6477a401e4ffcf385b80f51e  E13_R2_001.fastq.gz
c0d20187ee3403b1a1e4687afbd2fd91  E14_R1_001.fastq.gz
80f122aca444d7ebdfc153d73666cf4b  E14_R2_001.fastq.gz
a7632155416f4de7c80f0f1466c21457  E15_R1_001.fastq.gz
9eec50018582ad56cc2ad05970da7b90  E15_R2_001.fastq.gz
1f782733b62f17958a813d09793922f0  E16_R1_001.fastq.gz
91f58b438ae17444c510a6636bf244ce  E16_R2_001.fastq.gz
bf0038545deef9b4729b8b322ab5f33b  E1_R1_001.fastq.gz
f3e8fab360a0986397ad1e6dd6a85afd  E1_R2_001.fastq.gz
13cc2a96b34654c62eae372ebd8109ea  E2_R1_001.fastq.gz
a55437f2ba10eac88e70994622c20cc2  E2_R2_001.fastq.gz
be0aa42447cfa4e52e47e4c1587f07f6  E3_R1_001.fastq.gz
0dd2317421ce84e9ec1d9a029752a9b3  E3_R2_001.fastq.gz
44b7219ff47447522e24bb6d10ec871c  E4_R1_001.fastq.gz
6a9249b714e0a1596bad222d94bcca09  E4_R2_001.fastq.gz
beedc3c1f5d7918166f8f26d867cdb12  E5_R1_001.fastq.gz
f45e1969d25e3e07b993672b3a6ef8f7  E5_R2_001.fastq.gz
3792a2a487030d0b3c4ad7fe95bff398  E6_R1_001.fastq.gz
004eb857e060728c619262f7287ec273  E6_R2_001.fastq.gz
fd8d4df782b0142bea08adf04edbe835  E7_R1_001.fastq.gz
c21ffd9c4d24181ff6f5e8bf02c83a95  E7_R2_001.fastq.gz
0852370e047f41a248ca7e7014ad88dc  E8_R1_001.fastq.gz
b6e5da4446b5bccb63e46661e9b1e293  E8_R2_001.fastq.gz
1ecdc2f728e6ce233c9495a2c3bb97f5  E9_R1_001.fastq.gz
be19f114b314fdc30dd39962aa12a3dc  E9_R2_001.fastq.gz
```

c) Cross-reference the checksum document from GENEWIZ with the data we have on our computer

```
### with a small amount of files, able to first cross-check that the sequences matched between both files on the desktop
### used the code below in terminal to cross-check the files and compare for sanity check

```

d) Count number of reads per file using the code after @ in fastq.gz files (e.g.,@GWNJ).

```
zgrep -c "@GWNJ" *.gz > raw_seq_counts

```

```
C17_R1_001.fastq.gz:25158606
C17_R2_001.fastq.gz:25158606
C18_R1_001.fastq.gz:22733345
C18_R2_001.fastq.gz:22733345
C19_R1_001.fastq.gz:24846067
C19_R2_001.fastq.gz:24846067
C20_R1_001.fastq.gz:24030431
C20_R2_001.fastq.gz:24030431
C21_R1_001.fastq.gz:16484060
C21_R2_001.fastq.gz:16484060
C22_R1_001.fastq.gz:22990550
C22_R2_001.fastq.gz:22990550
C23_R1_001.fastq.gz:20905338
C23_R2_001.fastq.gz:20905338
C24_R1_001.fastq.gz:22578178
C24_R2_001.fastq.gz:22578178
C25_R1_001.fastq.gz:29417106
C25_R2_001.fastq.gz:29417106
C26_R1_001.fastq.gz:23267238
C26_R2_001.fastq.gz:23267238
C27_R1_001.fastq.gz:24990687
C27_R2_001.fastq.gz:24990687
C28_R1_001.fastq.gz:23396439
C28_R2_001.fastq.gz:23396439
C29_R1_001.fastq.gz:17900262
C29_R2_001.fastq.gz:17900262
C30_R1_001.fastq.gz:26361873
C30_R2_001.fastq.gz:26361873
C31_R1_001.fastq.gz:24642131
C31_R2_001.fastq.gz:24642131
C32_R1_001.fastq.gz:27076800
C32_R2_001.fastq.gz:27076800
E10_R1_001.fastq.gz:29131006
E10_R2_001.fastq.gz:29131006
E11_R1_001.fastq.gz:18568153
E11_R2_001.fastq.gz:18568153
E12_R1_001.fastq.gz:27805087
E12_R2_001.fastq.gz:27805087
E13_R1_001.fastq.gz:24455094
E13_R2_001.fastq.gz:24455094
E14_R1_001.fastq.gz:22630044
E14_R2_001.fastq.gz:22630044
E15_R1_001.fastq.gz:23796710
E15_R2_001.fastq.gz:23796710
E16_R1_001.fastq.gz:29523059
E16_R2_001.fastq.gz:29523059
E1_R1_001.fastq.gz:25368402
E1_R2_001.fastq.gz:25368402
E2_R1_001.fastq.gz:23770610
E2_R2_001.fastq.gz:23770610
E3_R1_001.fastq.gz:25641209
E3_R2_001.fastq.gz:25641209
E4_R1_001.fastq.gz:21751368
E4_R2_001.fastq.gz:21751368
E5_R1_001.fastq.gz:16381619
E5_R2_001.fastq.gz:16381619
E6_R1_001.fastq.gz:24937261
E6_R2_001.fastq.gz:24937261
E7_R1_001.fastq.gz:24020166
E7_R2_001.fastq.gz:24020166
E8_R1_001.fastq.gz:23675842
E8_R2_001.fastq.gz:23675842
E9_R1_001.fastq.gz:25068848
E9_R2_001.fastq.gz:25068848
=======
/data/putnamlab/KITT/hputnam/20201209_Becker_RNASeq_combo/combo/md5sum_list.txt
```

# 3) Run FastQC

a) Make folders for raw FastQC results and scripts

b) Write script for checking quality with FastQC and submit as job on bluewaves

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/fastqc_raw.sh
```

```  
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=danielle_becker@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/raw
#SBATCH --error="script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script" #once your job is completed, any final job report comments will be put in this file

module load FastQC/0.11.8-Java-1.8

for file in /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/raw/*.gz
do
fastqc $file --outdir /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/raw/qc
done
```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/fastqc_raw.sh
```

### Submitted batch job 1816766 on 20210104



c) Make sure all files were processed

```
ls -1 | wc -l 
#64
```

## Combined QC output into 1 file with MultiQC

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/raw/qc

```

c) Copy MultiQC files to local computer

```
scp -r danielle_becker@bluewaves.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/raw/qc/*.html /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RNASeq/qc

```


# 4) Trim and clean reads 

a) Trimmed reads were already processed in the trimmed folder from previous workflow mapping to the host genome, copied trimmed directory contents to a new directory for C1 processing

```

mkdir data/trimmed_C1
cd trimmed_C1

cp -R /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/ .

```

### See below steps for how reads were trimmed and cleaned in previous workflow for host mapping

c) Write script for Trimming and run on bluewaves

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/trim.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=danielle_becker@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/raw
#SBATCH --error="script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script" #once your job is completed, any final job report comments will be put in this file

module load fastp/0.19.7-foss-2018b

array1=($(ls *.fastq.gz)) #Make an array of sequences to trim
for i in ${array1[@]}; do #Make a loop that trims each file in the array
fastp --in1 ${i} --in2 $(echo ${i}|sed s/_R1/_R2/) --out1 ../trimmed/${i} --out2 ../trimmed/$(echo ${i}|sed s/_R1/_R2/) --qualified_quality_phred 20 --unqualified_percent_limit 10 --length_required 100 --detect_adapter_for_pe --cut_right cut_right_window_size 5 cut_right_mean_quality 20
done

```
```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/trim.sh
Submitted batch job 1819574
```


# 5) Check quality of trimmed files 

a) Check number of files 

```
ls -1 | wc -l 
#64
```

b) Check number of reads

```
zgrep -c "@GWNJ" *.gz > trimmed_seq_counts

C21_R1_001.fastq.gz:12569436
C21_R2_001.fastq.gz:12569436
C23_R1_001.fastq.gz:15920237
C18_R1_001.fastq.gz:17448563
C23_R2_001.fastq.gz:15920237
C22_R1_001.fastq.gz:17783162
C20_R1_001.fastq.gz:18346189
C22_R2_001.fastq.gz:17783162
C26_R2_001.fastq.gz:17679573
C18_R2_001.fastq.gz:17448563
C24_R1_001.fastq.gz:17244418
C26_R1_001.fastq.gz:17679573
C19_R2_001.fastq.gz:18977336
C24_R2_001.fastq.gz:17244418
C20_R2_001.fastq.gz:18346189
C19_R1_001.fastq.gz:18977336
C17_R1_001.fastq.gz:19053464
C17_R2_001.fastq.gz:19053464
C25_R1_001.fastq.gz:22683683
C25_R2_001.fastq.gz:22683683
C29_R1_001.fastq.gz:13714128
C27_R1_001.fastq.gz:18756419
C29_R2_001.fastq.gz:13714128
C27_R2_001.fastq.gz:18756419
C28_R1_001.fastq.gz:18131021
E11_R1_001.fastq.gz:14185860
E11_R2_001.fastq.gz:14185860
C28_R2_001.fastq.gz:18131021
C31_R1_001.fastq.gz:18796208
C30_R1_001.fastq.gz:20242216
C31_R2_001.fastq.gz:18796208
C30_R2_001.fastq.gz:20242216
C32_R1_001.fastq.gz:20583523
C32_R2_001.fastq.gz:20583523
E10_R1_001.fastq.gz:22441661
E12_R1_001.fastq.gz:21608172
E10_R2_001.fastq.gz:22441661
E13_R1_001.fastq.gz:18961264
E12_R2_001.fastq.gz:21608172
E13_R2_001.fastq.gz:18961264
E14_R1_001.fastq.gz:17170729
E14_R2_001.fastq.gz:17170729
E15_R1_001.fastq.gz:17920244
E15_R2_001.fastq.gz:17920244
E1_R1_001.fastq.gz:19467393
E5_R1_001.fastq.gz:12558989
E1_R2_001.fastq.gz:19467393
E4_R1_001.fastq.gz:16639960
E16_R1_001.fastq.gz:22697805
E5_R2_001.fastq.gz:12558989
E2_R2_001.fastq.gz:18593178
E4_R2_001.fastq.gz:16639960
E3_R1_001.fastq.gz:19740214
E2_R1_001.fastq.gz:18593178
E16_R2_001.fastq.gz:22697805
E3_R2_001.fastq.gz:19740214
E7_R1_001.fastq.gz:18431790
E6_R1_001.fastq.gz:18968203
E6_R2_001.fastq.gz:18968203
E7_R2_001.fastq.gz:18431790
E8_R1_001.fastq.gz:17970740
E8_R2_001.fastq.gz:17970740
E9_R1_001.fastq.gz:19333946
E9_R2_001.fastq.gz:19333946


```

c) Run FastQC on trimmed data
```
mkdir trimmed_qc

```
```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/fastqc_trimmed.sh
```

```  
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/trimmed_qc
#SBATCH --error="script_error" 
#SBATCH --output="output_script" 

module load FastQC/0.11.8-Java-1.8

for file in /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/*.gz
do
fastqc $file --outdir /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/trimmed_qc
done
```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/fastqc_trimmed.sh
Submitted batch job 1834516
```


d) Run MultiQC on trimmed data 
```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/trimmed_qc
```
```
scp -r danielle_becker@bluewaves.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/trimmed_qc/*.html /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RNASeq/trimmed_qc

```

# 6) Align reads 

a) Generate genome build with Symbiodiniaceae C1 genome

### HiSat2 Align reads to refernece genome
[HiSat2](https://daehwankimlab.github.io/hisat2/main/)
[HiSat2 Github](https://github.com/DaehwanKimLab/hisat2)


```
# made a new directory for my scripts for Sym_C1 mapping

nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/Hisat2_genome_build_C1.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Sym_C1
#SBATCH --error="script_error" 
#SBATCH --output="output_script"

module load HISAT2/2.1.0-foss-2018b

hisat2-build -f /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Sym_C1/SymbC1.Genome.Scaffolds.fasta ./Sym_C1_ref

```
```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/Hisat2_genome_build_C1.sh
Submitted batch job 1912620
```

b) Align reads to genome

```
mkdir mapped_C1
```
```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/Hisat2_align2_C1.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed_C1
#SBATCH --cpus-per-task=3
#SBATCH --error="script_error" 
#SBATCH --output="output_script"


module load HISAT2/2.1.0-foss-2018b

#Aligning paired end reads
#Has the R1 in array1 because the sed in the for loop changes it to an R2. SAM files are of both forward and reverse reads

array1=($(ls *_R1_001.fastq.gz)) 
for i in ${array1[@]}; do 
hisat2 -p 48 --rna-strandness RF --dta -q -x /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Sym_C1_ref -1 /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed_C1/${i} \
-2 /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed_C1/$(echo ${i}|sed s/_R1/_R2/) -S /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/${i}.sam
done
```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/Hisat2_align2_C1.sh
Submitted batch job 1912654

```

## Sort and convert sam to bam

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/SAMtoBAM_C1.sh

#There will be lots of .tmp file versions in your folder, this is normal while this script runs and they should delete at the end to make one sorted.bam file
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1
#SBATCH --cpus-per-task=3
#SBATCH --error="script_error" 
#SBATCH --output="output_script"

module load SAMtools/1.9-foss-2018b

array1=($(ls *.sam))  
for i in ${array1[@]}; do
samtools sort -o /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/${i}.sorted.bam /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/${i}  
done
```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/SAMtoBAM_C1.sh
Submitted batch job 1912675

```

### Remove Sam files to save space
```
rm /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/*.sam

```


### Check number of mapped reads

Explanation for SAMtools functions for checking the mapped reads in a paired-end dataset found here: https://www.biostars.org/p/138116/

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/SAMtools_mapped.sh

```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1
#SBATCH --cpus-per-task=3
#SBATCH --error="script_error" 
#SBATCH --output="output_script"

module load SAMtools/1.9-foss-2018b

array1=($(ls *.bam))  
for i in ${array1[@]}; do
samtools view -F 0x4 /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/${i} | cut -f 1 | sort | uniq | wc -l > mapped_read_counts_C1
done
```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/SAMtools_mapped.sh
Submitted batch job 1912697
```



# 7) Perform gene counts with stringTie


b) Assemble and estimate reads 

### makes count directory inside script

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/StringTie_Assemble_C1.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1
#SBATCH --cpus-per-task=3

module load StringTie/2.1.4-GCC-9.3.0

array1=($(ls *.bam)) 
for i in ${array1[@]}; do 
stringtie -p 48 --rf -e -G /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Sym_C1/SymbC1.Gene_Models.GFF3 -o /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts_C1/${i}.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/${i}
done
```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/StringTie_Assemble_C1.sh
Submitted batch job 1912768

```

c) Merge stringTie gtf results 

#in this step we are making a file with all the gtf names and stringtie will merge them all together for a master list for your specific genes

```
ls *gtf > mergelist.txt_C1
cat mergelist.txt_C1

module load StringTie/2.1.4-GCC-9.3.0

stringtie --merge -p 8 -G /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Sym_C1/SymbC1.Gene_Models.GFF3 -o stringtie_merged.gtf mergelist.txt_C1

```

d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Sym_C1/SymbC1.Gene_Models.GFF3 -o merged stringtie_merged.gtf

```

e) Re-estimate assembly

```
nano re_estimate.assembly_C1.sh
```
```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1
#SBATCH --cpus-per-task=3

module load StringTie/2.1.4-GCC-9.3.0


array1=($(ls *.bam)) 
for i in ${array1[@]}; do 
stringtie -e -G /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Sym_C1/SymbC1.Gene_Models.GFF3 -o ${i}.merge.gtf ${i}
echo "${i}"
done

```
```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/re_estimate.assembly_C1.sh
Submitted batch job 1913016
```
```
# move merged GTF files to their own folder 
mv *merge.gtf ../GTF_merge_C1

```


f) Create gene matrix


```
#making a sample txt file with all gtf file names

F=/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/

array2=($(ls *merge.gtf))
for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_C1.txt
done

```
```
#sample_list.txt document

C17_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C17_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C18_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C18_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C19_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C19_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C20_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C20_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C21_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C21_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C22_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C22_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C23_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C23_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C24_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C24_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C25_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C25_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C26_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C26_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C27_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C27_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C28_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C28_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C29_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C29_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C30_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C30_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C31_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C31_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
C32_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/C32_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E10_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E10_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E11_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E11_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E12_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E12_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E13_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E13_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E14_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E14_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E15_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E15_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E16_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E16_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E1_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E1_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E2_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E2_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E3_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E3_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E4_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E4_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E5_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E5_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E6_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E6_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E7_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E7_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E8_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E8_R1_001.fastq.gz.sam.sorted.bam.merge.gtf
E9_R1_001.fastq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/E9_R1_001.fastq.gz.sam.sorted.bam.merge.gtf

```
```
#create gene matrix
```

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts_C1/GTFtoCounts_C1.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1
#SBATCH --cpus-per-task=3


module load StringTie/2.1.4-GCC-9.3.0
module load Python/2.7.18-GCCcore-9.3.0

python prepDE.py -g Sym_C1_gene_count_matrix.csv -i sample_list_C1.txt

```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts_C1/GTFtoCounts_C1.sh
Submitted batch job 1913031

```


g) Secure-copy gene counts onto local computer, make sure to open a seperate command shell outside of bluewaves on your own terminal

```
#copy Sym gene count matrix

scp danielle_becker@bluewaves.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped_C1/GTF_merge_C1/Sym_C1_gene_count_matrix.csv /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RAnalysis/Data/RNA-seq


#copy Sym transcript count matrix
scp danielle_becker@bluewaves.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/transcript_count_matrix.csv /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RAnalysis/Data/RNA-seq

```

