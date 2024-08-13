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

I am using the *Pocillopora verrucosa* genome since it is the closest complete genome for *Pocillopora meandrina and eydouxi*

[Buitrago-López et al 2020](https://academic.oup.com/gbe/article/12/10/1911/5898631)

Location on Andromeda, the HPC server for URI:
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

```
Core files	MD5 hash

Genome scaffolds    fb4d03ba2a9016fabb284d10e513f873
Gene models (CDS)   019359071e3ab319cd10e8f08715fc71
Gene models (proteins)    438f1d59b060144961d6a499de016f55
Gene models (GFF3)    614efffa87f6e8098b78490a5804c857

Miscellaneous files	MD5 hash
Full transcripts	            76b5d8d405798d5ca7f6f8cc4b740eb2

On Andromeda
Pver_genes_names_v1.0.fna.gz  019359071e3ab319cd10e8f08715fc71
```

### Functional annotation

[functional annotation xlsx file link](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/gbe/12/10/10.1093_gbe_evaa184/1/evaa184_supplementary_data.zip?Expires=1608305338&Signature=q7Y3SsLOTtlk5ZF6UMMLGy~HRtUFQqRbQt8ZiasnG3EeO11vmHcNgToGSowqYxQK1vibkmPEzMWDeS6u8qG~D20t7G31abz9zbpFrdW9T0cisHAQwY5g~lyK-WRFd-EDYW1eHFI4x~vU0G0xopva7kx1KlXdWxyZW86Fr7CDckFFvav78SAvZtmcvL8WuY4tWmEf33LK4ruuX7ZndqT8k~Kzag57phDdN1qleKWmeAf2wI-Wn8B4w-gV7UU4WQV1Ybs1wwdmexfPxH-DYEuSm-3T4sFq52FW1eRa8WD0V9XDUyysgajGh3sXRxHy-hEUUdnrzlVEk9~Doo9l9IIaUA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)


# make folder structure
```
mkdir data
cd data


ln -s ../../../../../KITT/hputnam/20201209_Becker_RNASeq_combo/combo/*.fastq.gz ./raw/
ln -s ../../../../../REFS/Pverr/ ./refs/

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
Submitted batch job 1816196 20201230
```

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

c) Verify data integrity with md5sum

Cross-reference the checksum document from GENEWIZ with the data we have on our computer

With a small amount of files, able to first cross-check that the sequences matched between both files on the desktop

Use the code below in terminal to cross-check the files and compare for sanity check
```
in directory: /data/putnamlab/KITT/hputnam/20201124_Becker_RNASeq

md5sum -c 20201124_Becker_RNASeq.md5
```
Should output 'OK' next to each file name

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

b) Write script for checking quality with FastQC and submit as job on Andromeda

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
Submitted batch job 1816766 on 20210104



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
scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/trimmed_qc/multiqc_report.html /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Bioinformatics/Output/RNASeq/trimmed.qc.multiqc

```


# 4) Trim and clean reads

a) Make trimmed reads folder in all other results folders

```

mkdir data/trimmed
cd trimmed

```

c) Write script for Trimming and run on Andromeda

#Run fastp on files
#Trims 20bp from 5' end of all reads
#Trims poly G, if present

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

array1=($(ls *R1*.fastq.gz)) #Make an array of sequences to trim
for i in ${array1[@]}; do
fastp --in1 ${i} --in2 $(echo ${i}|sed s/_R1/_R2/) --detect_adapter_for_pe --trim_poly_g --trim_front1 20 --trim_front2 20 --out1 ../trimmed/${i} --out2 ../trimmed/$(echo ${i}|sed s/_R1/_R2/)  
done

```
```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/trim.sh
Submitted batch job 235101 on 20220223 - started at 15:30 pm, ended at 01:30 am - 10 hours
```


# 5) Check quality of trimmed files

a) Check number of files in /trimmed directory

```
ls -1 | wc -l
#64
```


b) Check number of reads in /trimmed directory

```
zgrep -c "@GWNJ" *.gz > trimmed_seq_counts

C17_R1.fastp-trim.20230215.fq.gz:24075726
C17_R2.fastp-trim.20230215.fq.gz:24075726
C18_R1.fastp-trim.20230215.fq.gz:22164069
C18_R2.fastp-trim.20230215.fq.gz:22164069
C19_R1.fastp-trim.20230215.fq.gz:24108648
C19_R2.fastp-trim.20230215.fq.gz:24108648
C20_R1.fastp-trim.20230215.fq.gz:23308448
C20_R2.fastp-trim.20230215.fq.gz:23308448
C21_R1.fastp-trim.20230215.fq.gz:15900108
C21_R2.fastp-trim.20230215.fq.gz:15900108
C22_R1.fastp-trim.20230215.fq.gz:22276164
C22_R2.fastp-trim.20230215.fq.gz:22276164
C23_R1.fastp-trim.20230215.fq.gz:20264849
C23_R2.fastp-trim.20230215.fq.gz:20264849
C24_R1.fastp-trim.20230215.fq.gz:21739688
C24_R2.fastp-trim.20230215.fq.gz:21739688
C25_R1.fastp-trim.20230215.fq.gz:28400490
C25_R2.fastp-trim.20230215.fq.gz:28400490
C26_R1.fastp-trim.20230215.fq.gz:22356889
C26_R2.fastp-trim.20230215.fq.gz:22356889
C27_R1.fastp-trim.20230215.fq.gz:23741291
C27_R2.fastp-trim.20230215.fq.gz:23741291
C28_R1.fastp-trim.20230215.fq.gz:22662492
C28_R2.fastp-trim.20230215.fq.gz:22662492
C29_R1.fastp-trim.20230215.fq.gz:17264107
C29_R2.fastp-trim.20230215.fq.gz:17264107
C30_R1.fastp-trim.20230215.fq.gz:25439860
C30_R2.fastp-trim.20230215.fq.gz:25439860
C31_R1.fastp-trim.20230215.fq.gz:23852854
C31_R2.fastp-trim.20230215.fq.gz:23852854
C32_R1.fastp-trim.20230215.fq.gz:25904515
C32_R2.fastp-trim.20230215.fq.gz:25904515
E10_R1.fastp-trim.20230215.fq.gz:28094072
E10_R2.fastp-trim.20230215.fq.gz:28094072
E11_R1.fastp-trim.20230215.fq.gz:17954806
E11_R2.fastp-trim.20230215.fq.gz:17954806
E12_R1.fastp-trim.20230215.fq.gz:26996226
E12_R2.fastp-trim.20230215.fq.gz:26996226
E13_R1.fastp-trim.20230215.fq.gz:23732695
E13_R2.fastp-trim.20230215.fq.gz:23732695
E14_R1.fastp-trim.20230215.fq.gz:21831061
E14_R2.fastp-trim.20230215.fq.gz:21831061
E15_R1.fastp-trim.20230215.fq.gz:22614254
E15_R2.fastp-trim.20230215.fq.gz:22614254
E16_R1.fastp-trim.20230215.fq.gz:28776838
E16_R2.fastp-trim.20230215.fq.gz:28776838
E1_R1.fastp-trim.20230215.fq.gz:24465103
E1_R2.fastp-trim.20230215.fq.gz:24465103
E2_R1.fastp-trim.20230215.fq.gz:23120087
E2_R2.fastp-trim.20230215.fq.gz:23120087
E3_R1.fastp-trim.20230215.fq.gz:25025833
E3_R2.fastp-trim.20230215.fq.gz:25025833
E4_R1.fastp-trim.20230215.fq.gz:21164262
E4_R2.fastp-trim.20230215.fq.gz:21164262
E5_R1.fastp-trim.20230215.fq.gz:15836446
E5_R2.fastp-trim.20230215.fq.gz:15836446
E6_R1.fastp-trim.20230215.fq.gz:24285839
E6_R2.fastp-trim.20230215.fq.gz:24285839
E7_R1.fastp-trim.20230215.fq.gz:23378626
E7_R2.fastp-trim.20230215.fq.gz:23378626
E8_R1.fastp-trim.20230215.fq.gz:22972241
E8_R2.fastp-trim.20230215.fq.gz:22972241
E9_R1.fastp-trim.20230215.fq.gz:24270035
E9_R2.fastp-trim.20230215.fq.gz:24270035


scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/trimmed_seq_counts /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Bioinformatics/Output/RNASeq/

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
module load MultiQC/1.9-intel-2020a-Python-3.8.2
multiqc /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/trimmed_qc
```
```
scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/trimmed_qc/*.html /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Bioinformatics/Output/RNASeq/trimmed.qc.multiqc

```

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
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/Hisat2_genome_build.sh
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
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs
#SBATCH --error="script_error"
#SBATCH --output="output_script"

module load HISAT2/2.1.0-foss-2018b

hisat2-build -f /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0.fasta ./Pver_ref

```
```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/Hisat2_genome_build.sh
Submitted batch job 1834918
```

b) Align reads to genome

```
mkdir mapped
```
```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/Hisat2_align2.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3
#SBATCH --error="script_error"
#SBATCH --output="output_script"


module load HISAT2/2.1.0-foss-2018b

#Aligning paired end reads
#Has the R1 in array1 because the sed in the for loop changes it to an R2. SAM files are of both forward and reverse reads

array1=($(ls *R1*.fq.gz))
for i in ${array1[@]}; do
hisat2 -p 48 --rna-strandness RF --dta -q -x /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pver_ref -1 /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/${i} \
-2 /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/$(echo ${i}|sed s/_R1/_R2/) -S /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/${i}.sam
done
```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/Hisat2_align2.sh
#Submitted batch job 235530 20230223 at 17:14, ended at 21:18, 3 hours 44 minutes

#Download alignment statistics information from mapping, will be the output from your script above, even though it is a mapped output, you are in the trimmed folder here

scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/trimmed/script_error /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Bioinformatics/Output/RNASeq/

```

## Sort and convert sam to bam

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/SAMtoBAM.sh

#There will be lots of .tmp file versions in your folder, this is normal while this script runs and they should delete at the end to make one sorted.bam file
```
```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped
#SBATCH --cpus-per-task=3
#SBATCH --error="script_error"
#SBATCH --output="output_script"

module load SAMtools/1.9-foss-2018b

array1=($(ls *.sam))  
for i in ${array1[@]}; do
samtools sort -o /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/${i}.sorted.bam /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/${i}  
done
```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/SAMtoBAM.sh
#Submitted batch job 235853 on 20230224 at 09:02, finished at 06:51, was 18 hours total
```

### Remove Sam files to save space
```
rm /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/*.sam

```


### Check number of mapped reads
```
Explanation for SAMtools functions for checking the mapped reads in a paired-end dataset found here: https://www.biostars.org/p/138116/

nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/mapped_read_counts.sh

#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped
#SBATCH --cpus-per-task=3

module load SAMtools/1.9-foss-2018b

array1=($(ls *.bam))  
for i in ${array1[@]}; do
samtools view -F 0x4 /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/${i} | cut -f 1 | sort | uniq | wc -l > mapped_read_counts
done

sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/mapped_read_counts.sh

#Submitted batch job 236602, took 7 hours
#19421210 counts

```

# 7) Perform gene counts with stringTie

### Needed to modify Pverr_genome_assembly file, following steps described below

When I was trying to create my gene matrix using the prepDE.py script below, I was receiving an error relating to the transcript ID's in my .gtf files.

Error:
Traceback (most recent call last):
File "prepDE.py", line 158, in
t_id=RE_TRANSCRIPT_ID.search(v[8]).group(1)
AttributeError: 'NoneType' object has no attribute 'group'

Once the error was identified to be in my .gtf files, we saw that there was a tab in between the 1 and 3 and that true was within the parentheses below.

Example corruption files:

C17_R1_001.fastq.gz.sam.sorted.bam.merge.gtf, line 5925:
Pver_Sc0000004_size1560107 StringTie transcript 1538 1660 1000 - . gene_id "STRG.507"; transcript_id "Pver_g535.t1 3_prime_partial true"; cov "0.000013"; FPKM "0.000004"; TPM "0.000011";

Pver_xfSc0007555_size1079 . CDS 3 1079 . + 0 ID=cds.Pver_g25526.t1;Parent=Pver_g25526.t1 5_prime_partial true 3_prime_partial true

To check the frequency and occurrence of these errors within my files, we used:

grep transcript.*transcript_id.*true *gtf

It showed that all of my files contained these errors. The 3_prime_partial and 5_prime_partial labels that were being created when using StringTie/2.1.4-GCC-9.3.0 and our reference Pver_genome_assembly_v1.0.gff3 file to the assemble and estimate reads was the main issue.

Once identifying this, we used the command:

sed -r -e 's/(Parent=[^[:space:]])./\1/' Pver_genome_assembly_v1.0.gff3 > fixed_Pver_genome_assembly_v1.0.gff3

to edit the original Pver_genome_assembly_v1.0.gff3 file in my directory to remove the 5_prime_partial true, 3_prime_partial true, and the space causing the file corruption issue. This command keeps the text (non-space characters) after Parent= and discards the rest of the line.

After making this correction to the Pver_genome_assembly_v1.0.gff3 file, I tested the new fixed_Pver_genome_assembly_v1.0.gff3 reference file on one sample file and was able to work through the entire workflow with no issues to create a gene matrix.

Code authored by Polina Shpilker to amend issue:

```
input_file = "Pver_genome_assembly_v1.0.gff3"
output_file = "Pver_genome_assembly_v1.0_modified.gff3"
log_file = "annotation_repair.log"

## Notes:

.gff3 file is separated by gene blocks - i.e. until you reach a new gene, it is safe to assume that all current annotations are related to the last gene you saw.

I don't want to assume that there can be only 1 mrna per gene. There is also a similar assumption as to the gene
assumption; until you reach a new mRNA/gene annotation, all annotations are related to the last mRNA you saw.

The problem we are trying to solve is that for some reason these non-mRNA and non-gene annotations have the wrong parent ID. They are, in some cases, listing their own ID as the parent ID.

## Method to the madness:

- Scan the file for a 'mRNA' tag
- Until you reach another mRNA tag:
- read in every annotation
- check that its parent ID is set to the mRNA ID you last saw
- If it is: spit into the next file without comment.
- If it isn't: spit into the next file with the correct parent ID, take note.

## Unless otherwise stated, every line is copied into the output file exactly as-is.

## Helper function.

Takes in the string of the entirety of the final column, splits it into the individual attributes.
Then extracts the ID and Parent attributes, and returns them as a tuple.

def get_id_pid(datcol) :
     sub-columns within the final column are separated by a semicolon
    cols = datcol.split(';')

Get only the attributes 'ID' and 'Parent'
    i_id = None
    i_parent = None
    for col in cols :
        - Every attribute is stored as 'name=value', so by splitting on '=' we can get back the name of the attribute
        -   and the value of the attribute.
        (attr, val) = col.split("=")
        if attr == "ID" :
            i_id = val
        elif attr == "Parent" :
            i_parent = val

    - If we could not find values for this feature's ID and its parent's ID, throw an exception. Otherwise, return
    -   the values we found.
    if i_id == None or i_parent == None :
        raise ValueError("Could not find both ID and parent ID in data: " + datcol)
    else :
        return(i_id, i_parent)

## Helper function.

- Takens in the string of the entirety of the annotation column, and what to change the parent ID to.
- Creates an array that stores all of the attr=value strings of all attributes in the column, except for the parent
-   attribute. Once it reaches the parent attribute, it replaces it with an attr=value string containing the new_parent ID for the value.
- Then returns a string that contains all of these attr=value strings combined with semicolons separating them.
def adjust_pid(datcol, new_parent) :
    cols = datcol.split(';')
    output = []
    for col in cols :
        (attr, val) = col.split("=")
        if(attr == "Parent") :
            output.append("Parent=" + new_parent)
        else :
            output.append(col)
    return(";".join(output) + "\n")

##Helper function.

- Takes in an array of data; the first index is the correctly formatted data column that follows the gff standard.
-   Every following index contains a space-separated pair of attribute and its value that needs to be added to the correctly-formatted column.

- Example:
- INPUT     ["ID=cds.Pver_g25526.t1;Parent=Pver_g25526.t1", "5_prime_partial true", "3_prime_partial true"]
- OUTPUT    "ID=cds.Pver_g25526.t1;Parent=Pver_g25526.t1;5_prime_partial=true;3_prime_partial=true"
def fix_extra_columns(datcols) :
    - Get the correctly formatted first index
    output = datcols[0]
    - For every following index, split by space and add to the output string
    for col in datcols[1:] :
        split = col.split(" ")
        output = output + ";" + split[0] + "=" + split[1]
    return output

cur_mrna_id = "NA"
cur_id = None
cur_parent = None

with open(input_file, 'r') as inp :
    with open(output_file, 'w') as out :
        with open(log_file, 'w') as log:
            for lineno, line in enumerate(inp, 1):
                - If it starts with a '#', ignore it and spit it into the output file. These are just comments.
                if line[0] == "#" :
                    out.write(line)
                    continue

                columns = line.strip().split("\t")

                dat_type = columns[2]

                - If there are more columns than expected, that means we have the additional 5_prime_partial or
                -   3_prime_partial tags incorrectly appended. Fix the attribute column to store this data, and then
                -  get rid of the extra columns. Log that the change has been made.
                if(len(columns) > 9) :
                    columns[8] = fix_extra_columns(columns[8:])
                    columns = columns[:9]
                    line = "\t".join(columns) + "\n"
                    log.write("Modified line: " + str(lineno) + "\tFixed misplaced attributes.\n")

                attributes = columns[8]


                - if the 3rd column is not gene or mRNA -- and we have not yet seen an mrna -- throw an error.
                -   This is only because of the assumptions I made up in the notes section; I want to know if they're not
                -   valid! If I can't assume the block-like structure where every exon, CDS, etc is related to the mRNA
                -  annotation proceeding it, we flat out can't repair the file because we have no idea what the real
                -   mRNA parent is.
                if cur_mrna_id == "NA" and dat_type != "gene" and dat_type != "mRNA":
                    raise AssertionError("Assumption of file structure is incorrect. Repair cannot proceed.")

                - if the type is gene, just paste it into the next file without anything fancy.
                if dat_type == "gene" :
                    out.write(line)
                - If the type name is mRNA, update our stored mRNA name and paste it into the output file without any changes.
                elif dat_type == "mRNA":
                    cur_mrna_id, cur_parent = get_id_pid(attributes)
                    out.write(line)
                - If it is of any other type, check to make sure that the parent ID is set correctly.
                - If not, log it in our log file & write to the output file an adjusted line.
                else :
                    cur_id, cur_parent = get_id_pid(attributes)
                    if(cur_parent != cur_mrna_id) : #the parent is supposed to be the id of the mrna we saw last
                        log.write("Modified line: " + str(lineno) + "\tTo: " + cur_mrna_id + "\tFrom: " + cur_parent + "\n")
                        tmp = columns[:8]
                        tmp.append(adjust_pid(columns[8], cur_mrna_id))

                        out.write("\t".join(tmp))
                    else :
                        out.write(line)

```

```
##copy modified genome assembly file to Andromeda, enter this command into local computer shell

scp -r /Users/Danielle/Downloads/Pver_genome_assembly_v1.0_modified.gff3 danielle_becker@Andromeda.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/
```

```
mkdir counts
cd counts

```


b) Assemble and estimate reads

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/StringTie_Assemble.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped
#SBATCH --cpus-per-task=3

module load StringTie/2.2.1-GCC-11.2.0

array1=($(ls *.bam))
for i in ${array1[@]}; do
stringtie -p 48 --rf -e -G /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0_modified.gff3 -o /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/${i}.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/${i}
done
```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/StringTie_Assemble.sh
Submitted batch job 237016, 1 hour 20 minutes
```

c) Merge stringTie gtf results

#in this step we are making a file with all the gtf names and stringtie will merge them all together for a master list for your specific genes

```
ls *gtf > mergelist.txt
cat mergelist.txt

module load StringTie/2.2.1-GCC-11.2.0

stringtie --merge -p 8 -G /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0_modified.gff3 -o stringtie_merged.gtf mergelist.txt

```

d) Assess assembly quality

```

nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/gffcompare.sh

#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts
#SBATCH --cpus-per-task=3

module load GffCompare/0.12.6-GCC-11.2.0

gffcompare -r /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0_modified.gff3 -o merged stringtie_merged.gtf

sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/gffcompare.sh

#Submitted batch job 237924, took one minute

```

```

#= Summary for dataset: stringtie_merged.gtf
#     Query mRNAs :   27439 in   27401 loci  (24124 multi-exon transcripts)
#            (38 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   27439 in   27401 loci  (24124 multi-exon)
# Super-loci w/ reference transcripts:    27401
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |   100.0    |
        Exon level:   100.0     |   100.0    |
      Intron level:   100.0     |   100.0    |
Intron chain level:   100.0     |   100.0    |
  Transcript level:   100.0     |   100.0    |
       Locus level:   100.0     |   100.0    |

     Matching intron chains:   24124
       Matching transcripts:   27438
              Matching loci:   27400

          Missed exons:       0/208636  (  0.0%)
           Novel exons:       0/208633  (  0.0%)
        Missed introns:       0/181193  (  0.0%)
         Novel introns:       0/181193  (  0.0%)
           Missed loci:       0/27401   (  0.0%)
            Novel loci:       0/27401   (  0.0%)

 Total union super-loci across all input datasets: 27401
27439 out of 27439 consensus transcripts written in merged.annotated.gtf (0 discarded as redundant)

For matching transcripts we noticed that one transcript did not match the reference, so we identified the transcript as TCONS_00002577: "Pver_gene_novel_gene_109_5de57afd" and it had a class code of '-k' listed which indicates that this transcript is fully contained/enclosed within a reference transcript from the provided reference annotation

'=' Flag:
The '=' flag indicates that the query transcript has an identical intron chain as a reference transcript from the provided annotation. This means the query transcript has exactly the same exon-intron structure and coordinates as a reference transcript, except potentially different lengths for the terminal exons.

'k' Flag:
The 'k' flag means the query transcript is fully contained within (enclosed by) a reference transcript from the annotation. 'k' indicates it is a complete subset of a reference but with potential differences in terminal exon lengths or other minor variations from the reference model.

```

e) Re-estimate assembly (this step is redundant and not necessary unless you have an updated reference or another assembly you need to cross-reference)

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/re_estimate.assembly.sh
```
```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped
#SBATCH --cpus-per-task=3

module load StringTie/2.2.1-GCC-11.2.0

array1=($(ls *.bam))
for i in ${array1[@]}; do
stringtie -e -G /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0_modified.gff3 -o ${i}.merge.gtf ${i}
echo "${i}"
done

```
```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/re_estimate.assembly.sh
Submitted batch job 237925, 1 hour 24 minutes
```

```
# move merged GTF files to their own folder
mkdir GTF_merge

mv *merge.gtf GTF_merge

```

f) Create gene matrix


```
#making a sample txt file with all gtf file names

F=/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/

array2=($(ls *merge.gtf))
for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list.txt
done

```
```
#sample_list.txt document

C17_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C17_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C18_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C18_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C19_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C19_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C20_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C20_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C21_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C21_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C22_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C22_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C23_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C23_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C24_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C24_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C25_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C25_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C26_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C26_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C27_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C27_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C28_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C28_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C29_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C29_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C30_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C30_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C31_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C31_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
C32_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/C32_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E10_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E10_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E11_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E11_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E12_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E12_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E13_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E13_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E14_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E14_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E15_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E15_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E16_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E16_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E1_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E1_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E2_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E2_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E3_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E3_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E4_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E4_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E5_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E5_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E6_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E6_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E7_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E7_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E8_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E8_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf
E9_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/E9_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam.merge.gtf

```
```
#create gene matrix
```

```
nano /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/GTFtoCounts.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge
#SBATCH --cpus-per-task=3


module load StringTie/2.2.1-GCC-11.2.0
module load Python/2.7.18-GCCcore-9.3.0

python prepDE.py -g Poc_gene_count_matrix.csv -i sample_list.txt

```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/GTFtoCounts.sh
Submitted batch job 237934
```


g) Secure-copy gene counts onto local computer, make sure to open a seperate command shell outside of Andromeda on your own terminal

```
#copy gene count matrix

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/Poc_gene_count_matrix.csv /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RAnalysis/Data/RNA-seq/Host


#copy transcript count matrix
scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/GTF_merge/transcript_count_matrix.csv /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RAnalysis/Data/RNA-seq/Host

```


# 8) Compare modified.gff3 and fixed.gff3

To ensure that our modified.gff3 did not change any necessary transcripts or information needed for downstream analysis, we are going to compare the counts assigned to each gene from one sample to a fixed.gff3 that was composed using the below Rscript:


This script adds transcript and gene id into GFF file for alignment.  

Here, I'll be adding transcript_id= and gene_id= to 'gene' column because we needs that label to map our RNAseq data  

Load libraries and data.
```{r}
#Load libraries
library(tidyverse)
library(R.utils)
```

Load  gene gff file

```{r}
gff <- read.csv(file="../../Pver_genome_assembly_v1.0.gff3", header=FALSE, sep="\t", skip=1) #gff <- na.omit(gff)

#gff.DB <- read.csv(file="data/Tagseq_zymo_genohub/Pver_genome_assembly_v1.0_modified.gff3", header=FALSE, sep="\t", skip=1)
gff.DB <- na.omit(gff.DB)
```

Rename columns ``{r}
colnames(gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(gff)

unique(gff$scaffold)
#shows that #PROT and # PASA_UPDATE # ORIGINAL: lines have been included but only have information in column 1 scaffold
```

Remove extra #PROT and # PASA_UPDATE and # ORIGINAL:
check the number of genes against the number in the genome
```{r}
#genes in the Pverrucosa genome 27,439

# remove any rows with NA
gff <- na.omit(gff)

# check the number of genes against the number in the genome
# Pverrucosa should return 27439
nrow(gff %>%filter(id=="gene"))

#genes <- gff %>%filter(id=="mRNA")
#genes$gene <- gsub("Pver_gene_split_gene_g", "Pver_split_gene_g", genes$gene) #remove extra "_gene_" from name

```

Remove cds. and duplication of _gene in the names
```{r}
#gff$ID <- sub(";.*", "", gff$gene) #move everything before the first ";" to a new column
#gff$ID <- gsub("cds.", "", gff$gene) #remove "cds." from the new column
gff$gene <- gsub("Pver_gene_g", "Pver_g", gff$gene) #remove "_gene" from the new column
gff$gene <- gsub("Pver_gene_split_gene_g", "Pver_split_gene_g", gff$gene) #remove extra "_gene_" from name
gff$gene <- gsub("Pver_gene_novel_model_", "Pver_novel_model_", gff$gene) #remove extra "_gene_" from name

head(gff)
```

Create a Parent column.
`{r}
gff$parent_id <- sub(".*Parent=", "", gff$gene)
gff$parent_id <- sub(";.*", "", gff$parent_id)
gff$parent_id <- gsub("ID=", "", gff$parent_id) #remove ID= ha(gff)
```

Create a transcript ID column
```{r}
gff$transcript_id <- sub(".*Name=", "", gff$gene)
gff$transcript_id <- gsub("ID=", "", gff$transcript_id) # remove ID
gff$transcript_id <- gsub("Parent=", "", gff$transcript_id) # remove Parent
gff$transcript_id <- sub(".*;", "", gff$transcript_id) #keep everything after semi colon
head(gff)
```

Now add these values into the gene column separated by semicolons.  

```{r}
gff <- gff %>%
 uate(test = ifelse(id != "gene", paste0("ID=", gff$transcript_id, ";parent_id=", gff$parent_id),  paste0(gene)))
head(gff)
```

Check which things are not matching between the gff and annotation file. ``r
#Read in annotation file to check gene names are the same in the gff and the "Query" column of the annotation
Annot <- readxl::read_xlsx("../../../Downloads/evaa184_supplementary_data/FileS2_Pver_gene_annot_May28.xlsx", skip=4)

gff_list<-gff$transcript_id
head(gff_list)

annot_list<-Annot$Query
head(annot_list)

missing1 <- setdiff(gff_list, annot_list) #things in gff that are not in annotation
missing1
length(missing1)
#the gff has 2543 extra things
msig2 <- setdiff(annot_list, gff_list) #things in annotation that are not in gff
missing2
length(missing2)

```

Now remove the transcript and parent ID separate columns.  

```{r}
gff<-gff %>%
  dplyr::select(!transcript_id)%>%
  dplyr::select(!parent_id)%>%
  dplyr::select(!gene)

head(gff)

gff<-gff%>%
  dplyr::rename(gene=test)

head(gff)
```

Save file. Gzip and then upload this to the server for use in bioinformatic steps.  
```{r}
write.table(gff, file="../../Pver_genome_assembly_v1.0_fixed.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

#gzip the file gi(./../Pver_genome_assembly_v1.0_fixed.gff3")
```

# Now we have a fixed.gff3 file edited from the original published reef genomics .gff3 that should match our modified.gff3 for all downstream analyses, we will also test the original published reef genomics .gff3

# We will be using sample E9's .bam file to compare between the three gff3 files

a) Assemble and estimate reads for each of the different .gff3s

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/StringTie_Assemble_gff3_compare.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped
#SBATCH --cpus-per-task=3

module load StringTie/2.2.1-GCC-11.2.0

# Run the first stringtie command in the background
stringtie -p 48 --rf -e -G /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0_fixed.gff3 -o /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/gff3.count.compare/E9_R1_fixed.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/E9_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam &

# Run the second stringtie command in the background
stringtie -p 48 --rf -e -G /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0_modified.gff3 -o /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/gff3.count.compare/E9_R1_modified.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/E9_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam &

# Run the third stringtie command in the background
stringtie -p 48 --rf -e -G /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0.gff3 -o /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/gff3.count.compare/E9_R1_original.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/E9_R1.fastp-trim.20230215.fq.gz.sam.sorted.bam &

# Wait for both commands to finish
waitsqueu

```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/StringTie_Assemble_gff3_compare.sh

Submitted batch job 322266
```




f) Create gene matrix


```
#making a sample txt file with all gtf file names

F=/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/gff3.count.compare/

array2=($(ls *.gtf))
for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list.txt
done

```
```
#sample_list.txt document

E9_R1_fixed.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/gff3.count.compare/E9_R1_fixed.gtf
E9_R1_modified.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/gff3.count.compare/E9_R1_modified.gtf
E9_R1_original.gtf /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/gff3.count.compare/E9_R1_original.gtf

```



```
# make gene matrix script using updated prepDY.py pyton script

wget https://raw.githubusercontent.com/gpertea/stringtie/master/prepDE.py3


nano GTFtoCounts_gff3_compare.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/gff3.count.compare
#SBATCH --cpus-per-task=3

module load StringTie/2.2.1-GCC-11.2.0
module load Python/3.9.6-GCCcore-11.2.0

# Output directory for count matrices
OUTPUT_DIR="/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/gff3.count.compare"

# Path to prepDE.py script
PREPDE_SCRIPT="/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/counts/gff3.count.compare/prepDE.py3"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through each line in sample_list.txt
while IFS= read -r line; do
    # Extract the GTF filename and its path from the line
    filename=$(echo "$line" | cut -d ' ' -f 1)
    filepath=$(echo "$line" | cut -d ' ' -f 2)

    # Extract the sample name (without extension) from the filename
    sample_name=$(basename "$filename" .gtf)

    # Create a temporary sample list for the current sample
    echo "$filename $filepath" > temp_sample_list.txt

    # Run prepDE.py for each GTF file
    python "$PREPDE_SCRIPT" -i temp_sample_list.txt -g "${OUTPUT_DIR}/${sample_name}_gene_count_matrix.csv" -t "${OUTPUT_DIR}/${sample_name}_transcript_count_matrix.csv"

    # Remove the temporary sample list
    rm temp_sample_list.txt

done < sample_list.txt

echo "Count matrices generated successfully."


```

```
sbatch GTFtoCounts_gff3_compare.sh

Submitted batch job 322279

```

# 9) Run RNAseq SNP calling and genetic relatedness

In our data we see some samples (n=5) that are in different blocks but the same enriched treatment clumping in our heatmap. These are the same haplotype from our mtORF and PocHIstone molecular analyses so we want to check the genetic relateness. We must first use SNP calling from the RNAseq data, then calculate genetic relatedness of our samples.

First, create a script for SNP calling with bcftools and samtools functions:

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/SNP_calling.sh

```

```
#!/bin/bash
#SBATCH --job-name=snp_calling
#SBATCH --output=snp_calling.out
#SBATCH --error=snp_calling.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=96:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH --mem=120G

# Load required modules
module load BCFtools/1.17-GCC-12.2.0
module load SAMtools/1.9-foss-2018b

# Set input and output directories
INPUT_DIR="/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped"
OUTPUT_DIR="/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/SNP_calling"

# Reference genome
REF_GENOME="/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0.fasta"

# Call SNPs
bcftools mpileup -Ou -f $REF_GENOME $INPUT_DIR/*.bam | \
bcftools call -mv -Oz -o $OUTPUT_DIR/variants.vcf.gz

# Index the VCF file
bcftools index $OUTPUT_DIR/variants.vcf.gz

```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/SNP_calling.sh

Submitted batch job 334327
```

Second, filter SNPs:

```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/SNP_filtering.sh

```

```
#!/bin/bash
#SBATCH --job-name=snp_filtering
#SBATCH --output=snp_filtering.out
#SBATCH --error=snp_filtering.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=96:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH --mem=120G

# Load required modules
module load BCFtools/1.17-GCC-12.2.0

# Set input and output files
INPUT_VCF="/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/SNP_calling/variants.vcf.gz"
OUTPUT_VCF="/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/SNP_calling/filtered_variants.vcf.gz"

# Filter SNPs
bcftools filter -i 'QUAL>20 && DP>10' $INPUT_VCF -Oz -o $OUTPUT_VCF

# Index the filtered VCF file
bcftools index $OUTPUT_VCF

```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/SNP_filtering.sh

```

Third, calculate genetic relatedness using PLINK :


```
nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/genetic_relatedness.sh

```

```
#!/bin/bash
#SBATCH --job-name=genetic_relatedness
#SBATCH --output=genetic_relatedness.out
#SBATCH --error=genetic_relatedness.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=96:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH --mem=120G

# Load required modules
module load PLINK/2.00a3.7-gfbf-2023a

# Set input and output files
INPUT_VCF="/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/SNP_calling/filtered_variants.vcf.gz"
OUTPUT_PREFIX="/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/mapped/SNP_calling/relatedness"

# Calculate genetic relatedness
plink --vcf $INPUT_VCF --double-id --allow-extra-chr --make-grm-gz --out $OUTPUT_PREFIX

```

```
sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/scripts/genetic_relatedness.sh

Submitted batch job 334351
```

Download relatedness.grm.gz to Desktop to visualize results of relatedness using this [Rscript](https://github.com/hputnam/Becker_E5/blob/master/Bioinformatics/Workflows/genetic_relatedness.Rmd).
