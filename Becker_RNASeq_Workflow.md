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

symlink files in RNASeq_Becker_E5

```
mkdir data
cd data
mkdir raw
mkdir refs

ln -s ../../../../KITT/hputnam/20201124_Becker_RNASeq/*.fastq.gz ./raw/
ln -s ../../../../REFS/Pverr/ ./ref/
```

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

md5sum /data/putnamlab/KITT/hputnam/20201124_Becker_RNASeq/*.gz > URIcheckmd5.md5

```

```
sbatch ../../scripts/check_transfer.sh 
```

### Checksum from Genewiz
```

34d81ff7934fb0b46e0eae693b48e917  C17_R1_001.fastq.gz
e4b8193cc41a76fc6804d158ca391145  C17_R2_001.fastq.gz
9642a7849e9687b697a06624641bedd2  C18_R1_001.fastq.gz
ff551943373f1e66102ccb5cbac94977  C18_R2_001.fastq.gz
ca0e0391351b0c41eed5edf15c8661ea  C19_R1_001.fastq.gz
6e1360fe186f2700747fcd2113b681b9  C19_R2_001.fastq.gz
22f7d0701bf574e95f8b84f5f8a08f8b  C20_R1_001.fastq.gz
8f63713efe5128c7101a3cdad365d684  C20_R2_001.fastq.gz
0c71965c534a2ebf1b38cb97b7f879a0  C21_R1_001.fastq.gz
b7df9377b56654d0e8d8c1ab42b60853  C21_R2_001.fastq.gz
cb1bf0af9c614c4ecdf9936ab038b89c  C22_R1_001.fastq.gz
6fe6527212de3b8306d13f8222302f57  C22_R2_001.fastq.gz
91ead175cd247cd5b8dc530c6abcc9e5  C23_R1_001.fastq.gz
ec668bb851766ef36d3c7ba7943bcf8e  C23_R2_001.fastq.gz
4d117bcfc72763a7db319a8fd13ed908  C24_R1_001.fastq.gz
b86a886e1dd61ea6f385898623e7af9e  C24_R2_001.fastq.gz
932b47049b4e8e90280f299428a597a4  C25_R1_001.fastq.gz
d2c313c69347dbb255c6106abc4d22f8  C25_R2_001.fastq.gz
98552602cac6bdbed1bb2008991f727b  C26_R1_001.fastq.gz
56600a80d327a943c1ecda670bc71ed3  C26_R2_001.fastq.gz
d540c45f1148373b99c1db4e43f60850  C27_R1_001.fastq.gz
cb8dda30d5900a8d42a90c7d1dda9744  C27_R2_001.fastq.gz
7ad2c8e05101f77a822715b79783db72  C28_R1_001.fastq.gz
15932d696eb051ba8af9f95381b2c0ee  C28_R2_001.fastq.gz
67da8fa379000428ee2b726d5b6d2c3f  C29_R1_001.fastq.gz
fc32dc17927c41680c7153b6578b0705  C29_R2_001.fastq.gz
e1971e0cb3bf6d888feac28d53cd8706  C30_R1_001.fastq.gz
1ebf88060ca850b8e77bf664ea97b10e  C30_R2_001.fastq.gz
8da01ecfc3147c2b64c9dfba6aa5738a  C31_R1_001.fastq.gz
99d44c69480270acdbd12c432a51ec33  C31_R2_001.fastq.gz
9b21a38a31156fdee7dc1958d0600bbe  C32_R1_001.fastq.gz
41d7fa676d75b0bb57e8fc7a01e3ec13  C32_R2_001.fastq.gz
387531bdfa3d8f7930a1f0c25548e429  E10_R1_001.fastq.gz
b8484541ab8dae01cc75e5d289a49c61  E10_R2_001.fastq.gz
1c428b3f05ff5c1e02a8c5769bf29ccd  E11_R1_001.fastq.gz
3e6c5eb4320363328a9be320cf885873  E11_R2_001.fastq.gz
956104a0993692ff076a054b19a24961  E12_R1_001.fastq.gz
edc41bfcc4daa59e89dfabb909ca4e4d  E12_R2_001.fastq.gz
43d09374f4005122d7639b0b6853fb52  E13_R1_001.fastq.gz
ecb847106b11106cc82225c8901f3377  E13_R2_001.fastq.gz
92c4bf857a315cf84e00e2b2c749a188  E14_R1_001.fastq.gz
1901fbdf138e13782842899a01692139  E14_R2_001.fastq.gz
f5758d5b2ca705bb6367a254da0f39a0  E15_R1_001.fastq.gz
087911a4cdcb1cf776a38c2a590cf8d1  E15_R2_001.fastq.gz
6c5eb3478cbedee31e0d3263b685ffd7  E16_R1_001.fastq.gz
1fea0b437baa7c034cf188045678eec9  E16_R2_001.fastq.gz
1da4e491aaab9f8dd31ec1cceb9df6e4  E1_R1_001.fastq.gz
f648c57469094634c5db06399de71adf  E1_R2_001.fastq.gz
bdb29d094b4e99accb4901129f5632cc  E2_R1_001.fastq.gz
7a901ba3ffb15757fa1f35d6fdd77977  E2_R2_001.fastq.gz
57003d309200f5e00bc1ae0ef7efc569  E3_R1_001.fastq.gz
945adabf4bb06a71750c8a7dcb2ffad4  E3_R2_001.fastq.gz
9c1a16d8efaf9e8b393c57e59cf582ed  E4_R1_001.fastq.gz
3ec5d310d5ca71be11fb61da8f89e16b  E4_R2_001.fastq.gz
bcbd03903528b4d23694aeac58cba29b  E5_R1_001.fastq.gz
729188fdac7ca7d5b984cd92c7cacf80  E5_R2_001.fastq.gz
03188b855cddf9c4310f2068b43461b6  E6_R1_001.fastq.gz
0a4f9de84fe32baea6f7cab4e3fa67a4  E6_R2_001.fastq.gz
a8f4c26ff35730b89ea72913a56721a7  E7_R1_001.fastq.gz
91502a50cf91652e37f305edd986774b  E7_R2_001.fastq.gz
8a6e9fc499a8fca32614cd0e0e86fb1d  E8_R1_001.fastq.gz
cd7de23f56172a399144a2dc6ae7911b  E8_R2_001.fastq.gz
e0e620f5822277987ec718839df370de  E9_R1_001.fastq.gz
871c3316046af6bbdbef3cf80ce90d91  E9_R2_001.fastq.gz
```

### Checksum from files on Bluewaves

```
34d81ff7934fb0b46e0eae693b48e917  C17_R1_001.fastq.gz
e4b8193cc41a76fc6804d158ca391145  C17_R2_001.fastq.gz  
9642a7849e9687b697a06624641bedd2  C18_R1_001.fastq.gz  
ff551943373f1e66102ccb5cbac94977  C18_R2_001.fastq.gz
ca0e0391351b0c41eed5edf15c8661ea  C19_R1_001.fastq.gz
6e1360fe186f2700747fcd2113b681b9  C19_R2_001.fastq.gz
22f7d0701bf574e95f8b84f5f8a08f8b  C20_R1_001.fastq.gz
8f63713efe5128c7101a3cdad365d684  C20_R2_001.fastq.gz
0c71965c534a2ebf1b38cb97b7f879a0  C21_R1_001.fastq.gz
b7df9377b56654d0e8d8c1ab42b60853  C21_R2_001.fastq.gz
cb1bf0af9c614c4ecdf9936ab038b89c  C22_R1_001.fastq.gz
6fe6527212de3b8306d13f8222302f57  C22_R2_001.fastq.gz
91ead175cd247cd5b8dc530c6abcc9e5  C23_R1_001.fastq.gz
ec668bb851766ef36d3c7ba7943bcf8e  C23_R2_001.fastq.gz
4d117bcfc72763a7db319a8fd13ed908  C24_R1_001.fastq.gz
b86a886e1dd61ea6f385898623e7af9e  C24_R2_001.fastq.gz
932b47049b4e8e90280f299428a597a4  C25_R1_001.fastq.gz
d2c313c69347dbb255c6106abc4d22f8  C25_R2_001.fastq.gz
98552602cac6bdbed1bb2008991f727b  C26_R1_001.fastq.gz
56600a80d327a943c1ecda670bc71ed3  C26_R2_001.fastq.gz
d540c45f1148373b99c1db4e43f60850  C27_R1_001.fastq.gz
cb8dda30d5900a8d42a90c7d1dda9744  C27_R2_001.fastq.gz
7ad2c8e05101f77a822715b79783db72  C28_R1_001.fastq.gz
15932d696eb051ba8af9f95381b2c0ee  C28_R2_001.fastq.gz
67da8fa379000428ee2b726d5b6d2c3f  C29_R1_001.fastq.gz
fc32dc17927c41680c7153b6578b0705  C29_R2_001.fastq.gz
e1971e0cb3bf6d888feac28d53cd8706  C30_R1_001.fastq.gz
1ebf88060ca850b8e77bf664ea97b10e  C30_R2_001.fastq.gz
8da01ecfc3147c2b64c9dfba6aa5738a  C31_R1_001.fastq.gz
99d44c69480270acdbd12c432a51ec33  C31_R2_001.fastq.gz
9b21a38a31156fdee7dc1958d0600bbe  C32_R1_001.fastq.gz
41d7fa676d75b0bb57e8fc7a01e3ec13  C32_R2_001.fastq.gz
387531bdfa3d8f7930a1f0c25548e429  E10_R1_001.fastq.gz
b8484541ab8dae01cc75e5d289a49c61  E10_R2_001.fastq.gz
1c428b3f05ff5c1e02a8c5769bf29ccd  E11_R1_001.fastq.gz
3e6c5eb4320363328a9be320cf885873  E11_R2_001.fastq.gz
956104a0993692ff076a054b19a24961  E12_R1_001.fastq.gz
edc41bfcc4daa59e89dfabb909ca4e4d  E12_R2_001.fastq.gz
43d09374f4005122d7639b0b6853fb52  E13_R1_001.fastq.gz
ecb847106b11106cc82225c8901f3377  E13_R2_001.fastq.gz
92c4bf857a315cf84e00e2b2c749a188  E14_R1_001.fastq.gz
1901fbdf138e13782842899a01692139  E14_R2_001.fastq.gz
f5758d5b2ca705bb6367a254da0f39a0  E15_R1_001.fastq.gz
087911a4cdcb1cf776a38c2a590cf8d1  E15_R2_001.fastq.gz
6c5eb3478cbedee31e0d3263b685ffd7  E16_R1_001.fastq.gz
1fea0b437baa7c034cf188045678eec9  E16_R2_001.fastq.gz
1da4e491aaab9f8dd31ec1cceb9df6e4  E1_R1_001.fastq.gz
f648c57469094634c5db06399de71adf  E1_R2_001.fastq.gz
bdb29d094b4e99accb4901129f5632cc  E2_R1_001.fastq.gz
7a901ba3ffb15757fa1f35d6fdd77977  E2_R2_001.fastq.gz
57003d309200f5e00bc1ae0ef7efc569  E3_R1_001.fastq.gz
945adabf4bb06a71750c8a7dcb2ffad4  E3_R2_001.fastq.gz
9c1a16d8efaf9e8b393c57e59cf582ed  E4_R1_001.fastq.gz
3ec5d310d5ca71be11fb61da8f89e16b  E4_R2_001.fastq.gz
bcbd03903528b4d23694aeac58cba29b  E5_R1_001.fastq.gz
729188fdac7ca7d5b984cd92c7cacf80  E5_R2_001.fastq.gz
03188b855cddf9c4310f2068b43461b6  E6_R1_001.fastq.gz
0a4f9de84fe32baea6f7cab4e3fa67a4  E6_R2_001.fastq.gz
a8f4c26ff35730b89ea72913a56721a7  E7_R1_001.fastq.gz
91502a50cf91652e37f305edd986774b  E7_R2_001.fastq.gz
8a6e9fc499a8fca32614cd0e0e86fb1d  E8_R1_001.fastq.gz
cd7de23f56172a399144a2dc6ae7911b  E8_R2_001.fastq.gz
e0e620f5822277987ec718839df370de  E9_R1_001.fastq.gz
871c3316046af6bbdbef3cf80ce90d91  E9_R2_001.fastq.gz
```

c) Count number of reads per file 

check for code after @ in fastq.gz files(e.g.,@GWNJ).

```
zgrep -c "@GWNJ" /data/putnamlab/KITT/hputnam/20201124_Becker_RNASeq/*.gz > raw_seq_counts

```

```
C17_R1_001.fastq.gz:11581540
C17_R2_001.fastq.gz:11581540
C18_R1_001.fastq.gz:10322369
C18_R2_001.fastq.gz:10322369
C19_R1_001.fastq.gz:11407222
C19_R2_001.fastq.gz:11407222
C20_R1_001.fastq.gz:10991646
C20_R2_001.fastq.gz:10991646
C21_R1_001.fastq.gz:7581776
C21_R2_001.fastq.gz:7581776
C22_R1_001.fastq.gz:10661421
C22_R2_001.fastq.gz:10661421
C23_R1_001.fastq.gz:9669778
C23_R2_001.fastq.gz:9669778
C24_R1_001.fastq.gz:10109664
C24_R2_001.fastq.gz:10109664
C25_R1_001.fastq.gz:13384360
C25_R2_001.fastq.gz:13384360
C26_R1_001.fastq.gz:10513807
C26_R2_001.fastq.gz:10513807
C27_R1_001.fastq.gz:11612160
C27_R2_001.fastq.gz:11612160
C28_R1_001.fastq.gz:10592999
C28_R2_001.fastq.gz:10592999
C29_R1_001.fastq.gz:8311723
C29_R2_001.fastq.gz:8311723
C30_R1_001.fastq.gz:12019128
C30_R2_001.fastq.gz:12019128
C31_R1_001.fastq.gz:11231196
C31_R2_001.fastq.gz:11231196
C32_R1_001.fastq.gz:12874230
C32_R2_001.fastq.gz:12874230
E10_R1_001.fastq.gz:13305449
E10_R2_001.fastq.gz:13305449
E11_R1_001.fastq.gz:8574604
E11_R2_001.fastq.gz:8574604
E12_R1_001.fastq.gz:13030583
E12_R2_001.fastq.gz:13030583
E13_R1_001.fastq.gz:11433773
E13_R2_001.fastq.gz:11433773
E14_R1_001.fastq.gz:10220049
E14_R2_001.fastq.gz:10220049
E15_R1_001.fastq.gz:11311009
E15_R2_001.fastq.gz:11311009
E16_R1_001.fastq.gz:13487378
E16_R2_001.fastq.gz:13487378
E1_R1_001.fastq.gz:11547975
E1_R2_001.fastq.gz:11547975
E2_R1_001.fastq.gz:10754697
E2_R2_001.fastq.gz:10754697
E3_R1_001.fastq.gz:11691716
E3_R2_001.fastq.gz:11691716
E4_R1_001.fastq.gz:9886838
E4_R2_001.fastq.gz:9886838
E5_R1_001.fastq.gz:7511858
E5_R2_001.fastq.gz:7511858
E6_R1_001.fastq.gz:11549219
E6_R2_001.fastq.gz:11549219
E7_R1_001.fastq.gz:11170722
E7_R2_001.fastq.gz:11170722
E8_R1_001.fastq.gz:10598636
E8_R2_001.fastq.gz:10598636
E9_R1_001.fastq.gz:11268722
E9_R2_001.fastq.gz:11268722
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
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/raw

module load all/FastQC/0.11.9-Java-11

for file in /data/putnamlab/KITT/hputnam/20201124_Becker_RNASeq/*.gz
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
`
rm /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/*.sam

`


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


python /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/prepDE.py -g /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/Poc_gene_count_matrix.csv -i /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/data/mapped/samplelist.txt
```

```
sbatch /data/putnamlab/hputnam/Becker_E5/RNASeq_Becker_E5/scripts/GTFtoCounts.sh
```


g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/GTF_merge/gene_count_ofav_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/
```

