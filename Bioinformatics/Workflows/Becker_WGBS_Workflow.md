Whole Genome Bisulfite Sequencing (WGBS) Workflow created by Hollie M. Putnam on 20201206 

WGBS Workflow edited by Danielle M. Becker 20210312

# Genewiz Data for WGBS

```
grep -E 'gene structure "(gene)";' input_file.txt

mkdir /data/putnamlab/hputnam/Becker_E5

cd /data/putnamlab/hputnam/Becker_E5

#Sample Information
Genewiz Sample Name | Putnam	Library Name | i7 Index Name |i7 Index Sequence |i5 Index Name | i5 Index Sequence 
---|---|---|---|---|---| ---|
1	|E3	|i7_UDI0001	|CAAGCAGAAGACGGCATACGAGATAACCGCGGGTGACTGGAGTTCAGACGTGT	|i5UDI0001	|AATGATACGGCGACCACCGAGATCTACACAGCGCTAGACACTCTTTCCCTACACGAC  
2	|E7	|i7_UDI0002	|CAAGCAGAAGACGGCATACGAGATGGTTATAAGTGACTGGAGTTCAGACGTGT	|i5UDI0002	|AATGATACGGCGACCACCGAGATCTACACGATATCGAACACTCTTTCCCTACACGAC
3	|E14	|i7_UDI0003	|CAAGCAGAAGACGGCATACGAGATCCAAGTCCGTGACTGGAGTTCAGACGTGT	|i5UDI0003	|AATGATACGGCGACCACCGAGATCTACACCGCAGACGACACTCTTTCCCTACACGAC
4	|C26	|i7_UDI0004	|CAAGCAGAAGACGGCATACGAGATTTGGACTTGTGACTGGAGTTCAGACGTGT	|i5UDI0004	|AATGATACGGCGACCACCGAGATCTACACTATGAGTAACACTCTTTCCCTACACGAC
5	|C31	|i7_UDI0005	|CAAGCAGAAGACGGCATACGAGATCAGTGGATGTGACTGGAGTTCAGACGTGT	|i5UDI0005	|AATGATACGGCGACCACCGAGATCTACACAGGTGCGTACACTCTTTCCCTACACGAC
6	|C29	|i7_UDI0006	|CAAGCAGAAGACGGCATACGAGATTGACAAGCGTGACTGGAGTTCAGACGTGT	|i5UDI0006	|AATGATACGGCGACCACCGAGATCTACACGAACATACACACTCTTTCCCTACACGAC
7	|E13	|i7_UDI0007	|CAAGCAGAAGACGGCATACGAGATCTAGCTTGGTGACTGGAGTTCAGACGTGT	|i5UDI0007	|AATGATACGGCGACCACCGAGATCTACACACATAGCGACACTCTTTCCCTACACGAC
8	|C20	|i7_UDI0008	|CAAGCAGAAGACGGCATACGAGATTCGATCCAGTGACTGGAGTTCAGACGTGT	|i5UDI0008	|AATGATACGGCGACCACCGAGATCTACACGTGCGATAACACTCTTTCCCTACACGAC
9	|C18	|i7_UDI0009	|CAAGCAGAAGACGGCATACGAGATCCTGAACTGTGACTGGAGTTCAGACGTGT	|i5UDI0009	|AATGATACGGCGACCACCGAGATCTACACCCAACAGAACACTCTTTCCCTACACGAC
10	|C24	|i7_UDI0010	|CAAGCAGAAGACGGCATACGAGATTTCAGGTCGTGACTGGAGTTCAGACGTGT	|i5UDI0010	|AATGATACGGCGACCACCGAGATCTACACTTGGTGAGACACTCTTTCCCTACACGAC
11	|E10	|i7_UDI0011	|CAAGCAGAAGACGGCATACGAGATAGTAGAGAGTGACTGGAGTTCAGACGTGT	|i5UDI0011	|AATGATACGGCGACCACCGAGATCTACACCGCGGTTCACACTCTTTCCCTACACGAC
12	|C27	|i7_UDI0012	|CAAGCAGAAGACGGCATACGAGATGACGAGAGGTGACTGGAGTTCAGACGTGT	|i5UDI0012	|AATGATACGGCGACCACCGAGATCTACACTATAACCTACACTCTTTCCCTACACGAC
13	|C22	|i7_UDI0013	|CAAGCAGAAGACGGCATACGAGATAGACTTGGGTGACTGGAGTTCAGACGTGT	|i5UDI0013	|AATGATACGGCGACCACCGAGATCTACACAAGGATGAACACTCTTTCCCTACACGAC
14	|E11	|i7_UDI0014	|CAAGCAGAAGACGGCATACGAGATGAGTCCAAGTGACTGGAGTTCAGACGTGT	|i5UDI0014	|AATGATACGGCGACCACCGAGATCTACACGGAAGCAGACACTCTTTCCCTACACGAC
15	|E4	|i7_UDI0015	|CAAGCAGAAGACGGCATACGAGATCTTAAGCCGTGACTGGAGTTCAGACGTGT	|i5UDI0015	|AATGATACGGCGACCACCGAGATCTACACTCGTGACCACACTCTTTCCCTACACGAC
16	|C28	|i7_UDI0016	|CAAGCAGAAGACGGCATACGAGATTCCGGATTGTGACTGGAGTTCAGACGTGT	|i5UDI0016	|AATGATACGGCGACCACCGAGATCTACACCTACAGTTACACTCTTTCCCTACACGAC
17	|E6	|i7_UDI0017	|CAAGCAGAAGACGGCATACGAGATCTGTATTAGTGACTGGAGTTCAGACGTGT	|i5UDI0017	|AATGATACGGCGACCACCGAGATCTACACATATTCACACACTCTTTCCCTACACGAC
18	|C21	|i7_UDI0018	|CAAGCAGAAGACGGCATACGAGATTCACGCCGGTGACTGGAGTTCAGACGTGT	|i5UDI0018	|AATGATACGGCGACCACCGAGATCTACACGCGCCTGTACACTCTTTCCCTACACGAC
19	|E8	|i7_UDI0019	|CAAGCAGAAGACGGCATACGAGATACTTACATGTGACTGGAGTTCAGACGTGT	|i5UDI0019	|AATGATACGGCGACCACCGAGATCTACACACTCTATGACACTCTTTCCCTACACGAC
20	|E12	|i7_UDI0020	|CAAGCAGAAGACGGCATACGAGATGTCCGTGCGTGACTGGAGTTCAGACGTGT	|i5UDI0020	|AATGATACGGCGACCACCGAGATCTACACGTCTCGCAACACTCTTTCCCTACACGAC
21	|E5	|i7_UDI0021	|CAAGCAGAAGACGGCATACGAGATAAGGTACCGTGACTGGAGTTCAGACGTGT	|i5UDI0021	|AATGATACGGCGACCACCGAGATCTACACAAGACGTCACACTCTTTCCCTACACGAC
22	|C30	|i7_UDI0022	|CAAGCAGAAGACGGCATACGAGATGGAACGTTGTGACTGGAGTTCAGACGTGT	|i5UDI0022	|AATGATACGGCGACCACCGAGATCTACACGGAGTACTACACTCTTTCCCTACACGAC
23	|E1	|i7_UDI0023	|CAAGCAGAAGACGGCATACGAGATAATTCTGCGTGACTGGAGTTCAGACGTGT	|i5UDI0023	|AATGATACGGCGACCACCGAGATCTACACACCGGCCAACACTCTTTCCCTACACGAC
24	|C23	|i7_UDI0024	|CAAGCAGAAGACGGCATACGAGATGGCCTCATGTGACTGGAGTTCAGACGTGT	|i5UDI0024	|AATGATACGGCGACCACCGAGATCTACACGTTAATTGACACTCTTTCCCTACACGAC
25	|E16	|i7_UDI0025	|CAAGCAGAAGACGGCATACGAGATATCTTAGTGTGACTGGAGTTCAGACGTGT	|i5UDI0025	|AATGATACGGCGACCACCGAGATCTACACAACCGCGGACACTCTTTCCCTACACGAC
26	|C19	|i7_UDI0026	|CAAGCAGAAGACGGCATACGAGATGCTCCGACGTGACTGGAGTTCAGACGTGT	|i5UDI0026	|AATGATACGGCGACCACCGAGATCTACACGGTTATAAACACTCTTTCCCTACACGAC
27	|C25	|i7_UDI0027	|CAAGCAGAAGACGGCATACGAGATATACCAAGGTGACTGGAGTTCAGACGTGT	|i5UDI0027	|AATGATACGGCGACCACCGAGATCTACACCCAAGTCCACACTCTTTCCCTACACGAC
28	|E2|	i7_UDI0028	|CAAGCAGAAGACGGCATACGAGATGCGTTGGAGTGACTGGAGTTCAGACGTGT	|i5UDI0028	|AATGATACGGCGACCACCGAGATCTACACTTGGACTTACACTCTTTCCCTACACGAC
29	|E9	|i7_UDI0029	|CAAGCAGAAGACGGCATACGAGATCTTCACGGGTGACTGGAGTTCAGACGTGT	|i5UDI0029	|AATGATACGGCGACCACCGAGATCTACACCAGTGGATACACTCTTTCCCTACACGAC
30	|C32	|i7_UDI0030	|CAAGCAGAAGACGGCATACGAGATTCCTGTAAGTGACTGGAGTTCAGACGTGT	|i5UDI0030	|AATGATACGGCGACCACCGAGATCTACACTGACAAGCACACTCTTTCCCTACACGAC
31	|E15	|i7_UDI0031	|CAAGCAGAAGACGGCATACGAGATAGAATGCCGTGACTGGAGTTCAGACGTGT	|i5UDI0031	|AATGATACGGCGACCACCGAGATCTACACCTAGCTTGACACTCTTTCCCTACACGAC
32	|C17	|i7_UDI0032	|CAAGCAGAAGACGGCATACGAGATGAGGCATTGTGACTGGAGTTCAGACGTGT	|i5UDI0032	|AATGATACGGCGACCACCGAGATCTACACTCGATCCAACACTCTTTCCCTACACGAC
```


# 1) Obtain Pocillopora verrucosa genome
[Buitrago-López et al 2020](https://academic.oup.com/gbe/article/12/10/1911/5898631)

## Location of Pocillopora verrucosa genome
```
cd /data/putnamlab/REFS/Pverr/
```
```
wget http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.fasta.gz
gunzip Pver_genome_assembly_v1.0.fasta.gz
```

# 2) Check Data Transfer Fidelity
```
raw data - read only
/data/putnamlab/KITT/hputnam/20201206_Becker_WGBS
```

```
md5sum *.gz > URI_md5sum.txt
```
```
Genewiz md5 | sample id | URI md5 
---|---|---|
de7c178629f722bbe869627bf580702a  |`1_R1_001.fastq.gz`|de7c178629f722bbe869627bf580702a
dea71bd12349e979cccdb18c634141fc  |`1_R2_001.fastq.gz`|dea71bd12349e979cccdb18c634141fc
c4aa994b3d3d3633a4704a183f983d03  |`2_R2_001.fastq.gz`|97e6e9b0905730b59a64104a1f41f1cc
97e6e9b0905730b59a64104a1f41f1cc  |`2_R1_001.fastq.gz`|c4aa994b3d3d3633a4704a183f983d03
f5de1954e3167522b964ca28144a1a80  |`3_R1_001.fastq.gz`|f5de1954e3167522b964ca28144a1a80
4ba29015cf612fd8c498f4d06c824424  |`3_R2_001.fastq.gz`|4ba29015cf612fd8c498f4d06c824424
e5bf0fff9077cf3ca35e39da68673df5  |`4_R1_001.fastq.gz`|e5bf0fff9077cf3ca35e39da68673df5
5a9cad6d0db9101b85f68a0da6e979ce  |`4_R2_001.fastq.gz`|5a9cad6d0db9101b85f68a0da6e979ce
0b9c1781483a7f258b2798d77375e75b  |`5_R1_001.fastq.gz`|0b9c1781483a7f258b2798d77375e75b
49444ba373b827bf01be0378268a5b64  |`5_R2_001.fastq.gz`|49444ba373b827bf01be0378268a5b64
6fd2b5a5a80dc71366d74dc8d6ff7ba7  |`6_R1_001.fastq.gz`|6fd2b5a5a80dc71366d74dc8d6ff7ba7
32e52c56f27445fe2b2c4826e86d4b1f  |`6_R2_001.fastq.gz`|32e52c56f27445fe2b2c4826e86d4b1f
5a4074e0705e44c2e6be63f7a8875563  |`7_R1_001.fastq.gz`|5a4074e0705e44c2e6be63f7a8875563
3ffb1c7d82d277b588977f091dc1257e  |`7_R2_001.fastq.gz`|3ffb1c7d82d277b588977f091dc1257e
63b0090ba62038db826b1533f0480d42  |`8_R1_001.fastq.gz`|63b0090ba62038db826b1533f0480d42
e4eb89a4197c0d83a17c7defcf0905e2  |`8_R2_001.fastq.gz`|e4eb89a4197c0d83a17c7defcf0905e2
04b9abf33a1094505d26fdce46c4f516  |`9_R1_001.fastq.gz`|04b9abf33a1094505d26fdce46c4f516
29fa0cb13b63e3670ac0a613dc9bffdb  |`9_R2_001.fastq.gz`|29fa0cb13b63e3670ac0a613dc9bffdb
5945c5f9e3a318dcb5deb39798868889  |`10_R1_001.fastq.gz`|5945c5f9e3a318dcb5deb39798868889
1bd20839daabd95cdef4672a7bb8e3bd  |`10_R2_001.fastq.gz`|1bd20839daabd95cdef4672a7bb8e3bd
42f8fb9d54ec859b5580cea0d8f49986  |`11_R1_001.fastq.gz`|42f8fb9d54ec859b5580cea0d8f49986
f12f044928f1d32a83eadd61d434879a  |`11_R2_001.fastq.gz`|f12f044928f1d32a83eadd61d434879a
cdf14cbb48a545f750d5c1d6dbfe838d  |`12_R1_001.fastq.gz`|cdf14cbb48a545f750d5c1d6dbfe838d
74c06455a1144bbac680f38414840c8b  |`12_R2_001.fastq.gz`|74c06455a1144bbac680f38414840c8b
7b4f92dec737564676962f535200f5a3  |`13_R1_001.fastq.gz`|7b4f92dec737564676962f535200f5a3
347f9c40db6a45abb2622c4e2999e02b  |`13_R2_001.fastq.gz`|347f9c40db6a45abb2622c4e2999e02b
60dbe8aaed56c4931914c7917cd077f5  |`14_R1_001.fastq.gz`|60dbe8aaed56c4931914c7917cd077f5
791f783dba458ce2d613108326bc5896  |`14_R2_001.fastq.gz`|791f783dba458ce2d613108326bc5896
6aef447ea8e0d3e45d8c652a2e925c08  |`15_R1_001.fastq.gz`|6aef447ea8e0d3e45d8c652a2e925c08
ab2cf6902b81a0a1f5c14c628be30248  |`15_R2_001.fastq.gz`|ab2cf6902b81a0a1f5c14c628be30248
208c03e26b5bbcf53b216263dfbcf5f0  |`16_R1_001.fastq.gz`|208c03e26b5bbcf53b216263dfbcf5f0
5403deb7a40c71476e45cb6228317b43  |`16_R2_001.fastq.gz`|5403deb7a40c71476e45cb6228317b43
c16aed82de1d33f97e178b92fccbc592  |`17_R1_001.fastq.gz`|c16aed82de1d33f97e178b92fccbc592
8913db9656dfe652ada9640cddf8bcf9  |`17_R2_001.fastq.gz`|8913db9656dfe652ada9640cddf8bcf9
1d1538ff62375d7bdc863b608b19296e  |`18_R1_001.fastq.gz`|1d1538ff62375d7bdc863b608b19296e
9a0b3623b76e206d2ea8417ec1fe84a3  |`18_R2_001.fastq.gz`|9a0b3623b76e206d2ea8417ec1fe84a3
5c21c1041b7e23a15b92519218dfb645  |`19_R1_001.fastq.gz`|5c21c1041b7e23a15b92519218dfb645
794b42ba3ef16abd1fd7db6ad9a63cb5  |`19_R2_001.fastq.gz`|794b42ba3ef16abd1fd7db6ad9a63cb5
5ae1238518d585c76d02008decd9e852  |`20_R1_001.fastq.gz`|5ae1238518d585c76d02008decd9e852
f624ad50bf05353dea0eee3d00828541  |`20_R2_001.fastq.gz`|f624ad50bf05353dea0eee3d00828541
1b6182d26cb3990a8c4c16704a74366f  |`21_R1_001.fastq.gz`|1b6182d26cb3990a8c4c16704a74366f
0a096ccce7a4074d0053b9a69543772b  |`21_R2_001.fastq.gz`|0a096ccce7a4074d0053b9a69543772b
cf56d401e81a7205e6a4cd20081dad66  |`22_R1_001.fastq.gz`|cf56d401e81a7205e6a4cd20081dad66
b404dbe889f20d8f8c156d06919fe379  |`22_R2_001.fastq.gz`|b404dbe889f20d8f8c156d06919fe379
fdee7108b5feb150117c7ceb3479c7ec  |`23_R1_001.fastq.gz`|fdee7108b5feb150117c7ceb3479c7ec
bb43c05f5fd1d26513d3bec42fc38c9d  |`23_R2_001.fastq.gz`|bb43c05f5fd1d26513d3bec42fc38c9d
0420ecbe17e3a347acc4de879db23caf  |`24_R1_001.fastq.gz`|0420ecbe17e3a347acc4de879db23caf
04375573c2dfca60efd104b3d533bc38  |`24_R2_001.fastq.gz`|04375573c2dfca60efd104b3d533bc38
8538fd3f3c5a0b28aac25a74880a0343  |`25_R1_001.fastq.gz`|8538fd3f3c5a0b28aac25a74880a0343
0c0f94d7849d7b880ae6d7a7ca1d78b9  |`25_R2_001.fastq.gz`|0c0f94d7849d7b880ae6d7a7ca1d78b9
f950fe6e468ebe288d6f3b499b726e1c  |`26_R1_001.fastq.gz`|f950fe6e468ebe288d6f3b499b726e1c
4d72a6412e2fd409017100455e0a5dac  |`26_R2_001.fastq.gz`|4d72a6412e2fd409017100455e0a5dac
88287eba1a92a9c694aeee9dfcf37e82  |`27_R1_001.fastq.gz`|88287eba1a92a9c694aeee9dfcf37e82
c9ff359677ded2aaa9f179681cc2cba5  |`27_R2_001.fastq.gz`|c9ff359677ded2aaa9f179681cc2cba5
9d9a8f0c4197243aac7e053857563058  |`28_R1_001.fastq.gz`|9d9a8f0c4197243aac7e053857563058
71c90daaddade7271bb793d92972510e  |`28_R2_001.fastq.gz`|71c90daaddade7271bb793d92972510e
a72a65f46f1c89a08872e8a7ef6f171c  |`29_R1_001.fastq.gz`|a72a65f46f1c89a08872e8a7ef6f171c
2f9a3bf3218b17160de14da31485a394  |`29_R2_001.fastq.gz`|2f9a3bf3218b17160de14da31485a394
1fb026a8aa9a5c22b94e858ff68b89fd  |`30_R1_001.fastq.gz`|1fb026a8aa9a5c22b94e858ff68b89fd
d6a96b7d73725da74ebfe76b865b47d3  |`30_R2_001.fastq.gz`|d6a96b7d73725da74ebfe76b865b47d3
7ed718e6214db23df2d415f887d35e25  |`31_R1_001.fastq.gz`|7ed718e6214db23df2d415f887d35e25
94a0c47fea951786bf06f93e3101f9a8  |`31_R2_001.fastq.gz`|94a0c47fea951786bf06f93e3101f9a8
0be0e18cebf3b0c52711a1cd15535bb8  |`32_R1_001.fastq.gz`|0be0e18cebf3b0c52711a1cd15535bb8
92afdc0db53a918b1bfd15c2c2ff4943  |`32_R2_001.fastq.gz`|92afdc0db53a918b1bfd15c2c2ff4943
```
## Count the number of reads per sample
```
zcat *fastq.gz | echo $((`wc -l`/4)) > rawread.counts.txt
This counts reads in goups of 4 lines per read
This should match with the Genewiz summary
```

```
cd /data/putnamlab/hputnam/Becker_E5
mkdir WGBS_Becker_E5
cd WGBS_Becker_E5
```
```
/data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0.fasta
```

# 3) Run NF Core MethylSeq

[NF Core Methylseq](https://github.com/nf-core/methylseq/)  
NEXTFLOW  ~  version 20.04.1  
nf-core/methylseq v1.5  

```
nano /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/scripts/methylseq.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed

module load Nextflow

# run nextflow methylseq

nextflow run nf-core/methylseq -profile singularity \
--aligner bismark \
--fasta /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0.fasta \
--save_reference \
--reads '/data/putnamlab/KITT/hputnam/20201206_Becker_WGBS/*_R{1,2}_001.fastq.gz' \
--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_r1 10 --three_prime_clip_r2 10 \
--non_directional \
--cytosine_report \
--relax_mismatches \
--unmapped \
--outdir Becker_WGBS \
-name Becker_WGBS_Poc
```

```
sbatch /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/scripts/methylseq.sh
```
- slurm-16249.out


### Run timed out

 - resumed run with longer time using code below in slurm-16594.out 

```
nano /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/scripts/methylseq_resume.sh
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed

module load Nextflow

# run nextflow methylseq

nextflow run nf-core/methylseq -resume \
-profile singularity \
--aligner bismark \
--fasta /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0.fasta \
--save_reference \
--reads '/data/putnamlab/KITT/hputnam/20201206_Becker_WGBS/*_R{1,2}_001.fastq.gz' \
--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_r1 10 --three_prime_clip_r2 10 \
--non_directional \
--cytosine_report \
--relax_mismatches \
--unmapped \
--outdir Becker_WGBS

```

```
sbatch /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/scripts/methylseq_resume.sh
```


## Library E8 - sample 19 failed
 - has weird library results from tapestation
 - has lower quality than all other libraries
 - has high GC content
 - has low duplication content
 - has very high adapter content
 
 It seems this library failed at the prep step. perhaps too little or poor quality DNA?
 We will move ahead with the other 31 libraries for analysis.


Library 2 (sample E7) and Library 16 (Sample C28) have low coverage and high duplication in MethylSeq MultiQc report

# 4) Merge strands 

The Bismark [coverage2cytosine](https://github.com/FelixKrueger/Bismark/blob/master/coverage2cytosine) command re-reads the genome-wide report and merges methylation evidence of both top and bottom strand.

```
nano /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/scripts/cov_to_cyto.sh
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH -D /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/
#SBATCH -p putnamlab
#SBATCH --cpus-per-task=3

# load modules needed

module load Bismark/0.20.1-foss-2018b

# run coverage2cytosine merge of strands

 find /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/Becker_WGBS/bismark_methylation_calls/methylation_coverage/*deduplicated.bismark.cov.gz \
 | xargs basename -s _R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz \
 | xargs -I{} coverage2cytosine \
 --genome_folder /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/Becker_WGBS/reference_genome/BismarkIndex \
 -o {} \
 --merge_CpG \
 --zero_based \
 /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/Becker_WGBS/bismark_methylation_calls/methylation_coverage/{}_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz

```

```
sbatch /data/putnamlab/hputnam/Becker_E5/WGBS_Becker_E5/scripts/cov_to_cyto.sh
```

## Scaffolds were not organized in same arrangemnet within tab files, needed to take merged files and sort them before for loop, tested on two files and it worked!


```
Test steps to make sure this worked on two files before writing larger for loop:

#make sorted test .cov files

bedtools sort -i 10.CpG_report.merged_CpG_evidence.cov > 10.CpG_report.merged_CpG_evidence_sorted.cov

bedtools sort -i 11.CpG_report.merged_CpG_evidence.cov > 11.CpG_report.merged_CpG_evidence_sorted.cov

#run loop to filter CpGs for 10x sorted test coverage files

for f in *merged_CpG_evidence_sorted.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence_sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_10x_sorted_test.tab
done

#Create a file with positions found in all samples at specified coverage for two test files 

module load BEDTools/2.27.1-foss-2018b

multiIntersectBed -i \
/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/10_10x_sorted_test.tab  /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/11_10x_sorted_test.tab   > /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/CpG.2samps.10x_sorted_test.bed

cat /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/CpG.2samps.10x_sorted_test.bed | awk '$4 ==2' > /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/CpG.filt.2samps.10x_sorted_test.bed 

The test worked! Going to sort all .tab files so scaffold order is consistent.

```

## Sorting the merged files so scaffolds are all in the same order and multiIntersectBed will run correctly

```
#run for loop using bedtools to sort all .tab files

nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/scripts/bedtools.sort.sh

#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/
#SBATCH --cpus-per-task=3

module load BEDTools/2.27.1-foss-2018b

for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  bedtools sort -i "${f}" \
  > "${STEM}"_sorted.cov
done

sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/bedtools.sort.sh
Submitted batch job 1893831

```

# 5) Create files for statistical analysis

## Above commands were run on Hollies server folder, further commands were run on mine
## Run loop to filter CpGs for 5x coverage, creating tab files with raw count for glms

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_5x_sorted.tab
done
```
## Run loop to filter CpGs for 10x coverage, creating tab files with raw count for glms

```

for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_10x_sorted.tab
done
```
```
wc -l *5x_sorted.tab

5650793 10_5x_sorted.tab
4445225 11_5x_sorted.tab
6214291 12_5x_sorted.tab
4537224 13_5x_sorted.tab
5792624 14_5x_sorted.tab
5310753 15_5x_sorted.tab
5966510 1_5x_sorted.tab
2247820 16_5x_sorted.tab
5538829 17_5x_sorted.tab
4581541 18_5x_sorted.tab
   6632 19_5x_sorted.tab
5177825 20_5x_sorted.tab
6711526 21_5x_sorted.tab
4912535 22_5x_sorted.tab
5606601 23_5x_sorted.tab
5173314 24_5x_sorted.tab
6130324 25_5x_sorted.tab
 541141 2_5x_sorted.tab
5705061 26_5x_sorted.tab
6983120 27_5x_sorted.tab
6258204 28_5x_sorted.tab
6021709 29_5x_sorted.tab
4558106 30_5x_sorted.tab
5705620 31_5x_sorted.tab
6790513 32_5x_sorted.tab
5687903 3_5x_sorted.tab
5580876 4_5x_sorted.tab
4123978 5_5x_sorted.tab
4254063 6_5x_sorted.tab
5736001 7_5x_sorted.tab
4843732 8_5x_sorted.tab
5138423 9_5x_sorted.tab
161932817 total


wc -l *10x_sorted.tab

2545720 10_10x_sorted.tab
3012855 1_10x_sorted.tab
1201637 11_10x_sorted.tab
3483529 12_10x_sorted.tab
1315708 13_10x_sorted.tab
2903726 14_10x_sorted.tab
2098508 15_10x_sorted.tab
 351245 16_10x_sorted.tab
2418593 17_10x_sorted.tab
1393271 18_10x_sorted.tab
   2389 19_10x_sorted.tab
2249149 20_10x_sorted.tab
  48301 2_10x_sorted.tab
4451413 21_10x_sorted.tab
1663741 22_10x_sorted.tab
2505402 23_10x_sorted.tab
2074085 24_10x_sorted.tab
3450513 25_10x_sorted.tab
2793280 26_10x_sorted.tab
5113183 27_10x_sorted.tab
3634844 28_10x_sorted.tab
3268323 29_10x_sorted.tab
1380396 30_10x_sorted.tab
2640988 3_10x_sorted.tab
2472008 31_10x_sorted.tab
4433221 32_10x_sorted.tab
2454172 4_10x_sorted.tab
1043789 5_10x_sorted.tab
1128028 6_10x_sorted.tab
2681849 7_10x_sorted.tab
1572925 8_10x_sorted.tab
1943521 9_10x_sorted.tab
73730312 total

```


### Samples 19, 16, and 2 have low data coverage, not moving forward with these samples in downstream steps


# 6) Create a file with positions found in all samples at specified coverage


```
nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/5x_intersect.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/
#SBATCH --cpus-per-task=3


# load modules needed

module load BEDTools/2.27.1-foss-2018b

multiIntersectBed -i \
/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/1_5x_sorted.tab  /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/3_5x_sorted.tab  /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/4_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/5_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/6_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/7_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/8_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/9_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/10_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/11_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/12_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/13_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/14_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/15_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/17_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/18_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/20_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/21_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/22_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/23_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/24_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/25_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/26_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/27_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/28_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/29_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/30_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/31_5x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/32_5x_sorted.tab > /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/CpG.29samps.5x_sorted.bed

cat /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/CpG.29samps.5x_sorted.bed | awk '$4 ==29' > /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/CpG.filt.29samps.5x_sorted.bed 

``` 

```
sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/5x_intersect.sh
Submitted batch job 1896854
```


```
nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/10x_intersect.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/
#SBATCH --cpus-per-task=3


# load modules needed

module load BEDTools/2.27.1-foss-2018b

multiIntersectBed -i \
/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/1_10x_sorted.tab  /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/3_10x_sorted.tab  /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/4_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/5_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/6_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/7_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/8_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/9_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/10_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/11_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/12_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/13_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/14_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/15_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/17_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/18_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/20_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/21_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/22_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/23_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/24_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/25_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/26_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/27_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/28_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/29_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/30_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/31_10x_sorted.tab /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/32_10x_sorted.tab > /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/CpG.29samps.10x_sorted.bed

cat /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/CpG.29samps.10x_sorted.bed | awk '$4 ==29' > /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/CpG.filt.29samps.10x_sorted.bed 

``` 

```
sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/10x_intersect.sh
Submitted batch job 1896853

```

 

# 7) Create bedgraphs post merge

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4}}' \
  > "${STEM}"_10x_sorted.bedgraph
done

```
```

for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${STEM}"_5x_sorted.bedgraph
done

```

# 8) Use intersectBed to find where loci and genes intersect, allowing loci to be mapped to annotated genes


### Modified Pver_genome_assembly file to only include genes in R markdown on desktop

```
#filtered gff3 to only include gene positions modified.gff3 > gene.gff3

awk '{if ($3 == "gene") {print}}' /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0_modified.gff3 > /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0.gene.gff3

```

# Use intersectBed to find where loci and genes intersect, allowing loci to be mapped to annotated genes

```
nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/5x_intersectBed.sh
```

```
## wb: Print all lines in the second file
## a: file that ends in pos Only
## b: annotated gene list


#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/
#SBATCH --cpus-per-task=3


# load modules needed

module load BEDTools/2.27.1-foss-2018b


for i in *5x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0.gene.gff3 \
  > ${i}_gene
done


sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/5x_intersectBed.sh
Submitted batch job 1897057

```
```
nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/10x_intersectBed.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/
#SBATCH --cpus-per-task=3


# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *10x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b /data/putnamlab/REFS/Pverr/Pver_genome_assembly_v1.0.gene.gff3 \
  > ${i}_gene
done

sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/10x_intersectBed.sh
Submitted batch job 1897058

```

## Intersect with file to subset only those positions found in all samples

```
nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/5x_intersect_final.sh
```
```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/
#SBATCH --cpus-per-task=3


# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *5x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.29samps.5x_sorted.bed \
  > ${i}_CpG_5x_enrichment.bed
done
```
```
sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/5x_intersect_final.sh
Submitted batch job 1900070
```


```
nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/10x_intersect_final.sh
```
```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu 
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/
#SBATCH --cpus-per-task=3


# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *10x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.29samps.10x_sorted.bed \
  > ${i}_CpG_10x_enrichment.bed
done
```

```
sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/10x_intersect_final.sh
Submitted batch job 1900071
```

```
wc -l *5x_enrichment.bed

477425 10_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 11_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 12_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 13_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 14_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 15_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 1_5x_sorted.tab_gene_CpG_5x_enrichment.bed
288120 16_5x_sorted.tab_gene_CpG_5x_enrichment.bed #not using for downstream analysis, low coverage 
477425 17_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 18_5x_sorted.tab_gene_CpG_5x_enrichment.bed
    1077 19_5x_sorted.tab_gene_CpG_5x_enrichment.bed #not using for downstream analysis, low coverage 
477425 20_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 21_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 22_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 23_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 24_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 25_5x_sorted.tab_gene_CpG_5x_enrichment.bed
    94912 2_5x_sorted.tab_gene_CpG_5x_enrichment.bed #not using for downstream analysis, low coverage 
477425 26_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 27_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 28_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 29_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 30_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 31_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 32_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 3_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 4_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 5_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 6_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 7_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 8_5x_sorted.tab_gene_CpG_5x_enrichment.bed
477425 9_5x_sorted.tab_gene_CpG_5x_enrichment.bed
  14229434 total



wc -l *10x_enrichment.bed

38637 10_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 1_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 11_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 12_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 13_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 14_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 15_10x_sorted.tab_gene_CpG_10x_enrichment.bed
25858 16_10x_sorted.tab_gene_CpG_10x_enrichment.bed #not using for downstream analysis, low coverage 
38637 17_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 18_10x_sorted.tab_gene_CpG_10x_enrichment.bed
    210 19_10x_sorted.tab_gene_CpG_10x_enrichment.bed #not using for downstream analysis, low coverage 
38637 20_10x_sorted.tab_gene_CpG_10x_enrichment.bed
    9866 2_10x_sorted.tab_gene_CpG_10x_enrichment.bed #not using for downstream analysis, low coverage
38637 21_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 22_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 23_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 24_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 25_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 26_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 27_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 28_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 29_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 30_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 3_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 31_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 32_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 4_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 5_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 6_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 7_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 8_10x_sorted.tab_gene_CpG_10x_enrichment.bed
38637 9_10x_sorted.tab_gene_CpG_10x_enrichment.bed
  1156407 total

```
# 9) Download final .bed files to desktop for statistical analysis! :)

```
scp -r danielle_becker@bluewaves.uri.edu:/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/CovtoCyto/*_enrichment.bed /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RAnalysis/Data/WGBS
```
