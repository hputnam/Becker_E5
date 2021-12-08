---
layout: post
title: Functional Annotation Pipeline
Author: Danielle Becker-Polinski
Last Updated: 2021/10/26
tags: [ Protocol, annotation, RNASeq, GO, KEGG ]
---


## Overview

Testing out new approaches for functional annotation of  *Pocillopora verrucosa*. Previous analyses (Erin and Jill) used the program [DIAMOND BLAST](http://www.diamondsearch.org/index.php) against the NCBI nr database. Like regular BLAST, DIAMOND is a sequence aligner for nucleotide and protein sequences; unlike BLAST, it is optimized for a higher performance capability of large datasets at 100x-20,000x speed of BLAST.  The output .xml file was then funneled to [BLAST2GO](https://www.blast2go.com) (B2G) which is a bioinformatics tools for functional annotation and analysis of gene or protein sequences. It was originally developed to provide a user-friendly interface for GO annotation and now hosts many different functional annotation tools. The B2G software makes it possible to generate annotation without requiring writing any code. While B2G can do an extensive array of functions, this analysis primarily utilizes the GO mapping and annotation functions.

Through a dense literature search, multiple databases should be used for BLAST to expand the hits possible for your sequences [Baumgarten et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4586855/), [Bhattacharya et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878875/), [Cunning et al. 2018](https://www.nature.com/articles/s41598-018-34459-8?proof=t+target%3D), and [Buitrago-López et al. 2020](https://academic.oup.com/gbe/article/12/10/1911/5898631). The main concensus from this literature search is that the protein sequences should be searched against the SwissProt, TrEMBL, NCBI nr databases using BLASTp (Basic Local Alignment Search Tool, e-value cut-off = 1e-05) and retaining annotations from databases in this order. Then, BLAST2GO should be used to provide GO annotations, and KEGG, Pfam, InterProScan, should be searched to further annotated gene sets.

Followed same approach for Funtional Annotation used by [Baumgarten et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4586855/), [Bhattacharya et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878875/), and [Buitrago-López et al. 2020](https://academic.oup.com/gbe/article/12/10/1911/5898631). I generally followed the same workflow found on the [Buitrago-López et al. 2020](https://academic.oup.com/gbe/article/12/10/1911/5898631) [GitHub](https://github.com/Carol-Symbiomics/Pocillopora-verrucosa-genome/blob/master/Scripts/03.gene.models.prediction.and.annotation.sh) page for using blastp on our protein sequences against multiple databases.


### Step 1: Obtain sequences of interest.

In order to conduct functional annotation steps, protein and transcript sequences are needed. There are two main hubs where coral genomic information is stored: [Reef Genomics](http://reefgenomics.org) and [NCBI](https://www.ncbi.nlm.nih.gov). Other researchers store their genomic infomation on their own personal webpages. Genomic information must be downloaded from one of these databases in order to proceed.

#### i) Identify species to work with.

For this project, the coral species of interest is *Pocillopora verrucosa*.

#### ii) Download genomic files for species of interest.

##### [*Pocillopora verrucosa* Full Transcripts ](http://pver.reefgenomics.org/download/)

##### [*Pocillopora verrucosa* Gene models (protein)](http://pver.reefgenomics.org/download/)

Ready to start annotating!

### Step 2: Identify homologous sequences

Homology refers to the similarity of structure or genes in different taxa due to shared ancestry. {example}

Sequence homology is the homology between DNA, RNA, and protein sequences in terms of shared ancestry. Sequence homology is usually inferred by the similarity of nucleotide or amino acid sequences. Strong sequence similarity (or percent homology) provides evidence that two or more sequences are related through shared ancestry.

Several software programs have the ability to compare sequences to assess homology, the most popular one being NCBI's [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (Basic Local Alignment Search Tool - what a great name). BLAST compares nucleotide or protein sequences of interest (called a query) to their sequence databases to find regions of similarity. If a nucleotide or protein sequence of interest significantly matches a sequence/sequences in the databases, BLAST will tag the sequence of interest with the information about the known sequence(s). For example, if a new transcript sequence is identified in a mouse, BLAST could be used to see if any other animals carry a similar sequence, and if they do, what the biological functions of the sequence are.

Other software programs include [SWISS-PROT](https://iop.vast.ac.vn/theor/conferences/smp/1st/kaminuma/SWISSPROT/access.html) which is a curated protein sequence database that provides a high level of annotation (such as the description of the function of a protein, its domain structure, post-translational modifications, variants, etc), a minimal level of redundancy and a high level of integration with other databases. Recent developments of the database include: an increase in the number and scope of model organisms; cross-references to seven additional databases; a variety of new documentation files; the creation of [TREMBL](http://www3.cmbi.umcn.nl/wiki/index.php/TrEMBL), an unannotated supplement to SWISS-PROT.

### Step 3: Download your databases

#### i) On the Andromeda server, download Swiss-Prot database from [UniProt/Swiss-Prot](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz).

I created a script on Andromeda to download the updated Swiss-Prot database, unzip it, and make it into a blast database. I followed the general instructions on how to create code to download database [here](https://www.biostars.org/p/354449/).  Current databases found [here](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/).

```
pwd /data/putnamlab/shared/sbatch_executables/download_swissprot_database.sh
```

Full script:

```
#!/bin/bash
#SBATCH --job-name="ref"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu # CHANGE EMAIL
#SBATCH -D /data/putnamlab/shared/databases/swiss_db

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b

cd databases/swiss_db

echo "Making swissprot database" $date
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -parse_seqids -dbtype prot -out swissprot_20211022
echo "STOP" $(date)

```

**It is normal to get multiple files per blast database. That is how makeblastdb is supposed to work. Just make sure files for a database stay together in the same directory and you use the "basename" for the database (a suggestion: name your database some thing other than your input file name) when you run your searches.**

**Options for makeblastdb which is an application that produces BLAST databases from FASTA files, other options [here](http://nebc.nerc.ac.uk/nebc_website_frozen/nebc.nerc.ac.uk/bioinformatics/documentation/blast+/user_manual.pdf)**

- ```makeblastdb``` - application produces BLAST databases from FASTA files
- ```in``` - database that will be used to search against
- ```parse_sequids``` - parse the seq-id(s) in the FASTA input provided
- ```dbtype``` - molecule type of target db ("nucl" or "prot")
- ```out``` - base output name for database file


#### ii) On the Andromeda server, download Trembl database from [UniProt/Trembl](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz).

I created a script on Andromeda to download the updated Trembl database, unzip it, and make it into a blast database. I followed the general instructions on how to create code to download database [here](https://www.biostars.org/p/354449/).  Current databases found [here](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/).

```
pwd /data/putnamlab/shared/sbatch_executables/download_trembl_database.sh
```

Full script:

```
#!/bin/bash
#SBATCH --job-name="ref"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu # CHANGE EMAIL
#SBATCH -D /data/putnamlab/shared/databases/trembl_db

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b

cd databases/trembl_db

echo "Making trembl database" $date
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gunzip uniprot_trembl.fasta.gz

makeblastdb -in uniprot_trembl.fasta -parse_seqids -dbtype prot -out trembl_20211022
echo "STOP" $(date)

```

**It is normal to get multiple files per blast database. That is how makeblastdb is supposed to work. Just make sure files for a database stay together in the same directory and you use the "basename" for the database (a suggestion: name your database some thing other than your input file name) when you run your searches.**

**Options for makeblastdb which is an application that produces BLAST databases from FASTA files, other options [here](http://nebc.nerc.ac.uk/nebc_website_frozen/nebc.nerc.ac.uk/bioinformatics/documentation/blast+/user_manual.pdf)**

- ```makeblastdb``` - application produces BLAST databases from FASTA files
- ```in``` - database that will be used to search against
- ```parse_sequids``` - parse the seq-id(s) in the FASTA input provided
- ```dbtype``` - molecule type of target db ("nucl" or "prot")
- ```out``` - base output name for database file


#### iii) On the Andromeda server, download updated [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz) nr.gz database and unzip. This is updated daily on NCBI's website, so always re-download before you run your scripts for most updated version.

**If you are in the Putnam Lab and on Andromeda you can use the /data/shared/ncbi-nr/nr database that is updated by Kevin Bryan whenever you need to use contact him at bryank@uri.edu**  

**If you are not in the Putnam Lab, follow instructions [here](https://danielbruzzese.wordpress.com/2018/12/08/how-to-update-or-install-your-local-ncbi-blast-database-in-a-unix-shell-using-update_blastdb-pl/) to download the NCBI database to your server**


### Step 4: Align query protein sequences against databases

Now that the reference databases have been properly generated, the sequences of interest can be aligned against them.

#### 1) BLAST the protein sequences against Swiss-Prot

```
pwd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/Swissprot/swissprot_blast.sh

```

Full script:

```
#!/bin/bash
#SBATCH --job-name="swissprot-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="swissprot_blastp_out_error"
#SBATCH --output="swissprot_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against swissprot database" $(date)

blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/REFS/Pverr/Pver_proteins_names_v1.0.faa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out PverGeneModels_vs_sprot_1e-5_max5.out

echo "STOP" $(date)

```

Blastp: align protein query sequences against protein reference database

**Options**

- ```db``` - path to nr database file
- ```query``` - path to query fasta file
- ```out``` - base output name
- ```outfmt``` - output format
    - 6 = .out format
- ```evalue``` - maximum expected value to report an alignment
    - 1e-05 is typical cutoff for sequence alignments
- ```max_target_seqs``` maximum top sequences
    - Set at 5 to report multiple sequence for each gene
- ```num_threads``` maximum top sequences
    - Use <integer> CPU cores on a multicore system, if they are available


#### i) Get the best hit for each Gene Model (protein) Swiss-Prot

```
#Sort by 1. query name, 2. bitscore, 3. evalue, 4. protein identity, and extract the best line for each query (bitscore more important than evalue, evalue more important than nucleotide identity).

cat PverGeneModels_vs_sprot_1e-5_max5.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PverGeneModels_vs_sprot_1e-5_besthit.out

wc -l PverGeneModels_vs_sprot_1e-5_besthit.out #19,540

```

#### View output file

This is an example of the output .out file:

qaccver | saccver | pident | length | mismatch | gapopen | qstart | qend | sstart | send | evalue | bitscore | qlen |
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Pver_g1.t2 | sp_P25291_GP2_CANLF  | 38.346 | 133 | 65 | 3 | 175 | 306 |  36 | 152  |  2.58e-19 | 91.7 |  317 |

**Column names**

- qaccver - query seq-id
- saccver - subject seq-id (from BLAST databases)
- pident - % identical matches
- length - alignment length
- mismatch - number of mismatches
- gapopen - number of gap openings
- qstart - start of alignment in query
- qend - end of alignment in query
- sstart - start of alignment in subject
- send - end of alignment in subject
- evalue - number of expected hits of similar quality that could be found by chance (similar to a pvalue). The smaller the evalue, the better the match!
- bitscore - required size of a sequence database in which the current match could be found just by chance. The higher the bitscore, the better the sequence similarity!
- qlen - query sequence length


#### ii) Select the gene model proteins without hits in Swiss-Prot

```
#first use awk to print a list of all the Gene Model names from besthits.out

awk '{print $1}' PverGeneModels_vs_sprot_1e-5_besthit.out > list_of_Pvergenemodelproteins_sprot.txt

#then exclude these Gene Model names from your original fasta/.faa/protein file
#needed to load the module that has the script with the -exclude command in it

#first, loaded the newest module for kentUtils/416-foss-2020b

module load kentUtils/416-foss-2020b

#second use module show command to see paths to certain scripts and softwares in the module

module show kentUtils/416-foss-2020b

-------------------------------------------------------------------
/opt/modules/all/kentUtils/416-foss-2020b:

module-whatis     Description: LiftOver, Blat and other utilities
module-whatis     Homepage: https://hgdownload.soe.ucsc.edu/
module-whatis     URL: https://hgdownload.soe.ucsc.edu/
conflict     kentUtils
prepend-path     CMAKE_PREFIX_PATH /opt/software/kentUtils/416-foss-2020b
prepend-path     PATH /opt/software/kentUtils/416-foss-2020b/bin
setenv         EBROOTKENTUTILS /opt/software/kentUtils/416-foss-2020b
setenv         EBVERSIONKENTUTILS 416
setenv         EBDEVELKENTUTILS /opt/software/kentUtils/416-foss-2020b/easybuild/kentUtils-416-foss-2020b-easybuild-devel
-------------------------------------------------------------------

#I selected to the prepend-path /opt/software/kentUtils/416-foss-2020b/bin to see if it took me to the 'faSomeRecords' script which it did

/opt/software/kentUtils/416-foss-2020b/bin/faSomeRecords

#I then ran the -exclude command to exclude the blasted Gene Models from the .faa file

/opt/software/kentUtils/416-foss-2020b/bin/faSomeRecords -exclude Pver_proteins_names_v1.0.faa list_of_Pvergenemodelproteins_sprot.txt Pver_proteins_names_v1.0.faa.prot4trembl

#check the number of Gene Models

wc -l Pver_proteins_names_v1.0.faa.prot4trembl #15,798

#using this file to blast against trembl
```

#### You will also need a .xml file in order to use Blast2Go later on, so use another script to also download the .xml file

```
pwd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/Swissprot/swissprot_blast_xml.sh

```

Full script:

```
#!/bin/bash
#SBATCH --job-name="swissprot-blastp-protein-xml"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="xml_blastp_out_error"
#SBATCH --output="xml_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against swissprot database with xml format out" $(date)
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/REFS/Pverr/Pver_proteins_names_v1.0.faa -evalue 1e-5 -outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out PverGeneModels_maxhit.xml

echo "STOP" $(date)

```


#### iii) Secure-copy output files to local computer

```
# After doing the Swissprot top hits to .xml file instead, .xml files take up a lot of room and cannot always be uploaded to GitHub so I save them to my computer
# From a new terminal window (ie not Andromeda or remote server)

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/Swissprot/PverGeneModels_maxhit.xml /Users/Danielle/Documents/URI/XML_files/Pocillopora_verrucosa/Swissprot

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/Swissprot/PverGeneModels_vs_sprot_1e-5_besthit.out /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Functional_Annotation/Swissprot


```


#### 2) BLAST the remaining protein sequences against Trembl

```
pwd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/Trembl/trembl_blastp.sh

```

Full Script:

```
#!/bin/bash
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="trembl_blastp_out_error"
#SBATCH --output="trembl_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against trembl database" $(date)
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/trembl_db/trembl_20211022 -query Pver_proteins_names_v1.0.faa.prot4trembl -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out PverGeneModels_vs_trembl_1e-5_max5.out

echo "STOP" $(date)

```

```
Submitted batch job 94823
```

#### i) Get the best hit for each Gene Model (protein) Trembl

```
#Sort by 1. query name, 2. bitscore, 3. evalue, 4. protein identity, and extract the best line for each query (bitscore more important than evalue, evalue more important than nucleotide identity).

cat PverGeneModels_vs_trembl_1e-5_max5.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PverGeneModels_vs_trembl_1e-5_besthit.out

wc -l PverGeneModels_vs_trembl_1e-5_besthit.out #7095


```

#### ii) Select the gene model proteins without hits in Trembl

```
#first use awk to print a list of all the Gene Model names from besthits.out

awk '{print $1}' PverGeneModels_vs_trembl_1e-5_besthit.out > list_of_Pvergenemodelproteins_trembl.txt

#load the newest module for kentUtils/416-foss-2020b

module load kentUtils/416-foss-2020b

#then exclude these Gene Model names from your original fasta/.faa/protein file

/opt/software/kentUtils/416-foss-2020b/bin/faSomeRecords -exclude Pver_proteins_names_v1.0.faa.prot4trembl list_of_Pvergenemodelproteins_trembl.txt Pver_proteins_names_v1.0.faa.prot4nr


#check the number of Gene Models

Pver_proteins_names_v1.0.faa.prot4nr #1,608

#using this file to blast against nr database
```

#### You will also need a .xml file in order to use Blast2Go later on, so use another script to also download the .xml file

```
pwd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/Trembl/xml_blastp.sh

```

Full script:

```
#!/bin/bash
#SBATCH --job-name="xml-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="xml_blastp_out_error"
#SBATCH --output="xml_blastp_out"
#SBATCH --exclusive
#SBATCH -c 36

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against trembl database for xml" $(date)
blastp -max_target_seqs 5 -num_threads $SLURM_CPUS_ON_NODE -db /data/putnamlab/shared/databases/trembl_db/trembl_20211022 -query Pver_proteins_names_v1.0.faa.prot4trembl -evalue 1e-5 -outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out Pver_protein_blastp_trembl.xml

echo "STOP" $(date)

```

#### iii) Secure-copy output files to local computer

```
# .xml files take up a lot of room and cannot always be uploaded to GitHub so I save them to my computer
# From a new terminal window (ie not Andromeda or remote server)

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/Trembl/Pver_protein_blastp_trembl.xml /Users/Danielle/Documents/URI/XML_files/Pocillopora_verrucosa/Trembl

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/Trembl/PverGeneModels_vs_trembl_1e-5_besthit.out /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Functional_Annotation/Trembl

```


#### 3) BLAST the remaining protein sequences against nr

```
pwd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/NCBI/ncbi_blastp.sh

```

Full Script:


```
#!/bin/bash
#SBATCH --job-name="ncbi-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="ncbi_blastp_out_error"
#SBATCH --output="ncbi_blastp_out"
#SBATCH --exclusive
#SBATCH -c 36

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against ncbi database" $(date)
blastp -max_target_seqs 5 -num_threads $SLURM_CPUS_ON_NODE -db /data/shared/ncbi-nr/nr -query Pver_proteins_names_v1.0.faa.prot4nr -evalue 1e-5 -outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out PverGeneModels_ncbi.xml

echo "STOP" $(date)

```

#### i) Get the best hit for each Gene Model (protein) NCBI

```
#Sort by 1. query name, 2. bitscore, 3. evalue, 4. protein identity, and extract the best line for each query (bitscore more important than evalue, evalue more important than nucleotide identity).

cat PverGeneModels_ncbi_max5hits.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PverGeneModels_vs_nr_1e-5_besthit.out


wc -l PverGeneModels_vs_nr_1e-5_besthit.out #237


```

#### You will also need a .xml file in order to use Blast2Go later on, so use another script to also download the .xml file

```
pwd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/NCBI/ncbi_blastp.sh

```

Full script:

```
#!/bin/bash
#SBATCH --job-name="ncbi-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="ncbi_blastp_out_error"
#SBATCH --output="ncbi_blastp_out"
#SBATCH --exclusive
#SBATCH -c 36

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against ncbi database" $(date)
blastp -max_target_seqs 5 -num_threads $SLURM_CPUS_ON_NODE -db /data/shared/ncbi-nr/nr -query Pver_proteins_names_v1.0.faa.prot4nr -evalue 1e-5 -outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out PverGeneModels_ncbi.xml

echo "STOP" $(date)

```


#### ii) Secure-copy output files to local computer

```
# From a new terminal window (ie not Andromeda or remote server)

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/NCBI/PverGeneModels_ncbi.xml /Users/Danielle/Documents/URI/XML_files/Pocillopora_verrucosa/NCBI

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/NCBI/PverGeneModels_vs_nr_1e-5_besthit.out /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Functional_Annotation/NCBI

```


### Step 5: Assign gene ontology terms to sequences

After blastp steps are complete, analysis can move to assigning gene ontology (GO) terms to sequences.

The [Gene Ontology](http://geneontology.org) is an extensive consortium that aims to understand gene function and provide functional annotation for organisms across the tree of life. It also maintains a controlled vocabulary of gene and gene attributes across species.

The Gene Ontology has a system to classify genes into terms. The terms are grouped into 3 categories:

1. Molecular funciton - molecular-level activities performed by gene products (ie RNA or proteins). Describes the activities rather than the entities (molecules or protein complexes)
    - Examples: Toll receptor binding, transcription regulator activity
2. Cellular component - cellular anatomy and/or locations where gene products perform function
    - Examples: mitochondrion, ribosome
3. Biological process - larger biological processes accomplished by interaction of multiple molecular activities
    - Examples: glucose transmembrane transport, DNA repair

With this information, genes can be described/annotated with multiple terms. These terms are called GO terms. Here is the basic structure of a GO term:

```
GOID: GO:0007165
Term: signal transduction
Ontology: BP
Definition: The cellular process in which a signal is conveyed to
    trigger a change in the activity or state of a cell. Signal
    transduction begins with reception of a signal (e.g. a ligand
    binding to a receptor or receptor activation by a stimulus such as
    light), or for signal transduction in the absence of ligand,
    signal-withdrawal or the activity of a constitutively active
    receptor. Signal transduction ends with regulation of a downstream
    cellular process, e.g. regulation of transcription or regulation of
    a metabolic process. Signal transduction covers signaling from
    receptors located on the surface of the cell and signaling via
    molecules located within the cell. For signaling between cells,
    signal transduction is restricted to events at and within the
    receiving cell.
Synonym: GO:0023033
Synonym: signaling pathway
Synonym: signalling pathway
Synonym: signaling cascade
Synonym: signalling cascade
Secondary: GO:0023033
```

**Elements**

- GOID - unique 7-digit identifier that is intended to be machine readable
- Term - human-readable name
- Ontology - molecular function (MF), cellular component (CC), biological process (BP)
- Definition - description of what the term means
- Synonym - alternate words closely related in meaning to term, relevant connections
- Secondary - ID created when 2 or more terms are identical in meaning and so are merged into a single term. Secondary ID preserves the excess GO terms

There is so much more information available with these terms, including relationships to other genes/gene products, graphical representation of related GO terms, and much more, but that is beyond the scope of this analysis.

There are many different methdods to assigning GO terms to genes/gene products of interest. This workflow will focus on InterProScan, BLAST2GO, and Uniprot. These tools were chosen because they all utilize different databases for annotation. Additionally, they will all be run in different places: InterProScan on HPC, BLAST2GO on local computer, and Uniprot online. Results can be compared across tools to assess how each performed.


### Step 6: Run InterProScan to obtain GO terms

[InterProScan](https://www.ebi.ac.uk/interpro/) is a website/software that provides functional analysis of proteins by using predictive models that look across protein databases to find homologies with query proteins. If a homology is identified in one of the databases, the information about that homology is used to assign the query protein a GO term.

Note: Many fasta files willl use an asterisk to denote a STOP codon. InterProScan does not accept special characters within the sequences, so I removed them prior to running the program using the code below:

```
cp Pver_proteins_names_v1.0.faa /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/InterProScan/
sed -i 's/*//g' Pver_proteins_names_v1.0.faa
```

InterProScan utilizes several member databases to enhance the chance of obtaining good protein info:

- Structural domains - Gene3D, Superfamily
- Functional annotation of families/domains - PIRSF, TIGR, Panther, Pfam, SMART, Prints, Hamap, ProSite
- Protein features - ProSite
- Prediction of conserved domains - ProDom

#### i) Run InterProScan on Andromeda

```
# On Andromeda

# Load module
module load InterProScan/5.52-86.0-foss-2021a
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i Pver_proteins_names_v1.0.faa -b ./Pver.interpro.20210927  -iprlookup -goterms -pa
interproscan.sh -mode convert -f GFF3 -i ./Pver.interpro.20210927.xml -b ././Pver.interpro.20210927

```

**Options**

- ```version``` - display version number of IPS and java
- ```f``` - output format
- ```i``` - path to input file
- ```b``` - base output name
- ```iprlookup``` - provides mapping from matched member databases to InterPro entries
- ```goterms``` - provides mapping to Gene Ontology
- ```pa``` - provides mapping from matches to pathway info, which is based on matched mantually curated InterPro enteries
- ```mode    ``` - convert - change file format

#### ii) View output file

Example of part of IPS output file:

seqid | source | type | start | end | score | strand | phase | attributes
--- | --- | --- | --- | --- | --- | --- | --- | --- |
Pver_evm.model.Segkk4293_pilon.5     | Pver | protein_match | 158 | 257 | 1.1E-19 | + | . | date=26-09-2020;Target=Pver_evm.model.Segkk4293_pilon.5 158 257;Ontology_term="GO:0006811","GO:0016021";ID=match$3_158_257;signature_desc=Neurotransmitter-gated ion-channel transmembrane region;Name=PF02932;status=T;Dbxref="InterPro:IPR006029"

**Column names**

- seqid - query sequence  id
- source - algorithm or software that generated feature
- type - type of feature
- start - start coordinates of feature
- end - end coordinates of feature
- score - similar to evalue or pvalue
- strand - + for positive strand, - for negative strand
- phase - ?
- attributes - list of feature attributes in format tag=value. There can be multiple tag=value attributes
    - Date
    - Target - target sequence for analysis and start/stop coordinates
    - Ontology - link between feature and ontology databses, like Gene Ontology (ie GO terms)
    - ID - feature ID, required to be unique for features that have 'children' (gene or mRNA)
    - Signature desc
    - Name - feature name, does not have to be unique
    - Status
    - Dbxref - link between feature and other databases, like InterPro

I noticed some inconsistencies with my output GFF3 output table and Jill's GFF3 output table.

On my output table I noticed there were a lot of extra pathway annotations that occurred after the Dbxref="InterPro:IPR000504" section of the output. Compared to Jill's that did not have this. While their is nothing wrong with this, it takes up a lot of space and makes the file much longer.

My Output GFF3 Format:

![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/interproscan.gff3.example.db.jpg)

Jill's Output GFF3 Format:

![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/interproscan.gff3.output.ja.jpg)

Since this information was not needed for downstream analysis, I removed it using the code below:

```
sed -e 's#^\(.*Dbxref="InterPro:[A-Z0-9]*"\).*#\1#' Pver.interpro.20210927.gff3 > Pver.interpro.20210927-smaller.gff3

```

```
I then checked to make sure nothing was deleted by this change:

1. looked at number of gene names

zgrep -c "^>" Pver.interpro.20210927.gff3

455630

zgrep -c "^>" Pver.interpro.20210927-smaller.gff3

455630

2. looked at number of Dbxref= references

grep -c 'Dbxref=' Pver.interpro.20210927.gff3

257193

grep -c 'Dbxref=' Pver.interpro.20210927-smaller.gff3

257193

3. looked at number of Ontology_term= references

grep -c 'Ontology_term=' Pver.interpro.20210927.gff3

127451

grep -c 'Ontology_term=' Pver.interpro.20210927-smaller.gff3

127451

4. looked at number of status=T references

grep -c 'status=T' Pver.interpro.20210927.gff3

429951

grep -c 'status=T' Pver.interpro.20210927-smaller.gff3

429951

```


**Full Andromeda Script:**
**Took ~16 hours to complete**

```
Pver_InterProScan.sh:

#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="interproscan_out_error"
#SBATCH --output="interproscan_out"
#SBATCH --exclusive

cd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/InterProScan/

echo "START $(date)"

# Load module
module load InterProScan/5.52-86.0-foss-2021a
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh --cpu $SLURM_CPUS_ON_NODE ...
interproscan.sh -version
interproscan.sh -f XML -i Pver_proteins_names_v1.0.faa -b Pver.interpro.20210927  -iprlookup -goterms -pa
interproscan.sh -mode convert -f GFF3 -i Pver.interpro.20210927.xml -b Pver.interpro.20210927

# -i is the input data
# -b is the output file base
# -f is formats
# -iprlookup enables mapping
# -goterms is GO Term
# -pa is pathway mapping
# -version displays version number

echo "DONE $(date)"

```
```
Script error:
ModuleCmd_Load.c(213):ERROR:105: Unable to locate a modulefile for 'InterProScan/5.46-81.0-foss-2019b'
openjdk version "11.0.2" 2019-01-15
OpenJDK Runtime Environment 18.9 (build 11.0.2+9)
OpenJDK 64-Bit Server VM 18.9 (build 11.0.2+9, mixed mode)
/var/spool/slurmd/job1931747/slurm_script: line 21: interproscan.sh: command not found
/var/spool/slurmd/job1931747/slurm_script: line 22: interproscan.sh: command not found
/var/spool/slurmd/job1931747/slurm_script: line 23: interproscan.sh: command not found

Error searching for module:

ModuleCmd_Load.c(213):ERROR:105: Unable to locate a modulefile for 'InterProScan/5.52-86.0-foss-2019b'
```

```
#Kevin had to download the background packages on bluewaves and suggested to use andromeda for faster running times, ran on andromeda and it worked

Submitted batch job 88966

#noticed in my output script that it mentioned a new version/module of InterProtScan available, had Kevin Bryan update this to InterProScan/5.52-86.0-foss-2021a

#also Kevin suggested adding these two additions to the above code to speed up the process:

#SBATCH --exclusive
interproscan.sh --cpu $SLURM_CPUS_ON_NODE ...

#added them above and submitted

Submitted batch job 89016
```

#### ii) Secure-copy output file to local computer

```
# From a new terminal window (ie not Andromeda or remote server)

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/InterProScan/Pver.interpro.20210927-smaller.gff3 /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Functional_Annotation/InterProScan/

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/InterProScan/Pver.interpro.20210927.xml /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Functional_Annotation/InterProScan/

```

### Step 7:  Run BLAST2GO to obtain GO terms

[BLAST2GO](https://www.blast2go.com) (B2G) is a bioinformatics tools for functional annotation and analysis of gene or protein sequences. It was originally developed to provide a user-friendly interface for GO annotation and now hosts many different functional annotation tools. The B2G software makes it possible to generate annotation without requiring writing any code. While B2G can do an extensive array of functions, this analysis primarily utilizes the GO mapping and annotation functions.

### i) Download BLAST2GO to personal computer and activate the Basic subscription plan.

The B2G application can be downloaded [here](https://www.blast2go.com/blast2go-pro/download-b2g). B2G is available for Mac, Windows, and Linux systems. 2GB of RAM is recommended. Additionally, Internet connection is required to use most application features.

Register for B2G Basic [here](https://www.blast2go.com/b2g-register-basic). B2G Basic is free and includes the necessary features for this analysis. Registering will generate an activation key, which be put into the B2G software. You must be a part of some research institution to obtain B2G Basic.

### ii) Load the XML files generated from DIAMOND BLAST, INTERPROSCAN, SWISSPROT, and TREMBL separately.

While B2G has the ability to run BLAST, this analysis prefers to use results from DIAMOND because of its high performance capability and senstivity. The XML file generated from DIAMOND BLAST will be used here.

To load the file, go to File<Load<Load Blast results<Load Blast XML (Legacy)

Once the file is loaded, a table loads with info about those BLAST results (nr, Tags, SeqName, Description, Length, Hits, e-Value, and sim mean). All of the cells should be orange with Tags that say BLASTED. This indicates that these sequences have only been blasted, not mapped or annotated.

### iii) Map GO terms

Mapping is the process of retrieving GO terms associated with the Description obtained by the DIAMOND BLAST search. Several mapping {steps} occur:

- BLAST result accessions are used to retrieve gene names from NCBI and those gene names are used to search the GO database.
- BLAST result accessions are used to retrieve protein information (with GO terms already annotated) through UniProt, which uses the databases SD, UniProt, Swiss-Prot, TrEMBL, RefSeq, GenPept and PDB.
- BLAST result accessions are directly searched in GO database.

To map results, select the mapping icon (white circle with green circle inside) at the top of the screen. Then select Run Mapping. A box will open up; don't change anything, click run. Depending on the number of BLAST results, mapping could take hours to days. B2G will continue running if the computer is in sleep mode. Mapping status can be checked under the Progress tab in the lower left box. If mapping retrieved a GO term that may be related to a certain sequence, that sequence row will turn light green.

### iv) Annotate GO terms

Now that GO terms have been retrieved, the annotation process will select GO terms from the GO pool obtained through mapping and assign them to query sequences. Several annotation steps occur:

- For all found GO terms, an annotation rule (AR) is applied. The rule seeks to find the most specific annotations with a certain level of reliability and can be adjusted in the settings prior to running annotation.
- If a candidate GO term is found, an annotation score is calculated that weights the highest hit similarity of that candidate GO term by its evidence code (from mapping step). The candidate GO term will only be assigned to a sequence if it passes a certain threshold based on calculations above. For more info about the math behind annotation, go [here](http://docs.blast2go.com/user-manual/gene-ontology-annotation/).

To annotate, select the annot icon (white circle with blue circle inside) at the top of the screen. Then select Run Annotation. A box will open up; don't change anything unless thresholds need to be adjusted. If no changes are necessary, click next through the boxes until the final one and click run. Depending on the mapping results, annotating could take hours to days. B2G will continue running if the computer is in sleep mode. Annotation status can be checked under the Progress tab in the lower left box. If a GO term has been successfully assigned to a sequence, that sequence row will turn blue.

### v) Export annotated sequences and info

To export the file with B2G annotated genes, go to File<Export<Export as Table. A box will open up,  and save as a text file. Re-save as csv as needed.

### vi) View output file

seqName | top_hit | length | evalue | simMean | GO.ID | GO_names
--- | --- | --- | --- | --- | --- | --- |
Pver_evm.model.Segkk0_pilon.11 | EDO40119.1 | 1779 | 1.9e-199 | 68.38 | P:GO:0007166; P:GO:0007275; C:GO:0016020; C:GO:0120025 | P:cell surface receptor signaling pathway; P:multicellular organism development; C:membrane; C:plasma membrane bounded cell projection

**Column names**

- seqName - query sequence
- top_hit - top match from DIAMOND BLAST
- length - alignment length
- evalue - number of expected hits of similar quality that could be found by chance (similar to a pvalue). The smaller the evalue, the better the match!
- simMean - mean similarity between query and top match
- GO.ID - gene ontology mapping results, GO terms and evidence codes
- GO_names - annotated GO terms



### Step 9: Merge all information for full annotation

Follow instructions/code [here](https://github.com/hputnam/Becker_E5/blob/master/Functional_Annotation/Scripts/pver_annot_compile.Rmd) to compile all annotation information in R.s
