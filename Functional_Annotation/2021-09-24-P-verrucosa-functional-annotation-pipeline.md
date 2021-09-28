---
layout: post
title: P. verrucosa functional annotation pipeline 
Author: Danielle Becker-Polinski 
Last Updated: 2021/09/28
tags: [ Protocol, annotation, RNASeq, GO, KEGG ]
---

## Overview

All methods for this protocol were adapted following the workflows created by @echillie and @JillAshey GitHub's [here](https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/blob/main/1-BLAST-GO-KO/2020-10-08-M-capitata-functional-annotation-pipeline.md) and [here](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#b-navigate-to-uniprot-retrieveid-mapping-page-and-under-provide-your-identifiers-click-upload-your-own-file-and-upload-a-subsetted-tab-file).

# Table of Contents

1. [Obtain sequences of interest](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#step-1-obtain-sequences-of-interest)

    i) [Identify species to work with](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#i-identify-species-to-work-with)
    
    ii) [Download genomic files for species of interest](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#ii-download-genomic-files-for-species-of-interest)
    
2. [Identify homologous sequences](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#step-2-identify-homologous-sequences)

    i) [On a HPC server, download nr database from NCBI. Convert this database to be Diamond-readable](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#i-on-a-hpc-server-download-nr-database-from-ncbi-convert-this-database-to-be-diamond-readable)
    
    ii) [Align query protein sequences against database](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#ii-align-query-protein-sequences-against-database)
    
    iii) [Generate readable output files](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#iii-generate-readable-output-files) 
    
    iv) [View output file](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#iv-view-output-file)
    
    v) [Secure-copy output files to local computer](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#v-secure-copy-output-files-to-local-computer)

3. [Assign gene ontology terms to sequences](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#step-3-assign-gene-ontology-terms-to-sequences)

    i) [Run InterProScan to obtain GO terms](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#i-run-interproscan-to-obtain-go-terms)
    
    - a) [Run InterProScan on Andromeda](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#a-run-interproscan-on-Andromeda)
    - b) [View output file](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#b-view-output-file-copy-to-local-computer)
    - c) [Secure-copy output file to local computer](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#b-view-output-file-copy-to-local-computer)

    ii) [Run BLAST2GO to obtain GO terms](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#ii-run-blast2go-to-obtain-go-terms)

    - a) [Download BLAST2GO to personal computer and activate the Basic subscription plan](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#a-download-blast2go-to-personal-computer-and-activate-the-basic-subscription-plan)
    - b) [Load the XML files generated from DIAMOND BLAST](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#b-load-the-xml-files-generated-from-diamond-blast)
    - c) [Map GO terms](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#c-map-go-terms)
    - d) [Annotate GO terms](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#d-annotate-go-terms)
    - e) [Export annotated sequences and info](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#e-export-annotated-sequences-and-info)
    - f) [View output file](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#f-view-output-file)

    iii) [Run Uniprot to obtain GO terms](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#iii-run-uniprot-to-obtain-go-terms)
    
    - a) [Make a list of identifiers found by DIAMOND BLAST](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#a-make-a-list-of-identifiers-found-by-diamond-blast)
    - b) [Navigate to Uniprot Retrieve/ID mapping page and under Provide your identifiers, click 'upload your own file' and upload a subsetted tab file](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#b-navigate-to-uniprot-retrieveid-mapping-page-and-under-provide-your-identifiers-click-upload-your-own-file-and-upload-a-subsetted-tab-file)
    - c) [c) Under 'Select options', choose the place where the identifiers were generated (From) and what database to compare to (To)](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#c-under-select-options-choose-the-place-where-the-identifiers-were-generated-from-and-what-database-to-compare-to-to)
    - d) [Select appropiate columns to include in table](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#d-select-appropiate-columns-to-include-in-table)
    - e) [Save files of interest](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#e-save-files-of-interest)
    - f) [View output file](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#f-view-output-files)


4. [Merge all information for full annotation](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-24-P-verrucosa-functional-annotation-pipeline.md#step-4-merge-all-information-for-full-annotation)


This project aims to develop a functional genomic annotation workflow for non-model organisms (ie corals). The goal of functional annotation is to identify and tag genes in a refernce genome with known functions of homologous genes in other organisms. The following document is intended as a tutorial in understanding functional annotation in non-model organisms. The genomic information from the coral *Pocillopora verrucosa* is used in this workflow.

For this functional annotation workflow tutorial, you will need:

- Access to a high performance computing server (ie URI Andromeda) with the following programs:
**Make sure to check the modules and see if there are any new versions. If so, contact Kevin Bryant to update modules on Andromeda when needed**
    - DIAMOND v2.0.0-GCC-8.3.0
    - InterProScan v5.52-86.0-foss-2019b
    - Java v11.0.2
- Laptop with access to the Internet and the following programs installed: 
    - Blast2GO Basic 
    - R v4.0.2
    - RStudio v1.3.959


### Step 1: Obtain sequences of interest. 

In order to conduct functional annotation steps, protein and transcript sequences are needed. There are two main hubs where coral genomic information is stored: [Reef Genomics](http://reefgenomics.org) and [NCBI](https://www.ncbi.nlm.nih.gov). Other researchers store their genomic infomation on their own personal webpages. Genomic information must be downloaded from one of these databases in order to proceed. 

#### i) Identify species to work with. 

For this project, the coral species of interest is *Pocillopora verrucosa*. 
 
#### ii) Download genomic files for species of interest. 
 
##### [*Pocillopora verrucosa* Gene models (CDS) (mRNA) ](http://baumslab.org/research/data/)

##### [*Pocillopora verrucosa* Gene models (protein)](http://baumslab.org/research/data/)

Ready to start annotating!

### Step 2: Identify homologous sequences

Homology refers to the similarity of structure or genes in different taxa due to shared ancestry. {example}

Sequence homology is the homology between DNA, RNA, and protein sequences in terms of shared ancestry. Sequence homology is usually inferred by the similarity of nucleotide or amino acid sequences. Strong sequence similarity (or percent homology) provides evidence that two or more sequences are related through shared ancestry. 

Several software programs have the ability to compare sequences to assess homology, the most popular one being NCBI's [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (Basic Local Alignment Search Tool - what a great name). BLAST compares nucleotide or protein sequences of interest (called a query) to their sequence databases to find regions of similarity. If a nucleotide or protein sequence of interest significantly matches a sequence/sequences in the databases, BLAST will tag the sequence of interest with the information about the known sequence(s). For example, if a new transcript sequence is identified in a mouse, BLAST could be used to see if any other animals carry a similar sequence, and if they do, what the biological functions of the sequence are. 

In this analysis, the program [DIAMOND BLAST](http://www.diamondsearch.org/index.php) was used. Like regular BLAST, DIAMOND is a sequence aligner for nucleotide and protein sequences; unlike BLAST, it is optimized for a higher performance capability of large datasets at 100x-20,000x speed of BLAST. 

#### i) On a HPC server, download nr database from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz). Convert this database to be Diamond-readable.

**If you are in the Putnam Lab and on Andromeda**  

**I used this step for this functional annotation because using the below commands in bash took way too long**
This script, created by Erin Chille on August 6, 2020, downloads the most recent nr database in FASTA format from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz) and uses it to make a Diamond-formatted nr database.  This step was updated by Danielle Becker-Polinski on September 24th, 2021 because the scripts were not including the full CPUs to download and a couple other formatting errors. Go to the *sbatch_executables* subdirectory in the Putnam Lab *shared* folder and run the scripts, ```make_diamond_nr_db.sh```  and  ```make_diamond_nr_db.sh``` in this order:

```
$ sbatch download_nr_database.sh
Submitted batch job NNN
$ sbatch -d afterok:NNN make_diamond_nr_db.sh
```

**If you are not in the Putnam Lab**  
Download the nr database from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz). Then, use Diamond's ```makedb``` command to format the database in a Diamond-friendly format. You can also use the command ```dbinfo``` to find version information for the database.

The nr (non-redundant) database is a collection of non-identical protein sequences compiled by NCBI. It is updated on a daily basis.

```
# On Andromeda

# Download db from NCBI 
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

# Load Diamond module
module load DIAMOND/2.0.0-GCC-8.3.0 

diamond makedb --in nr.gz -d nr
diamond dbinfo -d nr.dmnd
```

**Options**

- ```makedb``` - create a Diamond binary database file 
- ```in``` - database that will be used to search against 
- ```d``` - base output name for database file 
- ```dbinfo``` - print information about database file

The file nr.dmnd is now Diamond-readable and can be used as a reference.

#### ii) Align query protein sequences against database 

Now that the reference database has been properly generated, the sequences of interest can be aligned against it.

```
# On Andromeda

# Load Diamond module 
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

#Run sequence alignment against the nr database
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/REFS/Pverr/Pver_mRNA_v1.0.fna -o Pver.annot.20210924 -f 100 -b 20 --more-sensitive -e 0.00001 -k1 --unal=1

```

Blastx: align translated DNA query sequences against protein reference database 

**Options**

- ```d``` - path to nr database file 
- ```q``` - path to query fasta file 
- ```o``` - base output name 
- ```f``` - output format
    - 100 = Diamond format
- ```b``` - Block size in billions of sequence letters to be processed at a time. Larger block sizes increase the use of memory and temporary disk space, but also improve performance. Set at 20. 20 is the highest recommended value. 
- ```more-sensitive```
- ```e``` - maximum expected value to report an alignment
    - 1e-05 is typical cutoff for sequence alignments
- ```k``` maximum top sequences
    - Set at 1 to only report top sequence for each gene 

#### iii) Generate readable output files

```
# On Andromeda

#Converting format to XML format for BLAST2GO
diamond view -a Pver.annot.20210924.daa -o Pver.annot.20210924.xml -f 5
diamond view -a Pver.annot.20210924.daa -o Pver.annot.20210924.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

```

View: generate formatted output files 

**Options**

- ```a``` - path to input DAA file
- ```o``` - base output name 
- ```f``` - output format 
    - 5 = XML output; 6 = TAB output 

The output files (XML and TAB) will both be used downstream in this workflow.

#### iv) View output file

This is an example of part of the DIAMOND output .tab file:

qseqid | sseqid | pident | length | mismatch | gapopen | qstart | qend | sstart | send | evalue | bitscore | qlen | slen
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Pver_evm.model.Segkk0_pilon.12 | XP_015753513.1 | 86.7 | 1103 | 15 | 1 | 1 | 2913 | 23 | 1125 | 0.0e+00 | 1807.7 | 2916 | 1125

**Column names**

- qseqid - query seq-id
- sseqid - subject seq-id (from BLAST databases)
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
- slen - subject sequence length

#### v) Secure-copy output files to local computer 

```
# From a new terminal window (ie not Andromeda or remote server)

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/Diamond/Pver.annot.20210924.xml /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RAnalysis/Genome/BLAST_GO_KO/

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/Diamond/Pver.annot.20210924.tab /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RAnalysis/Genome/BLAST_GO_KO/
```

DIAMOND BLAST results can now be used in further analyses. 

**Full Andromeda Script:**

```
Pver_annot_diamond.sh:

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="diamond_blastx_out_error"
#SBATCH --output="diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Pver annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/REFS/Pverr/Pver_mRNA_v1.0.fna -o Pver_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Pver_annot.daa -o Pver_annot.xml -f 5
diamond view -a Pver_annot.daa -o Pver_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)
```
```
Submitted batch job 1931742
```

### Step 3: Assign gene ontology terms to sequences

After DIAMOND BLAST is completed, analysis can move to assigning gene ontology (GO) terms to sequences. 

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

#### i) Run InterProScan to obtain GO terms

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

##### a) Run InterProScan on Andromeda

```
# On Andromeda 
 
# Load module
module load InterProScan/5.52-86.0-foss-2019b
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

##### b) View output file

Example of part of IPS output file: 

seqid | source | type | start | end | score | strand | phase | attributes 
--- | --- | --- | --- | --- | --- | --- | --- | --- |
Pver_evm.model.Segkk4293_pilon.5     | Pver | protein_match | 158 | 257 | 1.1E-19 | + | . | date=26-09-2020;Target=Acerv_evm.model.Segkk4293_pilon.5 158 257;Ontology_term="GO:0006811","GO:0016021";ID=match$3_158_257;signature_desc=Neurotransmitter-gated ion-channel transmembrane region;Name=PF02932;status=T;Dbxref="InterPro:IPR006029"

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
    

**Full Andromeda Script:**

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

cd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/InterProScan/

echo "START $(date)"

# Load module
module load InterProScan/5.52-86.0-foss-2019b
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i Pver_proteins_names_v1.0.faa -b ./Pver.interpro.20210927  -iprlookup -goterms -pa
interproscan.sh -mode convert -f GFF3 -i ./Pver.interpro.20210927.xml -b ././Pver.interpro.20210927

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
```

##### c) Secure-copy output file to local computer 

```
# From a new terminal window (ie not Andromeda or remote server)

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/InterProScan/Pver.annot.20210924.gff3 /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RAnalysis/Genome/BLAST_GO_KO/InterProScan/

```

#### ii) Run BLAST2GO to obtain GO terms

[BLAST2GO](https://www.blast2go.com) (B2G) is a bioinformatics tools for functional annotation and analysis of gene or protein sequences. It was originally developed to provide a user-friendly interface for GO annotation and now hosts many different functional annotation tools. The B2G software makes it possible to generate annotation without requiring writing any code. While B2G can do an extensive array of functions, this analysis primarily utilizes the GO mapping and annotation functions. 

##### a) Download BLAST2GO to personal computer and activate the Basic subscription plan. 

The B2G application can be downloaded [here](https://www.blast2go.com/blast2go-pro/download-b2g). B2G is available for Mac, Windows, and Linux systems. 2GB of RAM is recommended. Additionally, Internet connection is required to use most application features.

Register for B2G Basic [here](https://www.blast2go.com/b2g-register-basic). B2G Basic is free and includes the necessary features for this analysis. Registering will generate an activation key, which be put into the B2G software. You must be a part of some research institution to obtain B2G Basic. 

##### b) Load the XML files generated from DIAMOND BLAST. 

While B2G has the ability to run BLAST, this analysis prefers to use results from DIAMOND because of its high performance capability and senstivity. The XML file generated from DIAMOND BLAST will be used here.

To load the file, go to File<Load<Load Blast results<Load Blast XML (Legacy)

Once the file is loaded, a table loads with info about those BLAST results (nr, Tags, SeqName, Description, Length, Hits, e-Value, and sim mean). All of the cells should be orange with Tags that say BLASTED. This indicates that these sequences have only been blasted, not mapped or annotated. 

##### c) Map GO terms

Mapping is the process of retrieving GO terms associated with the Description obtained by the DIAMOND BLAST search. Several mapping {steps} occur:

- BLAST result accessions are used to retrieve gene names from NCBI and those gene names are used to search the GO database. 
- BLAST result accessions are used to retrieve protein information (with GO terms already annotated) through UniProt, which uses the databases SD, UniProt, Swiss-Prot, TrEMBL, RefSeq, GenPept and PDB.
- BLAST result accessions are directly searched in GO database. 

To map results, select the mapping icon (white circle with green circle inside) at the top of the screen. Then select Run Mapping. A box will open up; don't change anything, click run. Depending on the number of BLAST results, mapping could take hours to days. B2G will continue running if the computer is in sleep mode. Mapping status can be checked under the Progress tab in the lower left box. If mapping retrieved a GO term that may be related to a certain sequence, that sequence row will turn light green.

##### d) Annotate GO terms

Now that GO terms have been retrieved, the annotation process will select GO terms from the GO pool obtained through mapping and assign them to query sequences. Several annotation steps occur:

- For all found GO terms, an annotation rule (AR) is applied. The rule seeks to find the most specific annotations with a certain level of reliability and can be adjusted in the settings prior to running annotation. 
- If a candidate GO term is found, an annotation score is calculated that weights the highest hit similarity of that candidate GO term by its evidence code (from mapping step). The candidate GO term will only be assigned to a sequence if it passes a certain threshold based on calculations above. For more info about the math behind annotation, go [here](http://docs.blast2go.com/user-manual/gene-ontology-annotation/). 

To annotate, select the annot icon (white circle with blue circle inside) at the top of the screen. Then select Run Annotation. A box will open up; don't change anything unless thresholds need to be adjusted. If no changes are necessary, click next through the boxes until the final one and click run. Depending on the mapping results, annotating could take hours to days. B2G will continue running if the computer is in sleep mode. Annotation status can be checked under the Progress tab in the lower left box. If a GO term has been successfully assigned to a sequence, that sequence row will turn blue.

##### e) Export annotated sequences and info

To export the file with B2G annotated genes, go to File<Export<Export as Table. A box will open up,  and save as a text file. Re-save as csv as needed.

##### f) View output file

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


#### iii) Run Uniprot to obtain GO terms

[Uniprot](https://www.uniprot.org) (Universal Protein Resource) provides information for protein sequence and annotation data. It maintains several protein databases:

- UniProt Knowledgebase (UniProtKB) - central collection hub for functional protein info and annotation
- UniProt Reference Clusters (UniRef) - provides clustering from UniProtKB to ensure complete coverage of sequences while masking redundant sequences 
- UniProt Archive (UniParc) - database of publicly available protein sequences 

In this analysis, Uniprot uses BLAST input to search against its protein databases. 

##### a) Make a list of identifiers found by DIAMOND BLAST. 

Uniprot uses the .tab file generated from DIAMOND BLAST as input. The column 'sseqid' is the primary input to Uniprot, as it can use that BLAST input to find protein matches. Because UniProt is a website and lacks proper storage and RAM, it cannot handle an entire .tab DIAMOND BLAST file. To ensure that UniProt reads all inputs, subset files so that a file only has ~2000-3000 lines per file in Terminal. (there is probably a better/more efficient way to do this)

```
# Check how many lines in the full file
wc -l Pver_annot.20210924.tab
   29515 Pver_annot.20210924.tab

sed -n '1,2000p' Pver_annot.20210924.tab > Pver_annot.20210924_1.tab
sed -n '2001,4000p' Pver_annot.20210924.tab > Pver_annot.20210924_2.tab
sed -n '4001,6000p' Pver_annot.20210924.tab > Pver_annot.20210924_3.tab
.
.
.
.
.
.
.
.
.
.
.
.
.
sed -n '28001,29515p' Pver_annot.20210924.tab > Pver_annot.20210924.tab
```

##### b) Navigate to Uniprot [Retrieve/ID mapping page](https://www.uniprot.org/uploadlists/) and under Provide your identifiers, click 'upload your own file' and upload a subsetted tab file. 

UniProt will be able to process the smaller/subsetted files. However, each subsetted .tab file will need to be put in as input one at a time. 

##### c) Under 'Select options', choose the place where the identifiers were generated (From) and what database to compare to (To). 

In this analysis, the identifiers were generated from the EMBL/GenBank/DDBJ CDS option (aka NCBI/BLAST) and they will be compared to (UniProtKB). Once that is finished, hit submit. It may take a few minutes). If the mapping fails and the page displays 'Service Unavailable', go back and try again. 

##### d) Select appropiate columns to include in table. 

If mapping was successful, the screen should have something like ```X out of Y EMBL/GenBank/DDBJ CDS identifiers were successfully mapped to X UniProtKB IDs in the table below```. There will be a table below with the results. Select the 'Columns' tab right above the table. This will open a window to pick more column options so more information can be included in the table. Many of the column options don't apply to this analysis. 

Under 'Names & Taxonomy', Entry name, Gene names, Organism, and Protein Names should already be selected. Under 'Sequences', Length should already be selected. Under 'Miscellanous', the gold paper with a star and the blue paper symbols should already be selected. If they are not, select those columns. Under 'Gene Ontology', select Gene ontology (GO) and Gene ontology IDs. Under 'Genome Annotation', select KEGG. Once selections are completed, scroll back to the top of the page and hit Save in the upper-right hand corner. 

##### e) Save files of interest

To save the main table, click Download in the top left of the table. Select Download all, Uncompressed, and Tab-seperated as Format. Click Go.

To save the unmapped identifiers, click on the 'Click here to download the 6281 unmapped identifiers'. To save UniParc results (if any are found), click on the UniParc link under the main header of ```X out of Y EMBL/GenBank/DDBJ CDS identifiers were successfully mapped to X UniProtKB IDs in the table below```. To finish downloading UniParc results, click Download in the top left of the table. Select Download all, Uncompressed, and Tab-seperated as Format. Click Go.

To see UniProt analysis for all coral species, go [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/Uniprot/Uniprot.md) 

##### f) View output files 

This is an example of Uniprot output .tab file:

my_list | top_hit | uniprotkb_entry | status | protein_names | gene_names | organism | length | go_ids | gene_ontology | ko | kegg 
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
EDO36971.1 | A7SH22 | A7SH22_NEMVE | unreviewed | Predicted protein (Fragment) | v1g118173 | Nematostella vectensis (Starlet sea anemone) | 239 | GO:0004888; GO:0005230; GO:0005887; GO:0007165; GO:0007268; GO:0030594; GO:0034220; GO:0042391; GO:0043005; GO:0045202; GO:0050877 | integral component of plasma membrane [GO:0005887]; neuron projection [GO:0043005]; synapse [GO:0045202]; extracellular ligand-gated ion channel activity [GO:0005230]; neurotransmitter receptor activity [GO:0030594]; transmembrane signaling receptor activity [GO:0004888]; chemical synaptic transmission [GO:0007268]; ion transmembrane transport [GO:0034220]; nervous system process [GO:0050877]; regulation of membrane potential [GO:0042391]; signal transduction [GO:0007165] | K05181 | nve:5508437

**Column names**

- my_list - query id (generated from DIAMOND BLAST)
- top_hit - top match from Uniprot databases
- uniprotkb_entry - official Uniprot entry for top_hit
- status - entry status; indicates if entry has been manually annotated and reviewed (reviewed = SwissProt section of Uniprot, unreviewed = computer-annotated TrEMBL section)
- protein_names - list of all names for a particular protein (Predicted protein = unreviewed = omputer-annotated TrEMBL section)
- gene_names - name of the genes that code for particular protein sequence 
- organism - name of organism that is source of protein information
- length - alignment length
- go_ids - gene ontology mapping results, GO terms and evidence codes
- gene_ontology - annotated GO terms 
- ko - K number used to reconstruct biological pathways in KEGG mapper (outside scope of this analysis)
- kegg - gene identifier used to search biological pathways in KEGG mapper (outside scope of this analysis)

### Step 4: Merge all information for full annotation

Follow instructions/code [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/RAnalysis/acerv_annot_compile.R) to compile all annotation information in R. 
