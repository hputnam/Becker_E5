# Title: Pver Annotation Compilation
# Author: Jill Ashey
# Updated by: Danielle Becker-Polinski
# Date: 20211015


#This script takes the results of functional annotation services and combines the results. Nucleotide sequences from ReefGenomics (http://pver.reefgenomics.org/download/) were annotated using DIAMONDSEARCH BLASTX, resulting in 21606 hits. These hits were used as input into:
#  - Uniprot
#  - Blast2GO

#Additional annotation was provided by
#  - InterProScan



## Load libraries
library(tidyverse)
library(dplyr)



## Diamond BLAST

# Load DIAMOND BLAST results
pver_blast <- read_tsv("Functional_Annotation/Diamond/Pver_annot.tab", col_names = FALSE)
colnames(pver_blast) <- c("seqName", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",  "qlen", "slen")
head(pver_blast)
dim(pver_blast) #21606 x 14



## Uniprot

# Uniprot mapping occurred on 20211015.
#EMBL/GenBank/DDBJ CDS protein IDs generated from Diamond BLAST were mapped to UniProtKB database IDs.

#Because there were many Diamond BLAST hits for Pver, I had to break up the tab file into chunks that Uniprot could handle. So there are 11 Pver_Uniprot files to be read and then rbind() together
u1 <- read_tsv("Functional_Annotation/Uniprot/uniprot_tabaa.tab", col_names = TRUE)
colnames(u1) <- c("uniprotkb_entry", "entry_name", "status", "protein_names", "gene_names", "organism", "length", "gene_ontology", "go_ids", "kegg", "my_list")
head(u1)
dim(u1) # 495 x 12


u2 <- read_tsv("Functional_Annotation/Uniprot/uniprot_tabab.tab", col_names = TRUE)
colnames(u2) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u2)
dim(u2) # 230 x 12

u3 <- read_tsv("Functional_Annotation/Uniprot/uniprot_tabac.tab", col_names = TRUE)
colnames(u3) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u3)
dim(u3) # 216 x 12

u4 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/pver_Uniprot_4.tab", col_names = TRUE)
colnames(u4) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u4)
dim(u4) # 238 x 12

u5 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/pver_Uniprot_5.tab", col_names = TRUE)
colnames(u5) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u5)
dim(u5) # 218 x 12

u6 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/pver_Uniprot_6.tab", col_names = TRUE)
colnames(u6) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u6)
dim(u6) # 238 x 12

u7 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/pver_Uniprot_7.tab", col_names = TRUE)
colnames(u7) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u7)
dim(u7) # 232 x 12

u8 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/pver_Uniprot_8.tab", col_names = TRUE)
colnames(u8) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u8)
dim(u8) # 219 x 12

u9 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/pver_Uniprot_9.tab", col_names = TRUE)
colnames(u9) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u9)
dim(u9) # 106 x 12

u10 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/pver_Uniprot_10.tab", col_names = TRUE)
colnames(u10) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u10)
dim(u10) # 86 x 12

u11 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/pver_Uniprot_11.tab", col_names = TRUE)
colnames(u11) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u11)
dim(u11) # 107 x 12


#Compile the Uniprot files 
uniprot_results <- bind_rows(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11)
uniprot_results <- unique(uniprot_results)
head(uniprot_results)
dim(uniprot_results) # 1844 x 12
uniprot_results <- filter(uniprot_results, grepl("GO",go_ids)) #Select only gnes with GO terms
dim(uniprot_results) # 1095 x 12

# Generate a list of GO terms - UniProt
uniprot_GO <- select(uniprot_results, my_list, go_ids)
splitted <- strsplit(as.character(uniprot_GO$go_ids), ";") #split into multiple GO ids
uniprot_GO <- data.frame(v1 = rep.int(uniprot_GO$my_list, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
uniprot_GO <- unique(uniprot_GO)
colnames(uniprot_GO) <- c("gene_id", "GO.ID")
uniprot_GO$GO.ID <- gsub(" ", "", uniprot_GO$GO.ID)
nrow(uniprot_GO) # 2875 total GO terms from Uniprot
length(unique(uniprot_GO$GO.ID)) # 824 unique GO terms from Uniprot
length(unique(uniprot_GO$gene_id)) # 1095 unique gene ids from Uniprot 
# Not technically gene ids. Uniprot has no gene id info--actually ids from Uniprot. When I combine the Uniport and B2G files, the uniprot ids will then be associated with gene ids 




## Blast2GO

#Nucleotide CDS sequences were annotated using DIAMONDSEARCH BLASTX, resulting in 21606. These hits were used as input into Blast2GO to obtain GO terms using the 20211015 obo database.

B2G_results <- read.csv("Functional_Annotation/Blast2GO/pver_blast2go_table.csv")
B2G_results <- select(B2G_results, c("SeqName", "Description", "Length", "e.Value", "sim.mean", "GO.IDs", "GO.Names"))
colnames(B2G_results) <- c("seqName", "top_hit", "length", "eValue", "simMean", "GO.ID", "GO_names")
head(B2G_results)
dim(B2G_results) # 21606 x 7
B2G_results <- filter(B2G_results, grepl("GO",GO.ID)) #Genes with GO terms - 136
dim(B2G_results) # 2092 x 7

# Generate a list of GO terms - B2G
B2G_results_GO <- select(B2G_results, top_hit, GO.ID)
splitted <- strsplit(as.character(B2G_results_GO$GO.ID), ";") #split into multiple GO ids
B2G_results_GO <- data.frame(v1 = rep.int(B2G_results_GO$top_hit, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
B2G_results_GO <- unique(B2G_results_GO)
colnames(B2G_results_GO) <- c("gene_id", "GO.ID")
B2G_results_GO$GO.ID <- gsub("F:", "", B2G_results_GO$GO.ID)
B2G_results_GO$GO.ID <- gsub("C:", "", B2G_results_GO$GO.ID)
B2G_results_GO$GO.ID <- gsub("P:", "", B2G_results_GO$GO.ID)
B2G_results_GO$GO.ID <- gsub(" ", "", B2G_results_GO$GO.ID)
head(B2G_results_GO)
nrow(B2G_results_GO) # 4822 total GO terms from B2G
length(unique(B2G_results_GO$GO.ID)) # 1130 unique GO terms from B2G
length(unique(B2G_results_GO$gene_id)) # 1458 unique gene ids from B2G



# Find intersections and unique results for each method (Uniprot and B2G)

## GO
# Intersection between GO terms for B2G and Uniprot
BU_GO <- intersect(B2G_results_GO$GO.ID, uniprot_GO$GO.ID)
length(unique(BU_GO)) # 626 similar GO terms between B2G and Uniprot

# Difference in GO terms for B2G and Uniprot - B2G
BUunique_GO <- setdiff(B2G_results_GO$GO.ID, uniprot_GO$GO.ID)
length(unique(BUunique_GO)) # 504 GO terms unique to B2G

# Difference in GO terms for B2G and Uniprot - Uniprot
UBunique <- setdiff(uniprot_GO$GO.ID, B2G_results_GO$GO.ID)
length(unique(UBunique)) # 198 GO terms unique to Uniprot

## gene id
# Intersection between gene id for B2G and Uniprot
BU_gene <- intersect(B2G_results_GO$gene_id, uniprot_GO$gene_id)
length(unique(BU_gene)) # 1042 similar gene ids between B2G and Uniprot

# Difference in gene ids for B2G and Uniprot
BUunique_gene <- setdiff(B2G_results_GO$gene_id, uniprot_GO$gene_id)
length(unique(BUunique_gene)) # 416 gene ids unique to B2G

# Difference in gene ids for B2G and Uniprot
UBunique_gene <- setdiff(uniprot_GO$gene_id, B2G_results_GO$gene_id)
length(unique(UBunique_gene)) # 53 gene ids unique to uniprot



## Merge Annotations - uniprot + b2g
pver_annot <- left_join(pver_blast, B2G_results, by="seqName")
pver_annot <- select(pver_annot, seqName, top_hit.x, length.x, evalue, bitscore, simMean, GO.ID, GO_names)
colnames(pver_annot) <- c("gene_id", "top_hit", "length", "evalue", "bitscore", "simMean", "GO.ID", "GO_names")
uniprot_results_GO <- select(uniprot_results_GO, -top_hit)
uniprot_results_GO <- rename(uniprot_results_GO, "top_hit"="my_list")
pver_annot <- merge(pver_annot, uniprot_results_GO, by="top_hit", all.x = T)
pver_annot$GO.IDs <- paste(pver_annot$GO.ID, pver_annot$go_ids, sep=';') #generate new column with concatenated GO IDs
pver_annot$GO_terms <- paste(pver_annot$GO_names, pver_annot$gene_ontology, sep=';') #generate new column with concatenated GO IDs
pver_annot <- select(pver_annot, c("top_hit", "gene_id", "length.x", "evalue", "bitscore", "simMean", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "ko", "kegg", "GO.IDs", "GO_terms"))
pver_annot <- rename(pver_annot, "GO.ID"="GO.IDs")
names(pver_annot)
head(pver_annot)
tail(pver_annot)
dim(pver_annot) # 21606 x 16
write.csv(pver_annot, "Functional_Annotation/Final_Annotations/pver_FuncAnn_UniP_B2G.csv")




# IPS
IPS_GO <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOseq/pver_GOterms.csv", header=TRUE)
IPS_GO <- select(IPS_GO, -X)
colnames(IPS_GO)[1] <-"gene_id"
#IPS_GO$gene_id <- gsub(".m1", "", IPS_GO$gene_id)
#IPS_GO$gene_id <- gsub("model", "TU", IPS_GO$gene_id)
length(unique(IPS_GO$gene_id)) # 12828
IPS_GO <- select(IPS_GO, c("gene_id", "GO_term")) # taking out source and score tho I want to leave them in--just for comparison though
splitted <- strsplit(as.character(IPS_GO$GO), ",") #split into multiple GO ids
IPS_GO <- data.frame(v1 = rep.int(IPS_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their 
colnames(IPS_GO) <- c("gene_id", "GO.ID")
IPS_GO <- unique(IPS_GO)
head(IPS_GO)
nrow(IPS_GO) # 30085 total GO terms from IPS
length(unique(IPS_GO$GO.ID)) # 1945 unique GO terms from IPS
length(unique(IPS_GO$gene_id)) # 12828 unique gene ids from IPS

# From B2G and Uniprot
FuncAnn_UniP_B2G <- read.csv("Functional_Annotation/Final_Annotations/pver_FuncAnn_UniP_B2G.csv")
length(unique(FuncAnn_UniP_B2G$gene_id)) # 21606
FuncAnn_UniP_B2G_GO <- filter(FuncAnn_UniP_B2G, grepl("GO",GO.ID)) #Select only gnes with GO terms
length(unique(FuncAnn_UniP_B2G_GO$gene_id)) # 2172
FuncAnn_UniP_B2G_GO <- select(FuncAnn_UniP_B2G_GO, c("gene_id", "GO.ID")) # taking out source and score tho I want to leave them in--just for comparison though
splitted <- strsplit(as.character(FuncAnn_UniP_B2G_GO$GO.ID), ";") #split into multiple GO ids
FuncAnn_UniP_B2G_GO <- data.frame(v1 = rep.int(FuncAnn_UniP_B2G_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their 
colnames(FuncAnn_UniP_B2G_GO) <- c("gene_id", "GO.ID")
FuncAnn_UniP_B2G_GO$GO.ID <- gsub("F:", "", FuncAnn_UniP_B2G_GO$GO.ID)
FuncAnn_UniP_B2G_GO$GO.ID <- gsub("C:", "", FuncAnn_UniP_B2G_GO$GO.ID)
FuncAnn_UniP_B2G_GO$GO.ID <- gsub("P:", "", FuncAnn_UniP_B2G_GO$GO.ID)
FuncAnn_UniP_B2G_GO$GO.ID <- gsub(" ", "", FuncAnn_UniP_B2G_GO$GO.ID)
FuncAnn_UniP_B2G_GO <- unique(FuncAnn_UniP_B2G_GO)
head(FuncAnn_UniP_B2G_GO)
nrow(FuncAnn_UniP_B2G_GO) # 8390 total GO terms from B2G + Uniprot
length(unique(FuncAnn_UniP_B2G_GO$GO.ID)) # 1329 unique GO terms from B2G + Uniprot
length(unique(FuncAnn_UniP_B2G_GO$gene_id)) # 2172 unique gene ids from B2G + Uniprot

# Find intersections and unique results for each method (B2G + Uniprot and IPS)
## GO
# Intersection between GO terms for B2G + Uniprot and IPS
IF_GO <- intersect(IPS_GO$GO.ID, FuncAnn_UniP_B2G_GO$GO.ID)
length(unique(IF_GO)) # 747 similar GO terms between B2G + Uniprot and IPS

# Difference in GO terms for B2G + Uniprot and IPS - FuncAnn first
FIunique_GO <- setdiff(FuncAnn_UniP_B2G_GO$GO.ID, IPS_GO$GO.ID)
length(unique(FIunique_GO)) # 582 GO terms unique to B2G + Uniprot

# Difference in GO terms for B2G + Uniprot and IPS - IPS
IFunique_GO <- setdiff(IPS_GO$GO.ID, FuncAnn_UniP_B2G_GO$GO.ID)
length(unique(IFunique_GO)) # 1198 GO terms unique to Uniprot

## gene id
# Intersection between GO terms for B2G + Uniprot and IPS
IF_gene<- intersect(IPS_GO$gene_id, FuncAnn_UniP_B2G_GO$gene_id)
length(unique(IF_gene)) # 1140 similar gene ids between B2G + Uniprot and IPS

# Difference in gene ids for B2G + Uniprot and IPS - B2G + Uniprot
FIunique_gene <- setdiff(FuncAnn_UniP_B2G_GO$gene_id, IPS_GO$gene_id)
length(unique(FIunique_gene)) # 1032 gene ids unique to B2G + Uniprot

# Difference in gene ids for B2G + Uniprot and IPS - IPS
IFunique_gene <- setdiff(IPS_GO$gene_id, FuncAnn_UniP_B2G_GO$gene_id)
length(unique(IFunique_gene)) # 11688 gene ids unique to IPS


# Aggregate results 
agg <- aggregate(IPS_GO$GO.ID, list(IPS_GO$gene_id), paste, collapse = ",")
colnames(agg) <- c("gene_id", "GO.ID")

full_annot <- merge(agg, FuncAnn_UniP_B2G, by = "gene_id", all.x = T)
full_annot$GO.IDs <- paste(full_annot$GO.ID.x, full_annot$GO.ID.y, sep=';') #generate new column with concatenated GO IDs
full_annot <- select(full_annot, -c(X, GO.ID.x, GO.ID.y))
write.csv(full_annot, file = "Functional_Annotation/Final_Annotations/pver_fullAnnot.csv")




