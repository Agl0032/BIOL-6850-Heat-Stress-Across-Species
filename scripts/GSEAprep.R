thisDir <- dirname(parent.frame(2)$ofile)
thisDir <- "C:/Users/schne/OneDrive - Auburn University/college/10 spring 2023/biol 6850/HeatStressConsistency/scripts"
setwd(thisDir)
setwd("../Results/ThaEle")

phenodata_te <- read.csv("PHENO_DATA.csv", row.names=1)

####  Preparing DGE result for Functional Enrichment Analysis ###
      # Functional Genomics BIOL: 6850
      # This file uses the dds object produced by the DGEseq2, and the statistical results file that was written.
  
############ Preparing Data for GSEA and Cytoscape.  #############

### Merge 'gene names' with DGE results by Gene Model (gene_id)

# ## Import Annotation file with results from Blast to databases
# Anno <- read.csv("Daphnia_pulex.annotations_Name.csv", stringsAsFactors = FALSE, na.strings=c(""))
# summary(Anno)
# dim(Anno)
  
## Import the DGE results file make sure the gene model name is 'gene_id' to match annotation file
DGEresults_te <- read.csv("DGESeq_results.csv", stringsAsFactors = FALSE)
names(DGEresults_te)[1] <- "gene_id"
summary(DGEresults_te)
dim(DGEresults_te)

## Merge anno with DGE results. Merging by gene_id (the first column). 
# The all.y=FALSE will result in only the rows that are in common between the two files
#DGE_Anno <- merge(Anno,DGEresults,by="gene_id",all.y = FALSE)
#dim(DGE_Anno)
#summary(DGE_Anno)

DGEresults_te$Name <- gsub("^[^-]+-", "", DGEresults_te$gene_id)

#write.csv(as.data.frame(DGEresults_te), file="GSEA/DGE_results_GeneName.csv", row.names=FALSE)  



############################# Make ranked list for GSEA ####################
## Here we are calculating a rank for each gene, and adding that column
## The rank is based on the log2fold change (how up or down regulated they are) and the pvalue (significance) 
DGErank_te <-  within(DGEresults_te, rank <- sign(log2FoldChange) * -log10(pvalue))
DGErank_te
 
#subset the results so only Gene Name and rank
DGErank_te$Name <- gsub("^[^-]+-", "", DGErank_te$gene_id)
DGErank_te = subset(DGErank_te, select = c(Name,rank) )
DGErank_te

#subset the results so only rows with a Name and rank
DGErank_withName_te <- na.omit(DGErank_te)
DGErank_withName_te
dim(DGErank_withName_te)

#write.csv(as.data.frame(DGErank_withName), file="DGErankName.csv", row.names=FALSE) 
## for use in GSEA preranked analysis, the file must be tab delimitied and end in .rnk
write.table(DGErank_withName_te, file="GSEA/DGErankName.rnk", sep = "\t", row.names=FALSE, quote=FALSE)  

####  We also need the normalized count data. Here we are going back to the dds object
nt_te <- normTransform(dds_te) # defaults to log2(x+1)
head(assay(nt_te))
# compare to original count data
head(countdata_te)

# make it a new dataframe
NormTransExp_te <- assay(nt_te)
summary(NormTransExp_te)
head(NormTransExp_te)

gene_id_te <- gsub("^[^-]+-", "", rownames(NormTransExp_te))
#Description <- vector(mode = "character", length=length(gene_id_te))
Description_te <- rep("na", length=length(gene_id_te))
NormTransExpIDs_te  <- cbind(gene_id_te, Description_te, NormTransExp_te)
colnames(NormTransExpIDs_te)[1] <- "Name"
head(NormTransExpIDs_te)

# #merge with name
# Name_NormTransExp <- merge(Anno, NormTransExpIDs, by="gene_id", all.y = FALSE)
# dim(Name_NormTransExp_te)
# summary(Name_NormTransExp_te)

write.table(as.data.frame(NormTransExpIDs_te), file="GSEA/NormTransExpressionData.txt", row.names=FALSE, sep = "\t", quote = FALSE)  

########## Make class file

rownames(phenodata_te) == colnames(dds_te)

### Make .cls file for ThaEle
line1 <- paste(length(phenodata_te$treatment), length(unique(sort(phenodata_te$treatment))), "1")
line2 <- paste("#", "control", "heat")
line3 <- paste(phenodata_te[rownames(phenodata_te) %in% colnames(dds_te), "treatment"],
               sep = '',
               collapse = ' ')

cat(line1, "\n", line2, "\n", line3, file="GSEA/Treatments.cls")






############################### Repeat for SceUnd

setwd(thisDir)
setwd("../Results/SceUnd")

phenodata_su <- read.csv("PHENO_DATA.csv", row.names=1)

## Import the DGE results file make sure the gene model name is 'gene_id' to match annotation file
DGEresults_su <- read.csv("DGESeq_results.csv", stringsAsFactors = FALSE)
names(DGEresults_su)[1] <- "gene_id"
summary(DGEresults_su)
dim(DGEresults_su)

## Merge anno with DGE results. Merging by gene_id (the first column). 
# The all.y=FALSE will result in only the rows that are in common between the two files
#DGE_Anno <- merge(Anno,DGEresults,by="gene_id",all.y = FALSE)
#dim(DGE_Anno)
#summary(DGE_Anno)

DGEresults_su$Name <- gsub("^[^-]+-", "", DGEresults_su$gene_id)

#write.csv(as.data.frame(DGEresults_su), file="DGE_results_GeneName.csv", row.names=FALSE, quote=FALSE)  



############################# Make ranked list for GSEA ####################
## Here we are calculating a rank for each gene, and adding that column
## The rank is based on the log2fold change (how up or down regulated they are) and the pvalue (significance) 
DGErank_su <-  within(DGEresults_su, rank <- sign(log2FoldChange) * -log10(pvalue))
DGErank_su <- DGErank_su[DGErank_su$rank != Inf,]
DGErank_su
 
#subset the results so only Gene Name and rank
DGErank_su$Name <- gsub("^[^-]+-", "", DGErank_su$gene_id)
DGErank_su = subset(DGErank_su, select = c(Name,rank) )
DGErank_su

#subset the results so only rows with a Name and rank
DGErank_withName_su <- na.omit(DGErank_su)
DGErank_withName_su
dim(DGErank_withName_su)

#write.csv(as.data.frame(DGErank_withName), file="DGErankName.csv", row.names=FALSE) 
## for use in GSEA preranked analysis, the file must be tab delimitied and end in .rnk
write.table(DGErank_withName_su, file="GSEA/DGErankName.rnk", sep = "\t", row.names=FALSE, quote=FALSE)  

####  We also need the normalized count data. Here we are going back to the dds object
nt_su <- normTransform(dds_su) # defaults to log2(x+1)
head(assay(nt_su))
# compare to original count data
head(countdata_su)

# make it a new dataframe
NormTransExp_su <- assay(nt_su)
summary(NormTransExp_su)
head(NormTransExp_su)

gene_id_su <- gsub("^[^-]+-", "", rownames(NormTransExp_su))
#Description <- vector(mode = "character", length=length(gene_id_su))
Description_su <- rep("na", length=length(gene_id_su))
NormTransExpIDs_su  <- cbind(gene_id_su, Description_su, NormTransExp_su)
colnames(NormTransExpIDs_su)[1] <- "Name"
head(NormTransExpIDs_su)

# #merge with name
# Name_NormTransExp <- merge(Anno, NormTransExpIDs, by="gene_id", all.y = FALSE)
# dim(Name_NormTransExp_su)
# summary(Name_NormTransExp_su)

write.table(as.data.frame(NormTransExpIDs_su), file="GSEA/NormTransExpressionData.txt", row.names=FALSE, sep = "\t", quote = FALSE)  

########## Make class file

rownames(phenodata_su) == colnames(dds_su)

### Make .cls file for ThaEle
line1 <- paste(length(phenodata_su$treatment), length(unique(sort(phenodata_su$treatment))), "1")
line2 <- paste("#", "control", "heat")
line3 <- paste(phenodata_su[rownames(phenodata_su) %in% colnames(dds_su), "treatment"],
               sep = '',
               collapse = ' ')

cat(line1, "\n", line2, "\n", line3, file="GSEA/Treatments.cls")
