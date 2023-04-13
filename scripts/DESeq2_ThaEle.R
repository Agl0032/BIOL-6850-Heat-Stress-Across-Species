####  DESeq2 Script for Differential Gene Expression Analysis in
      # Functional Genomics BIOL: 6850
### Resources and Citations:
# Love et al 2016 DESeq2 GenomeBiology
# https://bioconductor.riken.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

### Set working directory
thisDir <- dirname(parent.frame(2)$ofile) # automatic but only works on source
#thisDir <- "/path/to/script/folder" 
setwd(thisDir)
setwd("../Results/ThaEle")


#### Install the DESeq2 package if you have not already
## try http:// if https:// URLs are not supported
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

## Load the DESeq2 library
library(DESeq2)


##########   1.3 Input data   ##############
### Input the count data, the gene(/transcript) count matrix and labels
  ### How you import this will depend on what your final output was from the mapper/counter that you used.
    #  *** Important *** if your gene matrix doesn't have "gene_id" as the first column header then you need to add it for this import command to work.
    ## this works with output from PrepDE.py from Ballgown folder, if you add the gene_id to the first column.
countdata_te <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
dim(countdata_te)
head(countdata_te)


### Input the meta data or phenotype data
  # NOTE: Created by the merge-ID-treatment.R script.
  # Note: The PHENO_DATA file contains information on each sample, e.g., sex or
  # population. The exact way to import this depends on the format of the file.
coldata_te <- read.csv("PHENO_DATA.csv", header=TRUE, row.names=1)
dim(coldata_te)
head(coldata_te)

dim(countdata_te)

rownames(coldata_te)
colnames(countdata_te)

##  Make sure the individual names match between the count data and the metadata
#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata_te) %in% colnames(countdata_te))
coldata_te <- coldata_te[rownames(coldata_te) %in% colnames(countdata_te),]
countdata_te <- countdata_te[, rownames(coldata_te)]
all(rownames(coldata_te) == colnames(countdata_te))


## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
dds_te <- DESeqDataSetFromMatrix(countData = countdata_te, colData=coldata_te,  design = ~treatment)
dds_te

#####   Prefiltering    Manual - starting at  1.3.6
# Here we perform a minimal pre-filtering to remove rows that have less than 1 read mapped.
dds_te <- dds_te[ rowSums(counts(dds_te)) > 1, ]
  ## 23443 (original) - 20537 (post-filtering) = 2906 (genes removed)
dds_te


####### set factors for statistical analyses
###### Note you need to change condition to treatment (to match our design above)
#  and levels to our treatment names in the PHENO_DATA: Ad_lib is the control, Caloric_Restriction is the treatment group
# example:
#dds$condition <- factor(dds$condition, levels=c("untreated","treated"))
dds_te$condition <- factor(dds_te$treatment, levels=c("control","HS"))


################     1.4 Differential expression analysis
dds_te <- DESeq(dds_te)
res_te <- results(dds_te)
res_te

# We can order our results table by the smallest adjusted p value:
  resOrdered_te <- res_te[order(res_te$padj),]
  resOrdered_te
# We can summarize some basic tallies using the summary function the default is p<0.1.
  summary(res_te)
#How many adjusted p-values were less than 0.1?
  sum(res_te$padj < 0.1, na.rm=TRUE)
#If the adjusted p value will be a value other than 0.1, alpha should be set to that value:
  res05_te <- results(dds_te, alpha=0.05)
  summary(res05_te)
  sum(res05_te$padj < 0.05, na.rm=TRUE)





###    1.5.1 MA-plot
  ## Question 4
  ## plotMA shows the log2 fold changes attributable to a given variable over the meanof normalized counts.
  ## Points will be colored blue if the adjusted p value is less than 0.1.
  ## Points which fall out of the window are plotted as open triangles pointing either up or down
  plotMA(res_te, main="DESeq2", ylim=c(-8,8))

  #After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot.
  # One can then recover the gene identiers by saving the resulting indices:
  #idx <- identify(res$baseMean, res$log2FoldChange)
    # after selecting a gene. You need to press escape to move on
  #rownames(res_te)[idx]


##  1.5.2 Plot counts - sanity check!
  ### Question 5
  # You can select the gene to plot by rowname or by numeric index.
  #plotCounts(dds, gene="FUN_021931", intgroup="treatment") # plotCounts(dds, gene="FUN_022644", intgroup="treatment")
  # You can plot the gene with th lowest adjuated P-value
  plotCounts(dds_te, gene=which.min(res_te$padj), intgroup="treatment")
  dds_te

  ##  Write your results to a file
  write.csv(as.data.frame(resOrdered_te), file="DGESeq_results.csv")

  ## 2.1.2 Extracting transformed values
  rld_te <- rlog(dds_te)
  vsd_te <- varianceStabilizingTransformation(dds_te)
  head(assay(rld_te), 3)

  ### Heatmap of the count matrix
  #library("genefilter")
  topVarGenes_te <- head(order(rowVars(assay(vsd_te)), decreasing = TRUE), 50)

  library("pheatmap")
  mat_te  <- assay(vsd_te)[ topVarGenes, ]
  mat_te  <- mat_te - rowMeans(mat_te)
  anno_te <- as.data.frame(colData(vsd_te)[, c("treatment", "type")])
  df_te <- as.data.frame(colData(dds_te)[,c("treatment","type")])
    pheatmap(mat_te, annotation_col = anno_te)


 ## 2.2.1 Heatmap of the count matrix
#  library("pheatmap")
#  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
#  nt <- normTransform(dds) # defaults to log2(x+1)
#  log2.norm.counts <- assay(nt)[select,]
#   df <- as.data.frame(colData(dds)[,c("treatment","type")])
#  pheatmap(assay(vsd)[mat,], cluster_rows=TRUE, show_rownames=TRUE,
#          cluster_cols=TRUE, annotation_col=df)

# pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
#  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,  cluster_cols=FALSE, annotation_col=df)


  #2.2.2 Heatmap of the sample-to-sample distances
  sampleDists_te <- dist(t(assay(rld_te)))
  library("RColorBrewer")
  sampleDistMatrix_te <- as.matrix(sampleDists_te)
  rownames(sampleDistMatrix_te) <- paste(rld_te$treatment)
  colnames(sampleDistMatrix_te) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix_te,
           clustering_distance_rows=sampleDists_te,
           clustering_distance_cols=sampleDists_te,
           col=colors)


 # 2.2.3 Principal component plot of the samples
  plotPCA(rld_te, intgroup=c("treatment"))
