## Alignment pipeline
- [0_1_downloadSRA_QC.sh](0_1_downloadSRA_QC.sh): Download FASTA files from the Sequence Read Archive (SRA) using SRA-toolkit and run FastQC.
- [2_cleanTrim_QC.sh](2_cleanTrim_QC.sh): Trim reads using Trimmomatic and run FastQC.
- [3_hisat2_index.sh](3_hisat2_index.sh): Create a reference index to use for mapping with HiSat2.
- [4_mapper_hisat2.sh](4_mapper_hisat2.sh): Map reads to indexed reference using HiSat2.

## Differential gene expression analysis
- [DESeq2_SceUnd.R](DESeq2_SceUnd.R): Analyze differentially expressed genes in _S. undulatus_ using DESeq2.
- [DESeq2_ThaEle.R](DESeq2_ThaEle.R): Analyze differentially expressed genes in _T. elegans_ using DESeq2.

## Utility scripts
- [calcQuals.sh](calcQuals.sh): Calculate average sequence quality for each FASTQ file.
- [runMultiQC.sh](runMultiQC.sh): Analyze FastQC output grouped by species.
- [GSEAprep.R](GSEAprep.R): Create correctly formatted input files for GSEA.
  - Normalized expression data file (.txt)
  - Ranked expression file (.rnk)
  - Phenotype data file (.cls)
