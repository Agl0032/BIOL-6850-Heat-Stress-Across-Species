# BIOL-6850-Heat-Stress-Across-Species

Final Project for the Functional Genomics course (BIOL 6850) at Auburn University. 

### Starting directory structure

```
HeatStressConsistency/
├── ReferenceGenome/
│   ├── SceUnd_RefGenome.fasta
│   ├── SceUnd_RefGenome.gff 
│   ├── ThaEle_RefGenome.fasta
│   └── ThaUnd_RefGenome.gff
└── scripts/
|   ├── 0_1_downloadSRA_QC.sh
|   ├── 2_cleanTrim_QC.sh
|   ├── 3_hisat2_index.sh
|   ├── 4_mapper_hisat2_SceUnd.sh
|   └── utilities/
|       └── calcQuals.sh
|       └── runMultiQC.sh
└── SraAccList_SceUnd.txt
└── SraAccList_ThaEle.txt
```


### Sample SRA IDs and treatments

#### _S. undulatus_
| Sample | Treatment |
| ------ | --------- |
| SRR6825935 | control |
| SRR6825934 | control |
| SRR6825939 | heat stress |
| SRR6825938 | control |
| SRR6825941 | heat stress |
| SRR6825940 | heat stress |
| SRR6825933 | heat stress |
| SRR6825932 | heat stress |
| SRR6825926 | control |
| SRR6825925 | control |
| SRR6825930 | control |
| SRR6825928 | heat stress |



#### _T. elegans_
| Sample | Treatment | Ecotype |
| ------ | --------- | ------- |
| SRR629599 | control | meadow |
| SRR629651 | heat stress | lakeshore |
| SRR629652 | control | meadow |
| SRR629653 | heat stress | lakeshore |
| SRR629654 | control | lakeshore |
| SRR629655 | heat stress | lakeshore |
| SRR629656 | control | meadow |
| SRR629657 | heat stress | lakeshore |
| SRR629658 | control | meadow |
| SRR629659 | heat stress | lakeshore |
| SRR629661 | control | meadow |
| SRR629660 | heat stress | meadow |
| SRR629662 | control | lakeshore |
| SRR629664 | heat stress | meadow |
| SRR629663 | control | lakeshore |
| SRR629665 | control | lakeshore |
| SRR629666 | heat stress | meadow |
| SRR629667 | control | lakeshore |
| SRR629668 | heat stress | meadow |
| SRR629669 | heat stress | meadow |
