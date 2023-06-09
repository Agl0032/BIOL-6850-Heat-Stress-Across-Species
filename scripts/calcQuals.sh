#!/bin/bash

## Uses the "fastqc_data.txt" files generated by FastQC to extract count of sequences and calculate the average
## sequence quality for each sample before and after cleaning. Only considers paired sequences for cleaned data.
## Assumes the presence of both a pre- and post-trimming subfolder with the following names:
##	PreCleanQuality - the FastQC outputs for the samples before cleaning (i.e., with Trimmomatic)
##	PostClean Quality - the FastQC outputs for the samples after cleaning
## Also assumes the presence of a species subfolder for each species in speciesArray in each of the pre and post trim
## directories. Samples present in preclean folder must also be present in postclean folder. Outputs a table of values
## for each species in a csv file starting with "quality-table", which can be changed by changing OutPrefix.

## If all *_fastqc files are zipped, run on the command line before running the script
## "echo *CleanQuality/*/*.zip | xargs -n 1 unzip"

OutPrefix="quality-table"
speciesArray=("ThaEle" "SceUnd")

for species in ${speciesArray[@]}; do
  echo "${species}"
  OUTFILE="${OutPrefix}-${species}.csv"

  echo "Sample,Pre Seq Count,Post Seq Count,Pre Quality,Post Quality" > "${OUTFILE}"
  echo -e "Sample\tPre Seq Count\tPost Seq Count\tPre Quality\tPost Quality"

  for preFile in PreCleanQuality/${species}/*/*fastqc_data.txt; do 
    ## EX: preFile="PreCleanQuality/ThaEle/SRR6819014_1_fastqc/fastqc_data.txt"

    sampleID=${preFile#*/*/} 			  #=> "SRR6819014_1_fastqc/fastqc_data.txt" 
    sampleID=${sampleID%%_fastqc/fastqc_data.txt} #=> "SRR6819014_1"

    postFile="PostCleanQuality/${species}/${sampleID}_paired_fastqc/fastqc_data.txt"

    preSeqCount=$(awk '/Total Sequences/ {print $3}' $preFile)
    postSeqCount=$(awk '/Total Sequences/ {print $3}' $postFile)

    preAvgQual=$(awk '/>>Per sequence quality scores/,/>>END_MODULE/' $preFile | awk '/[[:digit:]]/ {sum += $1 * $2; count +=$2} END {print sum/count}')
    postAvgQual=$(awk '/>>Per sequence quality scores/,/>>END_MODULE/' $postFile | awk '/[[:digit:]]/ {sum += $1 * $2; count +=$2} END {print sum/count}')

    echo "$sampleID,$preSeqCount,$postSeqCount,$preAvgQual,$postAvgQual" >> "${OUTFILE}"
    echo -e "$sampleID\t$preSeqCount\t$postSeqCount\t$preAvgQual\t$postAvgQual"

  done
done
