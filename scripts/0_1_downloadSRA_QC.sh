#! /bin/bash

##### Download FASTQ files for all samples from SRA using SRA-toolkit and analyze quality with FastQC

############## Load modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load sra/2.10.9
module load fastqc/0.11.9


############## Directories
DIR="/scratch/HeatStressConsistency"
SRADIR="$DIR/sra"
RAWDIR="$DIR/RawData/"
CLEANDIR="$DIR/CleanData/"
PREDIR="$DIR/PreCleanQuality/"
POSTDIR="$DIR/PostCleanQuality/"

speciesArray=("ThaEle" "SceUnd") # garter snake, fence lizard

for species in ${speciesArray[@]}
do
  mkdir -p ${RAWDIR}/${species}
  mkdir -p ${CLEANDIR}/${species}
  mkdir -p ${PREDIR}/${species}
  mkdir -p ${POSTDIR}/${species}


  ############## SRA toolkit to download sequences (in parallel using xargs)
  cd $SRADIR

  cat SraAccList${species}.txt | xargs --max-procs=4 --max-args=1 prefetch
  cat SraAccList${species}.txt | xargs --max-procs=4 --max-args=1 fastq-dump --outdir "${RAWDIR}/${species}" --readids --split-files

  cd ${RAWDIR}/${species}
  echo *.fastq | xargs --max-procs=4 --max-args=1 --verbose gzip 


  ############## FASTQC to assess quality of the sequence data
  ## FastQC: run on each of the data files that have 'All' to check the quality of the data
  ## The output from this analysis is a folder of results and a zipped file of results

  echo *.fastq.gz | xargs --max-procs=4 --max-args=1 --verbose fastqc --outdir=${PREDIR}/${species}

done


#######  Tarball the directory containing the FASTQC results
cd $DIR
tar cvzf PreCleanQuality.tar.gz PreCleanQuality/*
