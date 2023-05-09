#! /bin/bash

##### FASTQC to assess quality of the sequence data
## Only on files that haven't had it yet by checking for existing FastQC output files

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.11.9

speciesArray=("ThaEle" "SceUnd")    # not currently used
SPECIES="ThaEle"

DIR="/scratch/HeatStressConsistency"
CLEANDIR="$DIR/CleanData"
POSTDIR="$DIR/PostCleanQuality"
FILEDIR="$DIR/lists"


cd $CLEANDIR

## Make an array of all paired trimmed files
gzfiles=$(for file in $(echo $SPECIES/*_paired.fastq.gz); do file=${file%%.fastq.gz}; file=${file##*/}; echo $file; done)

## Initialize an array to store the files which have not been through FastQC
needfiles=()

## For every trimmed file, check if there is a corresponding FastQC output
## If there isn't, add that file to the list to be run
for file in ${gzfiles[@]}
do
  if [ ! -e $POSTDIR/$SPECIES/${file}_fastqc.zip ]; then
    needfiles+=(${file})
  fi
done

## Run FastQC where needed
for file in ${needfiles[@]}; do echo $file; done | \
  xargs -n 1 -I {} echo $SPECIES/{}.fastq.gz | \
  xargs --max-procs=4 --max-args=1 --verbose fastqc --outdir=$POSTDIR/$SPECIES

##### Tarball the directory containing the FastQC results
cd $DIR
tar cvzf PostCleanQuality.tar.gz PostCleanQuality/*
