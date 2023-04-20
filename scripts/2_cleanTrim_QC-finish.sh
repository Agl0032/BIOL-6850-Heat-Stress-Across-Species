#! /bin/bash

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load fastqc/0.11.9

speciesArray=("ThaEle" "SceUnd")
SPECIES="ThaEle"

DIR="/scratch/HeatStressConsistency"
RAWDIR="$DIR/RawData"
CLEANDIR="$DIR/CleanData"
POSTDIR="$DIR/PostCleanQuality"
FILEDIR="$DIR/lists"


############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results

cd $CLEANDIR

#qcfiles=$(for file in $(echo $POSTDIR/$SPECIES/*.zip); do file=${file%%_fastqc.zip}; file=${file##*/}; echo $file; done)
gzfiles=$(for file in $(echo $SPECIES/*.fastq.gz); do file=${file%%.fastq.gz}; file=${file##*/}; echo $file; done)
needfiles=()

for file in ${gzfiles[@]}
do
  if [ ! -e $POSTDIR/$SPECIES/${file}_fastqc.zip ]; then
    needfiles+=(${file})
  fi
done

for file in ${needfiles[@]}; do echo $file; done | xargs -n 1 -I {} echo $SPECIES/{}.fastq.gz | xargs --max-procs=4 --max-args=1 --verbose fastqc --outdir=$POSTDIR/$SPECIES

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
## when finished use scp or rsync to bring the tarballed .gz results file to your computer and open the .html file to evaluate the quality of your raw data.
cd $DIR
tar cvzf PostCleanQuality.tar.gz PostCleanQuality/*
