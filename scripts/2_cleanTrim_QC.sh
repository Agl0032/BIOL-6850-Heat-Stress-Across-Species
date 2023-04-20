#! /bin/bash

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load trimmomatic/0.39
module load fastqc/0.11.9

speciesArray=("ThaEle" "SceUnd")

DIR="/scratch/HeatStressConsistency"
RAWDIR="$DIR/RawData"
CLEANDIR="$DIR/CleanData"
POSTDIR="$DIR/PostCleanQuality"
FILEDIR="$DIR/lists"


############## Trimmomatic

for species in ${speciesArray[@]}
do
  cd ${RAWDIR}/${species}

  cat ${FILEDIR}/SraAccList${species}.txt | \
    xargs -I {} \
    java -jar /mnt/beegfs/home/aubmxa/.conda/envs/BioInfo_Tools/share/trimmomatic-0.39-1/trimmomatic.jar PE \
      -threads 6 \
      -phred33 \
      {}_1.fastq.gz {}_2.fastq.gz  \
      ${CLEANDIR}/${species}/{}_1_paired.fastq.gz ${CLEANDIR}/${species}/{}_1_unpaired.fastq.gz \
      ${CLEANDIR}/${species}/{}_2_paired.fastq.gz ${CLEANDIR}/${species}/{}_2_unpaired.fastq.gz \
      ILLUMINACLIP:${FILEDIR}/AdaptersToTrim_All.fa:2:35:10 \
      HEADCROP:10 \
      LEADING:30 \
      TRAILING:30 \
      SLIDINGWINDOW:6:30 \
      MINLEN:36
done


############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results

cd $CLEANDIR

echo */*.fastq.gz | xargs --max-procs=4 --max-args=1 --verbose fastqc --outdir=$POSTDIR

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
## when finished use scp or rsync to bring the tarballed .gz results file to your computer and open the .html file to evaluate the quality of your raw data.
cd $DIR
tar cvzf PostCleanQuality.tar.gz $POSTDIR/*
