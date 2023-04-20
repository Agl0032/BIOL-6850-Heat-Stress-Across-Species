#!/bin/sh
 
######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to 
##    Use HiSat2 to index your reference genome and then map your cleaned (paired) reads to the indexed reference
##              First need to use gffread to convert annotation file from .gff3 to .gft formate
##              Use Stringtie to count the reads mapped to genes and transcripts, defined in this case by the genome annotation file
##              use the python script to take the Stringtie results to make two counts matricies, one at the gene level and one at the transcript level
## HiSat2  Indexing  InPut: Reference genome file (.fasta), and annotation file (.gff3) (Optional)
##                    Output: Indexed genome 
## HiSat2 Mapping     Input: Cleaned read files, paired (.fasq); Indexed genome
##                    Output: Alignment .sam files  
## Samtools  Convert .sam to .bam and sort         Input: Alignment files  .sam
##                                                  Output: Sorted .bam files
## Stringtie  Counting reads  Input: sorted .bam file
##                            Output:  Directories of counts files for Ballgown (R program for DGE)
##              prepDE.py    Python script to create a counts matrics from the Stringtie output.  Inputs: Directory from Stringtie
##                                                                                                Output:  .csv files of counts matrix
## For running the script on the Alabama Super Computer.
##  For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
##  After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
##  then run the script by using "run_script [script name]"
##  suggested paramenters are below to submit this script.
##    queue: class
##    core: 6
##    time limit (HH:MM:SS): 04:00:00 
##    Memory: 12gb
##    run on dmc
###############################################

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load hisat2
module load stringtie/2.2.1
module load gffcompare
module load python/2.7.1
module load gcc/9.3.0
module load samtools
module load bcftools/1.2
module load gffread/
module load gffcompare/


#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x

#speciesArray=("ThaEle" "SceUnd")

DIR="/scratch/HeatStressConsistency"
CLEANDIR="$DIR/CleanData"
FILEDIR="$DIR/lists"
REFDIR="$DIR/ReferenceGenome"
MAPDIR="$DIR/Map_HiSat2"
COUNTDIR="$DIR/Counts_StringTie"
RESULTSDIR="$DIR/Results"
SCRIPTDIR="$DIR/scripts"

species="ThaEle"

#for species in ${speciesArray[@]}
#do
  ref="${REFDIR}/${species}_RefGenome"

  mkdir -p $MAPDIR/${species}
  mkdir -p $COUNTDIR/${species}
  mkdir -p $RESULTSDIR/${species}

  ##################  Prepare the Reference Index for mapping with HiSat2 #############################
  cd $REFDIR

  ###  Identify exons and splice sites
  gffread $ref.gff -T -o $ref.gtf               ## Convert annotation from .gff to .gft
  extract_splice_sites.py $ref.gtf > $ref.ss
  extract_exons.py $ref.gtf > $ref.exon

  gunzip ${ref}.fasta.gz

  #### Create a HISAT2 index for the reference genome
  hisat2-build -p 16 --ss $ref.ss --exon $ref.exon $ref.fasta ${ref}_index

