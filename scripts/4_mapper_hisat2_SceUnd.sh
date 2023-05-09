#!/bin/sh
 
##### Map the cleaned (paired) reads to the indexed reference
##    Use HiSat2 to map your cleaned (paired) reads to the indexed reference
##        Use Stringtie to count the reads mapped to genes and transcripts, defined by the genome annotation file
##        use the python script to take the Stringtie results to make two counts matricies, one at the gene level and one at the transcript level
## HiSat2 Mapping     Input: Cleaned read files, paired (.fastq); Indexed genome
##                    Output: Alignment .sam files  
## Samtools  Convert .sam to .bam and sort         Input: Alignment files  .sam
##                                                 Output: Sorted .bam files
## Stringtie  Counting reads  Input: sorted .bam file
##                            Output:  Directories of counts files for Ballgown (R program for DGE)
##              prepDE.py    Python script to create a counts matrics from the Stringtie output.  Inputs: Directory from Stringtie
##                                                                                                Output:  .csv files of counts matrix
## For running the script on the Alabama Super Computer.
## For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
## After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
## then run the script by using "run_script [script name]"
## suggested paramenters are below to submit this script.
##    queue: class
##    core: 6
##    time limit (HH:MM:SS): default
##    Memory: 12gb
##    run on dmc
###############################################

### Load modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load hisat2
module load stringtie/2.2.1
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

### Define directories
DIR="/scratch/HeatStressConsistency" # the project directory
CLEANDIR="$DIR/CleanData"            # contains trimmed paired read files
FILEDIR="$DIR/lists"                 # contains the file with the list of samples
REFDIR="$DIR/ReferenceGenome"        # contains the indexed reference genome
MAPDIR="$DIR/Map_HiSat2"             # destination for the mapped read (BAM) files
COUNTDIR="$DIR/Counts_StringTie"     # destination for the count matrices
RESULTSDIR="$DIR/Results"            # destination for key results files
SCRIPTDIR="$DIR/scripts"             # contains project scripts


species="ThaEle"

#speciesArray=("ThaEle" "SceUnd")

#for species in ${speciesArray[@]}
#do
  ref="${REFDIR}/${species}_RefGenome"
  
  samplelist="${FILEDIR}/SraAccList${species}.txt"
  
  ## Create output directories for this species if they do not already exist
  mkdir -p $MAPDIR/${species}
  mkdir -p $COUNTDIR/${species}
  mkdir -p $RESULTSDIR/${species}

  ##### Map and Count the Data using HiSAT2 and StringTie
  ## Use xargs to run each command for each sample listed in the samplelist file
  
  cd ${MAPDIR}/${species}

  ## Map reads with HiSat2
  cat ${samplelist} | \
    xargs -I {} --max-args=1 --verbose \
      hisat2 \
        -p 6 \
        --dta \
        --phred33 \
        -x ${ref}_index \
        -1 ${CLEANDIR}/${species}/{}_1_paired.fastq.gz \
        -2 ${CLEANDIR}/${species}/{}_2_paired.fastq.gz \
        -S {}.sam

  ## Convert SAM to BAM
  cat ${samplelist} | \
    xargs -I {} --max-args=1 --verbose \
      samtools view -@ 6 -bS {}.sam -o {}.bam

  ## Sort reads in BAM files
  cat ${samplelist} | \
    xargs -I {} --max-args=1 --verbose \
      samtools sort -@ 6 {}.bam {}_sorted

  ## Create an overview of BAM files
  cat ${samplelist} | \
    xargs -I {} --max-procs=6 --max-args=1 --verbose \
      bash -c 'samtools flagstat "${1}_sorted.bam" > "${1}_stats.txt"' bash {}

### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model. 
  ### The output from StringTie are counts folders in a directory that is ready to bring into the R program Ballgown to 
  ### Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons 

  ## Make a subdirectory for each sample
  cat ${samplelist} | \
    xargs -I {} --max-args=1 --verbose \
      mkdir -p ${COUNTDIR}/${species}/{}

  ## Count the reads mapped to each gene, exon, and transcript model
  cat ${samplelist} | \
    xargs -I {} --max-args=1 --verbose \
      stringtie \
        -p 6 \
        -e \
        -B \
        -G ${ref}.gtf \
        -o ${COUNTDIR}/${species}/{}/{}.gtf \
        -l {} \
        ${MAPDIR}/${species}/{}_sorted.bam


  ## Copy stats files from samtools to results directory
  cp *.txt ${RESULTSDIR}/${species}


  cd ${COUNTDIR}

  ## Convert files in Ballgown folder to count matrix
  python ${SCRIPTDIR}/prepDE.py -i ${species}

  ## Copy count matrices to results directory
  cp *.csv ${RESULTSDIR}/${species}

#done

cd $RESULTSDIR
tar cvzf $DIR/Results_${species}.tar.gz *
