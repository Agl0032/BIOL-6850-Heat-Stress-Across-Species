#!/bin/sh

##### Index the reference genome for alignment with HiSat2
## Run this script once for each species - change the species variable accordingly.
##   gffread: Convert annotation file from .gff3 to .gft format
##   HiSat2 Indexing  Input: Reference genome file (.fasta), and annotation file (.gff3) (Optional)
##                    Output: Indexed genome

## For running the script on the Alabama Super Computer.
## For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
## After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
## then run the script by using "run_script [script name]"
## suggested paramenters are below to submit this script.
##    queue: bigmem
##    core: 16
##    time limit (HH:MM:SS): default
##    Memory: 160gb
##    run on dmc
###############################################

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load gffread
module load hisat2
module load gffcompare
module load python/2.7.1
module load gcc/9.3.0

#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x

DIR="/scratch/HeatStressConsistency"
REFDIR="$DIR/ReferenceGenome"

### Comment out one of the species to run the other
species="ThaEle"
# species="SceUnd"

### OR use an array and a for loop
#speciesArray=("ThaEle" "SceUnd")

#for species in ${speciesArray[@]}
#do
  ref="${REFDIR}/${species}_RefGenome"

  ##### Prepare the Reference Index for mapping with HiSat2
  cd $REFDIR

  ###  Identify exons and splice sites
  gffread $ref.gff -T -o $ref.gtf               ## Convert annotation from .gff to .gft
  extract_splice_sites.py $ref.gtf > $ref.ss
  extract_exons.py $ref.gtf > $ref.exon
  
  # The reference FASTA needs to be unzipped for HiSat2
  gunzip ${ref}.fasta.gz

  #### Create a HISAT2 index for the reference genome
  hisat2-build -p 16 --ss $ref.ss --exon $ref.exon $ref.fasta ${ref}_index

