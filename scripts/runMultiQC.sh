#! /bin/bash

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/3-2019.10

cd /scratch/HeatStressConsistency

multiqc PreCleanQuality/ThaEle --outdir QCresults --title pre_ThaEle
multiqc PreCleanQuality/SceUnd --outdir QCresults --title pre_SceUnd

multiqc PostCleanQuality/ThaEle --ignore "*unpaired*" --outdir QCresults --title post_ThaEle
multiqc PostCleanQuality/SceUnd --ignore "*unpaired*" --outdir QCresults --title post_SceUnd
