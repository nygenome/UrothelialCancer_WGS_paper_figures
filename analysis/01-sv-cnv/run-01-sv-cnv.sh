#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper

################################################################# /COPYRIGHT ###
################################################################################


## This bash script contains scripts used to generate JaBbA CNV/SV calls, 
## followed by some postprocessing and AmpliconArchitect 
## Note: The steps prior to postprocessing require a Slurm cluster environment
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))



###################
## Preprocessing ##
###################

## Run fragcounter to get read depth profiles in 1KB bins 
bash $SRCDIR/jabba/preprocessing.sh fragcounter $SRCDIR/jabba_config.txt

## Generate dryclean PON from normal read depth profiles
bash $SRCDIR/jabba/preprocessing.sh dryclean-pon $SRCDIR/jabba_config.txt

## De-noise tumor read depth profiles with dryclean 
bash $SRCDIR/jabba/preprocessing.sh dryclean $SRCDIR/jabba_config.txt



###############
## Run JaBbA ##
###############

## Run JaBbA to infer junction-balanced genome graphs
bash $SRCDIR/jabba/jabba.sh $SRCDIR/jabba_config.txt

## Run calling to infer simple/complex SVs
bash $SRCDIR/jabba/calling.sh $SRCDIR/jabba_config.txt



####################
## Postprocessing ##
####################

## Get genes overlapping with events
while read tumor normal gender patient; do
    
    tn=$tumor--$normal

    in_file=$JABBA_DIR/jabba/$tn/jabba.events.footprints.rds
    out_file=$JABBA_DIR/jabba/$tn/jabba.events.genes.tab

    Rscript $SRCDIR/jabba/calling/intersect-footprints-with-genes.r \
                --gene_bed=$GENE_BED \
                --in_file=$in_file \
                --out_file=$out_file

done < $TNFILE


## Generate gene-level copy number profiles
while read tumor normal gender patient; do

    tn=$tumor--$normal
    bed=$JABBA_DIR/jabba/$tn/jabba.simple.cnv.bed
    out_file=$JABBA_DIR/gene-cn/$tn.jabba_gene_cn.bed

    Rscript $SRCDIR/jabba-gene-cn.r
                --in_file=$bed \
                --genes=$GENE_BED \
                --out_file=$out_file

done < $TNFILE



#######################
## AmpliconArchitect ##
#######################

while read tumor normal; do 

    tn=$tumor--$normal

    $AA_REPO/AmpliconSuite-pipeline.py \
    -o outputs/$tn/ \
    -s $tn \
    -t 4 \
    --cnv_bed $JABBA_DIR/jabba/$tn/jabba.simple.cnv.bed \
    --bam $PROJECT_DIR/Sample_$tumor/analysis/$tumor.final.bam \
    --run_AA \
    --run_AC \
    --cngain 4 \
    --cnsize_min 10000 \
    --downsample 10

done < $TNFILE