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


## This bash script contains scripts used to generate JaBbA CNV/SV calls 
## Note: The steps prior to postprocessing require a Slurm cluster environment
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))

## Dirs
JABBA_DIR=/path/to/jabba/output/dir

## Files 
TNFILE=tumor_normal_pairs.jabba_analysis.txt
GENE_BED=/gpfs/commons/resources/GRCh38_full_analysis_set_plus_decoy_hla/internal/ensembl_genes_unique_sorted.final.v93.short.chr.sorted.bed



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
while read tumor normal gender; do
    
    tn=$tumor--$normal

    in_file=$JABBA_DIR/jabba/$tn/jabba.events.footprints.rds
    out_file=$JABBA_DIR/jabba/$tn/jabba.events.genes.tab

    Rscript $SRCDIR/jabba/calling/intersect-footprints-with-genes.r \
                --gene_bed=$GENE_BED \
                --in_file=$in_file \
                --out_file=$out_file

done < $TNFILE


## Generate gene-level copy number profiles
while read tumor normal gender; do

    tn=$tumor--$normal
    bed=$JABBA_DIR/jabba/$tn/jabba.simple.cnv.bed
    out_file=$JABBA_DIR/gene-cn/$tn.jabba_gene_cn.bed

    Rscript $SRCDIR/jabba-gene-cn.r
                --in_file=$bed \
                --genes=$GENE_BED \
                --out_file=$out_file

done < $TNFILE
