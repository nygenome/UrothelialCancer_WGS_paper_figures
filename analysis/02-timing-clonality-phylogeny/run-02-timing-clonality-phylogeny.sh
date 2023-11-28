#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper & Timothy R. Chu 

################################################################# /COPYRIGHT ###
################################################################################

## Use "pileup method" for joint calling, use that to generate phylogenetic trees 
## with Lichee, as well as mutation timing and clonality information using MutationTimer 


TNFILE=
PURITY_PLOIDY=
WORKING_DIR=


###################
## Joint calling ##
###################

## TC: joint calling workflow 



##############################
## Build phylogenetic trees ##
##############################

## TC: lichee workflow 



#####################
## Mutation timing ##
#####################

while read tumor normal gender patient; do

    tn=$tumor--$normal

    ## Inputs
    vcf=$VARIANT_DIR/$patient/$tn.highconf.ffpefilter.union.vcf
    jba=$JABBA_DIR/jabba/$tn/jabba.simple.rds
    purity=$(awk -v FS=',' -v tumor=$tumor '($1 == tumor) {print $2}' $PURITY_PLOIDY)

    ## Case where this patient only has a single sample
    if [[ ! -f $vcf ]]; then
        vcf=$PROJECT_DIR/compbio/Somatic_analysis/$tn/$tn.snv.indel.final.v6.annotated.ffpePON_filtered.vcf
    fi

    [[ ! -f $jba ]] && continue

    ## Prepped input files
    ascn=$WORKING_DIR/cn/$tn-jabba-ascn.rds

    ## Output files
    out_file_vcf=$WORKING_DIR/mutationtimer/$tn-mutationtimer-snv.vcf
    out_file_bed=$WORKING_DIR/mutationtimer/$tn-mutationtimer-cnv.bed
    out_file_rds=$WORKING_DIR/mutationtimer/$tn-mutationtimer.rds
    out_file_txt=$WORKING_DIR/mutationtimer/$tn-mutationtimer.summary.txt
    out_file_png=$FIG_DIR/$tn-mutationtimer.png

    Rscript $SRCDIR/init-mutation-timer-input.r \
    --jba=$jba \
    --purity=$purity \
    --tumor=$tumor \
    --gender=$gender \
    --out_file_cn=$ascn

    Rscript $SRCDIR/run-mutation-timer.r \
    --cn=$ascn \
    --purity=$purity \
    --tumor=$tumor \
    --gender=$gender \
    --vcf=$vcf \
    --out_file=$out_file_rds

    Rscript $SRCDIR/postprocess-mutation-timer.r \
    --cn=$ascn \
    --mt=$out_file_rds \
    --tumor=$tumor \
    --purity=$purity \
    --gender=$gender \
    --vcf=$vcf \
    --out_file_vcf=$out_file_vcf \
    --out_file_bed=$out_file_bed \
    --out_file_summary=$out_file_txt \
    --out_file_pdf=$out_file_png

done < $TNFILE