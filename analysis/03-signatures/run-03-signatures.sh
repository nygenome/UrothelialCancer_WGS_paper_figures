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

## Compute signatures via deconstructSigs from joint-called SNV/INDELs
## Probabilistically assign signatures to individual mutations, incorporating 
## some knowledge-based filtering 

## Dirs
DECONSTRUCTSIGS_DIR=
MUTATIONTIMER_DIR= 
SIG_ASSIGN_DIR= 
TCGA_DIR=

## Files 
METADATA=
TCGA_MANIFEST=
SBS_REF=


#####################
## deconstructSigs ##
#####################

## TC: deconstructSigs workflow 


######################################
## Signature assignment + filtering ##
######################################

## If a mutation appears even once in the TCGA cohort, exclude
## it from chemo-assignment (non of these tumors are chemo treated)
SBS31_35_FILTER_LIST=$SIG_ASSIGN_DIR/tcga-blca-filter-list.txt
while read id filename md5 size state; do

  vcf=$TCGA_DIR/$id/$filename
  zgrep -v '^#' $vcf | awk '$7 == "PASS" && length($4) == length($5)' | cut -f1,2,4,5 

done < <(tail -n +2 $TCGA_MANIFEST) | sort -u > $SBS31_35_FILTER_LIST

SBS=

## Assign and apply single-sample filters
while read tumor normal gender patient; do

    tn=$tumor--$normal
    vcf=$MUTATIONTIMER_DIR/mutationtimer/$tn-mutationtimer-snv.vcf
    out_vcf=$SIG_ASSIGN_DIR/unfiltered/$tn.mutationtimer-snv.withSigProbs.vcf
    out_vcf_filtered=$OUT_DIR/unfiltered/$tn.mutationtimer-snv.withSigProbs.singleSampleFilter.vcf

    Rscript $SRCDIR/assign-signatures.r \
                --vcf=$vcf \
                --sbs=$SBS \
                --sample_id=$tumor \
                --sbs_ref=$SBS_REF \
                --out_file=$out_vcf

    Rscript $SRCDIR/filter-signatures.r \
            --in_file=$out_vcf \
            --sbs_ref=$SBS_REF \
            --sbs_31_35_filter=$SBS31_35_FILTER_LIST \
            --metadata=$METADATA \
            --tumor_id=$tumor \
            --out_file=$out_vcf_filtered
  
done < $TNFILE



## Apply a patient-level filter for chemotherapy assigned mutations
## I.e., if a mutation appears in a pre-chemo sample, it shouldn't be 
## assigned to SBS31/SBS35 in a post-chemo sample from the same patient
while read patient; do

    ## No need to filter if this patient only has a single sample
    if [[ $(awk -v p=$patient '$1 !~ /^#/ && $4 == p' $TNFILE | wc -l) == 1 ]]; then

        tn=$(awk -v p=$patient '$1 !~ /^#/ && $4 == p {print $1"--"$2}' $TNFILE)

        in_vcf=$SIG_ASSIGN_DIR/unfiltered/$tn.mutationtimer-snv.withSigProbs.singleSampleFilter.vcf
        out_vcf=$SIG_ASSIGN_DIR/$tn.mutationtimer-snv.withSigProbs.filtered.vcf

        cp $in_vcf $out_vcf

    else

        Rscript $SRCDIR/patient-level-chemo-signature-filter.r \
                    --in_dir=$SIG_ASSIGN_DIR/unfiltered \
                    --patient=$patient \
                    --tn_file=$TNFILE \
                    --out_dir=$SIG_ASSIGN_DIR \
                    --in_suffix=.mutationtimer-snv.withSigProbs.singleSampleFilter.vcf \
                    --out_suffix=.mutationtimer-snv.withSigProbs.filtered.vcf \
                    --metadata=$METADATA

    fi

done < <(grep -v '^#' $TNFILE | cut -f4 | sort -u)



## Count mutations assigned to each timing/signature combination 
## Used for fold-change boxplot
 while read tumor normal gender patient; do

  tn=$tumor--$normal
  vcf=$SIG_ASSIGN_DIR/$tn.mutationtimer-snv.withSigProbs.filtered.vcf
  
  Rscript $SRCDIR/vcf-to-counts.r \
        --in_file=$vcf \
        --sbs_ref=$SBS_REF \
        --out_file=$SIG_ASSIGN_DIR/counts/$tn-signature-clonality-counts.csv
  
done < $TNFILE



#########################
## Mutational velocity ##
#########################

## TC: mutational velocity workflow 



################################
## Normal urothelium analysis ##
################################

## TC: normal urothelium pileup workflow 