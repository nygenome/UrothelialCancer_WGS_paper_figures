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

## Run dndscv (SNV/INDEL) and FishHook (SV) to look for evidence of positive selection

## Dirs 
JABBA_DIR= 
ENCODE_DIR=
DNDSCV_DIR=
FISHHOOK_DIR=

## Files 
TNFILE=
CHR_LEN=
JABBA_BLOCK_LIST=
REPLICATION_TIMING=
FRAGILE_SITES=
CYTOBAND=
REPEATMASKER=
CGC=



############
## dndscv ##
############

## Pulled from https://github.com/im3sanger/dndscv_data/tree/master/data 
DNDSCV_REF=$RESOURCE_DIR/RefCDS_human_GRCh38_GencodeV18_recommended.rda
DNDSCV_COV=$RESOURCE_DIR/covariates_hg19_hg38_epigenome_pcawg.rda


## Build list of mutations
MUTATIONS=$WORKING_DIR/highconf-snvs.txt
echo -e "sampleID\tpatient\tchr\tpos\tref\tmut" > $MUTATIONS

while read tumor normal gender patient; do 
  
  tn=$tumor--$normal

  ## Grab union vcf where possible 
  vcf=$VARIANT_DIR/$patient/$tn.highconf.ffpefilter.union.vcf

  ## Case where this patient only has a single sample
  if [[ ! -f $vcf ]]; then
    vcf=$PROJECT_DIR/compbio/Somatic_analysis/$tn/$tn.snv.indel.final.v6.annotated.ffpePON_filtered.vcf
  fi

  ## Select high confidence / added mutations 
  grep -v '^#' $vcf | \
  grep 'HighConfidence\|TYPE=ADDED' | \
  awk -v OFS='\t' -v sample=$tn -v patient=$patient '{print sample, patient, $1, $2, $4, $5}' \
  >> $MUTATIONS

done < $TNFILE


## Run dNdScv
Rscript $SRCDIR/run-dndscv.r \
  --in_file=$MUTATIONS \
  --metadata=$METADATA \
  --ref=$DNDSCV_REF \
  --cov=$DNDSCV_COV \
  --out_file=$OUT_DIR/dndscv-results.rds
  


##############
## FishHook ##
##############

FDR=0.25
BIN_SIZE=1E5
DEDUP='all'

## Build allow list 
Rscript $SRCDIR/init-eligible-regions.r \
  --chr_len=$CHR_LEN \
  --mappability=$MAPPABILITY \
  --block=$JABBA_BLOCK_LIST \
  --out_file=$FISHHOOK_DIR/GRCh38-eligible-regions.bed

## Merge jabba junctions into a single file 
Rscript $SRCDIR/merge-jabba-junctions.r \
  --in_dir=$JABBA_DIR \
  --tn=$TNFILE \
  --out_file=$FISHHOOK_DIR/jabba-junctions.txt


## Assemble covariates
Rscript $SRCDIR/init-covariates.r \
          --chr_len=$CHR_LEN \
          --encode_dir=$ENCODE_DIR \
          --rep_timing=$REPLICATION_TIMING \
          --fragile_sites=$FRAGILE_SITES \
          --cytoband=$CYTOBAND \
          --repeatmasker=$REPEATMASKER \
          --bin_size=$BIN_SIZE \
          --step_size=$bin_size \
          --out_dir=$FISHHOOK_DIR/resources

if [[ $DEDUP == 'dedup' ]]; then
    dedup_flag='TRUE'   
else
    dedup_flag='FALSE'
fi

## Run FishHook
Rscript $SRCDIR/run-fishhook.r \
            --in_file=$FISHHOOK_DIR/jabba-junctions.txt \
            --chr_len=$CHR_LEN \
            --eligible=$FISHHOOK_DIR/GRCh38-eligible-regions.bed \
            --bin_size=$BIN_SIZE \
            --step_size=$BIN_SIZE \
            --one_patient_per_locus=$dedup_flag \
            --covariate_dir=$FISHHOOK_DIR/resources \
            --out_file=$FISHHOOK_DIR/fishhook-result.$BIN_SIZE.$dedup.rds



## Annotation 
in_file=$FISHHOOK_DIR/fishhook-result.$BIN_SIZE.$dedup.rds
out_file_1=$FISHHOOK_DIR/fishhook-result.$BIN_SIZE.$dedup.annotated.txt
out_file_2=$FISHHOOK_DIR/fishhook-result.$BIN_SIZE.$dedup.annotated.dispersion.txt
out_file_bed=$FISHHOOK_DIR/fishhook-hits.$BIN_SIZE.$dedup.bed

Rscript $SRCDIR/annotate-fishhook-results.r \
    --in_file=$in_file \
    --genes=$CGC \
    --cytoband=$CYTOBAND \
    --fdr=$FDR \
    --out_file=$out_file_1

Rscript $SRCDIR/init-dispersion-score.r \
    --in_file=$out_file_1 \
    --junctions=$WORKING_DIR/jabba-junctions.txt \
    --one_patient_per_locus=$dedup_flag \
    --out_file=$out_file_2

Rscript $SRCDIR/init-export-srb-to-bed.r \
    --in_file=$out_file_2 \
    --fdr=$FDR \
    --out_file=$out_file_bed
