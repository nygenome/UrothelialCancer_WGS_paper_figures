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

## Identify kataegis via SigProfilerClusters, nominate ecDNAs using JaBbA/AA data, 
## then identify kyklonas



#########################
## SigProfilerClusters ##
#########################
python run_Sigprofiler_SimClusters_vaf.py $project $SIGPROFILER_DIR



####################
## Postprocessing ##
####################

APOBEC_ASSOC_CUTOFF=0.75
dir='vcf_hc_vaf'

      
WORKING_DIR=$SIGPROFILER_DIR/$dir
FIG_DIR=$WORKING_DIR/fig
LOG_DIR=$WORKING_DIR/logs
VAF_DIR=$WORKING_DIR/output/ecdna-mutations-by-vaf

mkdir -p $FIG_DIR $LOG_DIR
mkdir -p $WORKING_DIR/output/clustered/kataegis-annotated \
        $WORKING_DIR/output/nonClustered/annotated \
        $FIG_DIR \
        $LOG_DIR \
        $VAF_DIR


## Export per-sample IMD cutoffs
python -mpickle $WORKING_DIR/output/simulations/data/imds.pickle > $WORKING_DIR/output/sample-imds.json 

## Annotate kataegis VCFs with signature info
while read tumor normal gender patient; do

    tn=$tumor--$normal
    in_file_ktg=$WORKING_DIR/output/clustered/kataegis/$tn-mutationtimer-snv.vcf
    in_file_nc=$WORKING_DIR/output/nonClustered/SBS/$tn-mutationtimer-snv.vcf
    in_file_sig=$SIG_ASSIGNMENT_DIR/$tn.mutationtimer-snv.withSigProbs.filtered.vcf
    in_file_jba=$JABBA_DIR/$tn/jabba.events.rds

    out_file_ktg=$WORKING_DIR/output/clustered/kataegis-annotated/$tn-mutationtimer-snv.annotated.txt
    out_file_nc=$WORKING_DIR/output/nonClustered/annotated/$tn-mutationtimer-snv.annotated.txt

    ## If no kataegis results, just skip
    [[ ! -f $in_file_ktg ]] && continue

    Rscript $SRCDIR/init-annotate-kataegis-results.r \
                --in_file_ktg=$in_file_ktg \
                --in_file_sig=$in_file_sig \
                --in_file_jba=$in_file_jba \
                --ktg_flag=$dir \
                --tumor=$tumor \
                --apobec_cutoff=$APOBEC_ASSOC_CUTOFF \
                --out_file=$out_file_ktg

    Rscript $SRCDIR/init-annotate-nonclustered-results.r \
                --in_file_nc=$in_file_nc \
                --in_file_sig=$in_file_sig \
                --in_file_jba=$in_file_jba \
                --tumor=$tumor \
                --apobec_cutoff=$APOBEC_ASSOC_CUTOFF \
                --out_file=$out_file_nc

    done < $TNFILE



## Cross kataegis results with jabba results
## init-amplicon-kataegis-merge.r generates putative ecDNA footprints
while read tumor normal gender patient; do

    tn=$tumor--$normal

    in_file_ktg=$WORKING_DIR/output/clustered/kataegis-annotated/$tn-mutationtimer-snv.annotated.txt
    in_file_nc=$WORKING_DIR/output/nonClustered/annotated/$tn-mutationtimer-snv.annotated.txt
    in_file_jba=$JABBA_DIR/$tn/jabba.events.footprints.txt
    
    in_file_aa=$AA_DIR/$tn.amplicons.10kb_4X_downsample10.bed

    ## This case failed with the 10X/4CN/10kb run, use original run 
    if [[ $tumor == 'PM636-Z2-1-Case-WGS' || $tumor == 'PM117-Z10-1-Case-WGS' ]]; then 
        in_file_aa=$AA_DIR/$tn.amplicons.5X_downsample10.bed
    fi

    out_file_ktg=$WORKING_DIR/output/clustered/kataegis-annotated/$tn-mutationtimer-snv.annotated.jabba.txt
    out_file_nc=$WORKING_DIR/output/nonClustered/annotated/$tn-mutationtimer-snv.annotated.jabba.txt
    out_file_jba=$JABBA_DIR/$tn/jabba.events.footprints.kataegis.$dir.txt
    

    Rscript $SRCDIR/init-amplicon-kataegis-merge.r \
                --in_file_ktg=$in_file_ktg \
                --in_file_jba=$in_file_jba \
                --in_file_aa=$in_file_aa \
                --out_file_ktg=$out_file_ktg \
                --out_file_jba=$out_file_jba


    Rscript $SRCDIR/init-amplicon-nonclustered-merge.r \
                --in_file_nc=$in_file_nc \
                --in_file_jba=$out_file_jba \
                --out_file_nc=$out_file_nc


done < $TNFILE  



## Kyklonas summary stats
Rscript $SRCDIR/summarize-kyklonas-stats.r \
            --in_dir_ktg=$WORKING_DIR/output/clustered/kataegis-annotated \
            --in_dir_nc=$WORKING_DIR/output/nonClustered/annotated \
            --tn_file=$TNFILE


## ecDNA size distribution
Rscript $SRCDIR/summarize-ecdna-sizes.r \
            --in_dir_jba=$JABBA_DIR \
            --ktg_flag=$dir \
            --tn_file=$TNFILE


## Look at how ecDNA features change pre/post-chemo 
Rscript $SRCDIR/summarize-ecdna-by-patient-metadata.r \
            --in_dir_jba=$JABBA_DIR \
            --ktg_flag=$dir \
            --cgc_bed=$CGC_BED \
            --tn_file=$TNFILE \
            --metadata=$METADATA \
            --out_file=$WORKING_DIR/output/ecDNA-metadata-breakdown.txt
    
    
## Extract SNVs by VAF
Rscript $SRCDIR/extract-ecdna-mutations-by-vaf.r \
            --in_dir_ktg=$WORKING_DIR/output/clustered/kataegis-annotated \
            --in_dir_nc=$WORKING_DIR/output/nonClustered/annotated \
            --tn_file=$TNFILE \
            --out_dir=$VAF_DIR


## TC: deconstructSigs for the ecDNA mutations
Rscript run_deconstructSigs_ecDNA.R -d $VAF_DIR --highconf -o $output_sbs

