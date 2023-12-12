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

## Attempt to harmonize JaBbA and AA calls for the purposes of visualization
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Collapse AA results
## AA decomposes amplicons, use each decomposition separately 
## PM636 Z2 failed the 10X/4CN/10kb run and is known to have a CCND1 
## ecDNA, use previous successful run
Rscript $SRCDIR/init-summarize-aa-results.r \
            --tn_file=$TNFILE \
            --amplicon_classification_profiles=$AA_CLASSIFICATIONS \
            --aa_run_flag=10kb_4X_downsample10 \
            --pm636_z2_profiles=$AA_CLASSIFICATIONS_PM636_Z2 \
            --pm636_z2_aa_run_flag=downsample10 \
            --aa_dir=$AA_DIR \
            --collapse_to_amplicons=FALSE \
            --out_file=$FIG_DIR/aa-results-summary.pdf \
            --out_file_txt=$WORKING_DIR/aa-results-summary.txt


## Compare AA results to jabba footprints
Rscript $SRCDIR/init-compare-aa-to-jabba.r \
            --tn_file=$TNFILE \
            --in_file_aa=$WORKING_DIR/aa-results-summary.txt \
            --in_dir_jabba=$JABBA_DIR \
            --out_file=$WORKING_DIR/aa-jabba-overlap-summary.txt


## Try calling ecDNA-like walks with the peel method
while read tumor normal gender patient; do 
    
    tn=$tumor--$normal
    walks=$WORKING_DIR/ecdna/$tn.amplicon.walks.rds
    walks_filtered=$WORKING_DIR/ecdna/$tn.amplicon.walks.filtered.rds

    Rscript $SRCDIR/init-decompose-amplicons.r \
                --sample=$tn \
                --in_file_comparison=$WORKING_DIR/aa-jabba-overlap-summary.txt \
                --in_file_aa=$WORKING_DIR/aa-results-summary.txt \
                --in_file_jabba=$JABBA_DIR/$tn/jabba.events.rds \
                --out_file=$walks

    [[ ! -f $walks ]] && continue

    Rscript $SRCDIR/init-filter-walks.r \
                --sample=$tn \
                --in_file_walks=$walks \
                --in_file_aa=$WORKING_DIR/aa-results-summary.txt \
                --out_file=$walks_filtered

done < $TNFILE



## Export "high-confidence" walks for visualization 
Rscript $SRCDIR/init-export-highconfidence-walks.r \
            --in_file_txt=$WORKING_DIR/putative-ecdna-calls.txt \
            --in_file_rds=$WORKING_DIR/putative-ecdna-calls.rds \
            --out_dir=$WORKING_DIR/ecdna/highconf_walks


## Export coverage to bedgraph
while read tumor normal gender patient; do 

    coverage=$JABBA_DIR/../dryclean/tumor-decomposition/$tumor.dryclean-decomp.rds

    Rscript $SRCDIR/dryclean2bedgraph.r \
                --in_file=$coverage \
                --yaml_template=$YAML_COV_TEMPLATE \
                --out_file_bedgraph=$WORKING_DIR/coverage/$tumor.dryclean-decomp.bedgraph \
                --out_file_yaml=$WORKING_DIR/coverage/$tumor.dryclean-decomp.yaml

done < $TNFILE



########################
## Plot with CycleViz ##
########################

. /gpfs/commons/groups/nygcfaculty/kancero/envs/miniconda3/etc/profile.d/conda.sh
conda activate cycleviz 

for f in $WORKING_DIR/ecdna/highconf_walks/*bed; do

    out_prefix=$FIG_DIR/highconf_reconstructions/$(basename $f .bed)
    
    tumor=$(basename $out_prefix | sed -e 's/.*intervals_//g' -e 's/--.*//g')
    coverage_bedgraph=$WORKING_DIR/coverage/$tumor.dryclean-decomp.bedgraph
    coverage_yaml=$WORKING_DIR/coverage/$tumor.dryclean-decomp.yaml
    
    [[ ! -f $coverage_bedgraph ]] && continue 

    ## Run CycleViz https://github.com/AmpliconSuite/CycleViz
    CycleViz.py \
        --ref hg38 \
        --structure_bed $f \
        --feature_yaml_list $FANTOM5_PROMOTER_YAML $FANTOM5_ENHANCER_YAML $SUPERENHANCER_YAML $coverage_yaml \
        --gene_highlight_list CCND1 \
        --rotate_to_min \
        --gene_subset_file "BUSHMAN" \
        --outname=$out_prefix \
        --center_hole 3 \
        --intertrack_spacing 0.2 \
        --feature_ref_offset 1 \
        --figure_size_style small

done 
