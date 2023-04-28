#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper, Rahul Singh, and Heather Geiger

################################################################# /COPYRIGHT ###
################################################################################

set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 4A - Comparison with TCGA data 
## TODO: Get code from Bishoy's team


## Figure 4B - CNV heatmap
Rscript $SRCDIR/plt-fig4b-cnv-heatmap.r \
    --in_dir=$JABBA_DIR \
    --tn_file=$TNFILE \
    --pp=$PURITY_PLOIDY \
    --metadata=$METADATA \
    --bed=$CNV_PLOT_BED \
    --gene_status=$RB1_SOMATIC_STATUS \
    --amp_status=$CCND1_AMP_STATUS \
    --out_file=$SRCDIR/fig4b-cnv-heatmap.svg


## Figure 4C - Pre/post-chemo comparison of JCN
Rscript $SRCDIR/plt-fig4c-ecdna-jcn-comparison-bxp.r \
    --in_file=$ECDNA_METADATA_SUMMARY \
    --metadata=$METADATA \
    --out_file=fig4c-ecdna-jcn-comparison-bxp.svg


## Figure 4D - Pre-chemo CCND1 amplification
tn=PM2079-Z1-1-Case-WGS--PM2079-EBC1-Ctrl-WGS
tumor=PM2079-Z1-1-Case-WGS
Rscript $SRCDIR/plt-fig4-ccnd1-gtrack.r \
    --jabba_gg=$JABBA_DIR/jabba/$tn/jabba.events.rds \
    --coverage=$JABBA_DIR/dryclean/tumor-decomposition/$tumor.dryclean-decomp.rds \
    --bed=$CCND1_BED \
    --aa=$AADIR/$tn.amplicons.10kb_4X_downsample10.bed \
    --padding=5E6 \
    --out_file=fig4d-ccnd1-gtrack.svg


## Figure 4E - Post-chemo CCND1 amplification
tn=PM2079-Z2-1-Case-WGS--PM2079-EBC1-Ctrl-WGS
tumor=PM2079-Z2-1-Case-WGS
Rscript $SRCDIR/plt-fig4-ccnd1-gtrack.r \
    --jabba_gg=$JABBA_DIR/jabba/$tn/jabba.events.rds \
    --coverage=$JABBA_DIR/dryclean/tumor-decomposition/$tumor.dryclean-decomp.rds \
    --bed=$CCND1_BED \
    --aa=$AADIR/$tn.amplicons.10kb_4X_downsample10.bed \
    --padding=5E6 \
    --out_file=fig4e-ccnd1-gtrack.svg


## Figure 4F - Pre-chemo CCND1 amplification
tn=PM911-Z5-1-Case-WGS--PM911-EBC2-1-Ctrl-WGS
tumor=PM911-Z5-1-Case-WGS
Rscript $SRCDIR/plt-fig4-ccnd1-gtrack.r \
    --jabba_gg=$JABBA_DIR/jabba/$tn/jabba.events.rds \
    --coverage=$JABBA_DIR/dryclean/tumor-decomposition/$tumor.dryclean-decomp.rds \
    --bed=$CCND1_BED \
    --aa=$AADIR/$tn.amplicons.10kb_4X_downsample10.bed \
    --padding=5E6 \
    --out_file=fig4f-ccnd1-gtrack.svg


## Figure 4G - Post-chemo CCND1 amplification
tn=PM911-Z1-1-Case-WGS--PM911-EBC2-1-Ctrl-WGS
tumor=PM911-Z1-1-Case-WGS
Rscript $SRCDIR/plt-fig4-ccnd1-gtrack.r \
    --jabba_gg=$JABBA_DIR/jabba/$tn/jabba.events.rds \
    --coverage=$JABBA_DIR/dryclean/tumor-decomposition/$tumor.dryclean-decomp.rds \
    --bed=$CCND1_BED \
    --aa=$AADIR/$tn.amplicons.10kb_4X_downsample10.bed \
    --padding=5E6 \
    --out_file=fig4g-ccnd1-gtrack.svg