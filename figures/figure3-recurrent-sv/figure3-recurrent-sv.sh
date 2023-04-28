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

set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 3A - FishHook manhattan plot 
Rscript $SRCDIR/plt-fig3a-manhattan-plot.r \
    --in_file=$FISHHOOK_HITS_ANNOTATED \
    --fdr=0.25 \
    --out_file=$SRCDIR/fig3a-manhattan-plot.svg


## Figure 3B - JaBbA oncoprint for FishHook hits
Rscript $SRCDIR/plt-fig3b-jabba-fishhook-oncoprint.r \
    --jba_dir=$JABBA_DIR \
    --tn_file=$TNFILE \
    --bed=$FISHHOOK_HITS_ANNOTATED_BED \
    --metadata=$METADATA \
    --out_file=$SRCDIR/fig3b-jabba-fishhook-oncoprint.svg