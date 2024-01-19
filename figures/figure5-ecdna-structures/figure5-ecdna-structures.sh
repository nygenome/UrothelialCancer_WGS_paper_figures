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


## Figure 5X - Circos plot of assembly
## TODO: Get code from Bishoy's team
 Rscript $SRCDIR/plt-assembly-circos.r \
        --flye_dir=$WORKING_DIR/flye \
        --ont_cov=$ont_cov \
        --sample=$ilmn_tumor \
        --regions=$jba_fp \
        --gtf=gencode.composite.collapsed.rds \
        --promoters=$PROMOTERS \
        --enhancers=$ENHANCERS \
        --superenhancers=$SUPERENHANCERS \
        --out_file=$ont_sample--$ilmn_tumor.circos.pdf


## Figure 5X - Scatterplot of junction support
Rscript $SRCDIR/plt-junction-read-length-distribution.r \
    --in_file=$ont_sample--$ilmn_tumor.pcr_free.readlength_distribution.txt \
    --out_file=$ont_sample--$ilmn_tumor.pcr_free.readlength_distribution.svg
