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


## Figure 2A - Junction burden heatmap
Rscript $SRCDIR/plt-2a-junction-burden-heatmap.r \
    --in_file=$JABBA_JUNCTION_BURDEN \
    --pp=$PURITY_PLOIDY \
    --metadata=$METADATA \
    --fga=$FGA \
    --tp53=$TP53_SOMATIC_STATUS \
    --out_file=$SRCDIR/fig2a-junction-burden-heatmap.svg


## Figure 2C - ecDNA/Kataegis co-localization barplot
Rscript $SRCDIR/plt-fig2c-ecdna-overlap-barplot.r \
    --in_dir_jba=$JABBA_DIR \
    --in_dir_ktg=$KATAEGIS_DIR \
    --ktg_flag=vcf_hc_vaf \
    --tn_file=$TNFILE \
    --out_file=$SRCDIR/fig2c-ecdna-overlap-barplot.svg


## Figure 2D - Kataegis distance to nearest junction ECDF
Rscript $SRCDIR/plt-fig2d-distance-to-junction-ecdf.r \
    --in_dir_ktg=$KATAEGIS_DIR \
    --in_dir_nc=$NONCLUSTERED_DIR \
    --tn_file=$TNFILE \
    --out_file=$SRCDIR/fig2d-distance-to-junction-ecdf.svg


## Figure 2E - Kyklonas/chemo violin plots
Rscript $SRCDIR/plt-fig2e-plot-kyklonas-chemo-violin.r \
    --in_dir_ktg=$KATAEGIS_DIR \
    --in_dir_nc=$NONCLUSTERED_DIR \
    --tn_file=$TNFILE \
    --out_file=$SRCDIR/fig2e-plot-kyklonas-chemo-violin.svg


## Figure 2F - ecDNA signatures barplot
Rscript $SRCDIR/plt-fig2f-ecdna-sbs-barplot.r \
    --in_file=$ECDNA_SBS_BY_VAF \
    --metadata=$METADATA \
    --out_file=$SRCDIR/fig2f-ecdna-sbs-barplot.svg


## Figure 2G - ecDNA / kyklonas / non-clustered genomic view
Rscript $SRCDIR/plt-fig2g-ecdna-kataegis-gtracks.r \
    --in_dir_jba=$JABBA_DIR \
    --in_dir_kataegis=$KATAEGIS_DIR \
    --ktg_flag=vcf_hc_vaf \
    --in_dir_nc=$NONCLUSTERED_DIR \
    --tn_file=$TNFILE \
    --out_file=$SRCDIR/fig2g-ecdna-kataegis-gtracks.pdf