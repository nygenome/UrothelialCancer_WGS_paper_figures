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

## This bash script contains figure 1 associated scripts
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 1A - Phylogenetic tree
## TODO: Tim will fill this in


## Figure 1B - Signature clonality/timing fold change boxplot
Rscript $SRCDIR/plt-fig1b-signature-clonality-bxp.r \
    --in_dir=$MUTATIONTIMER_DIR \
    --tn_file=$TNFILE \
    --metadata=$METADATA \
    --timing_summary=$TIMING_SUMMARY \
    --out_file=$SRCDIR/fig1b-signature-clonality-bxp.pdf


## Figure 1C - Mutations per month of exposure boxplot
## TODO: Tim will fill this in


## Figure 1D - Normal urothelium alluvial plot 1
## TODO: Tim will fill this in


## Figure 1E - Normal urothelium alluvial plot 2
## TODO: Tim will fill this in


## Figure 1F - dN/dS clonal selection scatter plot
Rscript $SRCDIR/plt-fig1f-dnds-gene-scatter.r \
    --tn_file=$TNFILE \
    --sig_genes_file=$DNDSCV_SIG_GENES \
    --in_dir_sig=$SIGNATURE_ASSIGNMENT_DIR \
    --in_dir_union=$UNION_VARIANT_DIR \
    --in_dir_singleton=$SINGLETON_DIR \
    --out_file=$SRCDIR/fig1f-dnds-gene-scatter.r

