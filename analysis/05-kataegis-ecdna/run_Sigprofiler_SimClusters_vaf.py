#!/nfs/sw/python/python-3.6.1/bin/python
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: Timothy R. Chu 

################################################################# /COPYRIGHT ###
################################################################################

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerSimulator import SigProfilerSimulator as sigSim
from SigProfilerClusters import SigProfilerClusters
import sys

sigSim.SigProfilerSimulator(sys.argv[1], 
    sys.argv[2], 
    "GRCh38", 
    contexts = ['288'], 
    chrom_based=True, 
    simulations=100
    )
SigProfilerClusters.analysis(sys.argv[1], 
    "GRCh38", 
    "96", 
    ["288"], 
    sys.argv[2], 
    analysis="all", 
    sortSims=True, 
    subClassify=True, 
    correction=True, 
    calculateIMD=True, 
    max_cpu = 1,
    TCGA=False, 
    sanger=True,
    standardVC = False,
    includedVAFs=True,
    includedCCFs=False
)
