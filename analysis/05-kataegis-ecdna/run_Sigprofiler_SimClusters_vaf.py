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
