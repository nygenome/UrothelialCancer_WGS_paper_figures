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

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('-o','--output_prefix', help='output_prefix', dest='output_prefix')
parser.add_argument('-r','--ref', dest='ref', default="GRCh38", help='reference')
parser.add_argument('-i','--input_dir', dest='input_dir', help='input dir')

args = parser.parse_args()

matrices = matGen.SigProfilerMatrixGeneratorFunc(args.output_prefix, 
                                                 args.ref, 
                                                 args.input_dir,
                                                 plot=False, 
                                                 exome=False, 
                                                 bed_file=None, 
                                                 chrom_based=False, 
                                                 tsb_stat=False, 
                                                 seqInfo=False, 
                                                 cushion=100)
