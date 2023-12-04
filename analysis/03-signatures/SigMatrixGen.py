#!/nfs/sw/python/python-3.6.1/bin/python

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('-o','--output_prefix', help='output_prefix', dest='output_prefix')
parser.add_argument('-r','--ref', dest='ref', default="GRCh38",
                    help='reference')
parser.add_argument('-i','--input_dir', dest='input_dir',
                    help='input dir')

args = parser.parse_args()

matrices = matGen.SigProfilerMatrixGeneratorFunc(args.output_prefix, args.ref, args.input_dir,plot=False, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
