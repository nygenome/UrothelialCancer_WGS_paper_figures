#!/nfs/sw/python/python-2.7.8/bin/python
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: Minita Shah

################################################################# /COPYRIGHT ###
################################################################################

import os
import sys
import argparse
import re
import pandas as pd
import numpy as np
import gzip
import io

########## Classes and Functions ##########

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))

################### Main ###################
parser = ArgumentParser(prog='calc_vaf.py', description='', epilog='')
parser.add_argument('-i', '--inp-file', help='Input directory', required=True)
parser.add_argument('-o', '--out-file', help='Output Directory', required=True)
parser.add_argument('-p', '--positions-file', help='Positions file', required=True)

args = parser.parse_args()
inp_file = args.inp_file
out_file = args.out_file
positions_file = args.positions_file

df_pos = pd.read_csv(positions_file, header=None, sep="\t")
df_pos.columns = ['CHROM', 'POS', 'REF', 'ALT']
df_pos['CHROM'] = df_pos.CHROM.astype(str)
df_pos['POS'] = df_pos.POS.astype(str)
print(df_pos.head())

df_allele = pd.read_csv(inp_file, header=0, sep="\t")
df_allele.columns = ['CHROM', 'POS', 'REF', 'RD', 'A', 'C', 'G', 'T', 'STRAND']
df_allele['CHROM'] = df_allele.CHROM.astype(str)
df_allele['POS'] = df_allele.POS.astype(str)

print(df_allele.head())

# Merge het vcf and sample het alleles 
df_combined = pd.merge(df_pos, df_allele, on=['CHROM', 'POS'])
print(df_combined.head())

# Get ref, alt allele count
for i in ['A', 'C', 'G', 'T']:
    df_combined.loc[df_combined['REF_x']==i,'REF_COUNT'] = df_combined[i].astype(float)
    df_combined.loc[df_combined['ALT']==i,'ALT_COUNT'] = df_combined[i].astype(float)

# Calc BAF
df_combined.loc[(df_combined['ALT_COUNT']==0) & (df_combined['REF_COUNT']==0), 'VAF'] = 0.0
df_combined.loc[(df_combined['ALT_COUNT']!=0) | (df_combined['REF_COUNT']!=0), 'VAF'] = df_combined['ALT_COUNT']/(df_combined['ALT_COUNT'] + df_combined['REF_COUNT'])

# Print output
df_combined.rename(columns={'REF_x':'REF'}, inplace=True)
#out_file = os.path.join(out_dir, sample_name+".vaf.tobacco.txt")
df_combined.to_csv(out_file, sep='\t', index=False, columns=['CHROM', 'POS', 'REF', 'ALT', 'REF_COUNT', 'ALT_COUNT', 'VAF'])
