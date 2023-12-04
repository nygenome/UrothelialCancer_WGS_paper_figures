#!/nfs/sw/python/python-2.7.11/bin/python
# coding=utf-8
########################################
# Minita Shah (mshah@nygenome.org)
# New York Genome Center
########################################

import os
import sys
import argparse
import re
import pandas as pd
import numpy as np
import gzip
import io

pd.set_option('display.max_colwidth', 1000)

########## Classes and Functions ##########

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.BytesIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str, normal_name: str, tumor_name: str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

################### Main ###################
parser = ArgumentParser(prog='get_tumor_pos_vafs.py', description='', epilog='')
parser.add_argument('-i', '--inp-vcf', help='Input vcf', required=True)
parser.add_argument('-o', '--out-prefix', help='Output file', required=True)
parser.add_argument('-t', '--tumor-name', help='Tumor name', required=True)
parser.add_argument('-n', '--normal-name', help='Normal name', required=True)

args = parser.parse_args()
inp_vcf = args.inp_vcf
out_prefix = args.out_prefix
tumor_name = args.tumor_name
normal_name = args.normal_name

df = read_vcf(inp_vcf)

signatures = ['SBS2','SBS13','SBS31','SBS35']
for i in signatures:
    df_tmp = df[df['INFO'].str.contains("assigned_signature="+i)].reset_index(drop=True)
    if df_tmp.shape[0] > 0:
        df_tmp[['AD', 'VAF', 'DP']] = df_tmp[tumor_name].str.split(':', expand=True)
        df_tmp['ASSIGNED_SIG'] = df_tmp['INFO'].str.split('[;=]').str[-1]

        df_tmp[['CHROM', 'POS', 'REF', 'ALT', 'VAF', 'ASSIGNED_SIG']].to_csv(out_prefix+"."+ i +".txt", sep="\t", index=False)

