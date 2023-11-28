#!/nfs/sw/R/R-4.0.0/bin/Rscript
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
## Convert dryclean RDS to a bedfile for CycleViz plotting
libs = c('optparse', 'GenomicRanges', 'rtracklayer', 'yaml')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),           type='character', help="Drycleaned coverage RDS"),
  make_option(c("-y", "--yaml_template"),     type='character', help="YAML template"),
  make_option(c("-o", "--out_file_bedgraph"), type='character', help="Output bedgraph"),
  make_option(c("-O", "--out_file_yaml"),     type='character', help="Output YAML file"))
opt = parse_args(OptionParser(option_list=option_list))



## Read coverage 
dta = readRDS(opt$in_file)

## Only include bins with coverage 
dta = dta[!is.na(dta$foreground)]
dta$score = round(dta$foreground, 3)
mcols(dta) = mcols(dta)[,'score', drop=F]

## Write result
rtracklayer::export(dta, opt$out_file_bedgraph, format='bedgraph')
message(opt$out_file_bedgraph)


## Read YAML template 
yml = read_yaml(opt$yaml_template)
yml$primary_feature_bedgraph = opt$out_file_bedgraph
write_yaml(yml, opt$out_file_yaml)
message(opt$out_file_yaml)
