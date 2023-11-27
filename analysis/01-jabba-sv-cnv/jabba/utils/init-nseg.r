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
## Create male- and female-specific GRanges objects specifying the normal copy number 
## for each segment (chr1-22,X,Y). This is passed to jabba's --nseg argument.
libs = c('optparse', 'GenomicRanges')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)

## Get arguments
option_list = list( 
  make_option(c("-c", "--chrlen"), type='character', help="Headerless tab-delmited file of chromosome lengths with columns chr start end (comment char '#' allowed"))
opt = parse_args(OptionParser(option_list=option_list))


## Read and convert to GRanges
x = read.csv(opt$chrlen, h=F, stringsAsFactors=F, sep='\t', comment.char='#')
colnames(x) = c('chr','start','end')
x = makeGRangesFromDataFrame(x)


## Male nseg 
male = x
male$cn = ifelse(seqnames(x) %in% c('chrX','chrY'), 1, 2)
male$ncn = male$cn

out.file = paste0(dirname(opt$chrlen),'/male_nseg.rds')
saveRDS(male, out.file)


## Female nseg
female = x
female$cn = ifelse(seqnames(x) %in% c('chrY'), 0, 2)
female$ncn = female$cn

out.file = paste0(dirname(opt$chrlen),'/female_nseg.rds')
saveRDS(female, out.file)
