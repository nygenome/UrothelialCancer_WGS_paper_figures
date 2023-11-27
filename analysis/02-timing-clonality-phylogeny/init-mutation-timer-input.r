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
## Extract allele-specific copy number estimates from jabba output 
libs = c('optparse', 'reshape2', 'GenomicRanges', 'gUtils', 'magrittr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Get arguments
option_list = list(
  make_option(c("-j", "--jba"),         type='character', help="Input jabba.simple.rds with slot agtrack (i.e., jabba was run with a het pileup supplied"),
  make_option(c("-p", "--purity"),      type='numeric',   help="Sample purity"),
  make_option(c("-t", "--tumor"),       type='character', help="Tumor ID"),
  make_option(c("-g", "--gender"),      type='character', help="Gender"),
  make_option(c("-o", "--out_file_cn"), type='character', help="Output RDS with per-segment allele-specific copy number"))
opt = parse_args(OptionParser(option_list=option_list))



#######################################
## Jabba allele-specific copy number ## 
#######################################

# Read data
jba = readRDS(opt$jba)
aseg = jba$asegstats
tseg = jba$segstats

aseg

## Reshape
seg = aseg %Q% (strand == '+') %>%
  as.data.table %>%
  dcast.data.table(seqnames + start + end ~ type, value.var='cn', fun.aggregate=sum) %>%
  setnames(c('high', 'low'), c('major_cn', 'minor_cn')) %>%
  makeGRangesFromDataFrame(keep.extra.columns=T)

seg

## Add sample purity
seg$clonal_frequency = opt$purity
tseg$clonal_frequency = opt$purity

## Kick out width=1 segments
seg = seg[width(seg) > 1]
tseg = tseg[width(tseg) > 1]




## No concept of major/minor copy number for male chrX, replace
## with mornal jabba segmentation
if (opt$gender == 'male') {

  message('Gender set to male -- using total copy number as major copy number for chrX')

  tseg = tseg[!duplicated(tseg$tile.id) & as.character(seqnames(tseg)) == 'chrX']
  strand(tseg) = '*'
  
  tseg$major_cn = tseg$cn
  tseg$minor_cn = 0
  mcols(tseg) = mcols(tseg)[, c('major_cn', 'minor_cn', 'clonal_frequency')]

  seg = seg[as.character(seqnames(seg)) != 'chrX']
  seg = c(seg, tseg)
  seqlengths(seg) = NA

}


## For now, kick out any segments that don't have allele-specific copy numbers assigned
# n.rm = sum(is.na(seg$major_cn) & is.na(seg$minor_cn))
# seg = seg[!is.na(seg$major_cn) & !is.na(seg$minor_cn)]
# message('Removed ', n.rm,' segments without allele-specific copy numbers')

## Write result
saveRDS(seg, opt$out_file_cn)
message(opt$out_file_cn)
