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
## Build list of eligible regions for FishHook analysis based on mappability and 
## the block list used for JaBbA
libs = c('optparse', 'GenomicRanges', 'gUtils', 'data.table')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


## Get arguments
option_list = list(
  make_option(c("-c", "--chr_len"),     type='character', help="Chromosome lengths"),
  make_option(c("-m", "--mappability"), type='character', help="Mappability BED"),
  make_option(c("-b", "--block"),       type='character', help="Block list"),
  make_option(c("-o", "--out_file"),    type='character', help="Output BED"))
opt = parse_args(OptionParser(option_list=option_list))



## Read chromosomes of interest
chr.len = read.csv(opt$chr_len, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','end'))
chr.len$start = 1
chr.len = chr.len[!chr.len$chr %in% c('chrM', 'chrY'), ]
chr.len = makeGRangesFromDataFrame(chr.len)


## Read mappability
map = data.table::fread(opt$mappability, 
                        header=F, 
                        stringsAsFactors=F, 
                        sep='\t', 
                        data.table=F, 
                        col.names=c('chr','start','end','map'))
map = makeGRangesFromDataFrame(map, keep.extra.columns=T)
map = map[map$map == 1 & as.character(seqnames(map)) %in% as.character(seqnames(chr.len))]


## Read block list
block = read.csv(opt$block, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','start','end'))
block = makeGRangesFromDataFrame(block)
block = block[as.character(seqnames(block)) %in% as.character(seqnames(chr.len))]


## First, remove blocked regions
allow = disjoin(grbind(chr.len, block))
allow = allow[!allow %^% block]

## Next, only take regions with mappability of 1 
allow = disjoin(grbind(allow, map))
allow = allow[allow %^% map]

## Report what we retained 
pct = round(sum(width(allow)) / sum(width(chr.len)), 4) * 100
message('\nRetaining ', pct,'% of chr1-22,X\n')

## Write result
res = data.frame(chr=as.character(seqnames(allow)), start=start(allow), end=end(allow))
write.table(res, opt$out_file, row.names=F, col.names=F, quote=F, sep='\t')
message(opt$out_file)
