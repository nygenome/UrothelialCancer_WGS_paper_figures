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
## Handle formatting GC%/mappability scores for fragcounter correction
libs = c('optparse', 'gUtils', 'GenomicRanges')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)

## Collect arguments
option_list = list(
  make_option(c("-s", "--score_file"),     type='character', help="3 column (chr,pos,score) or 4 column (chr,start,stop,score) input"),
  make_option(c("-c", "--chr_intervals"),  type='character', help="Chromosome intervals"),
  make_option(c("-w", "--window_size"),    type='numeric', help="Window size for tiling chromosomes with non-overlapping windows"),
  make_option(c("-r", "--convert_to_rds"), action="store_true", default=FALSE, help="Just convert to RDS"),  
  make_option(c("-o", "--out_file"),       type='character', help="Output RDS GenomicRanges object"))
opt = parse_args(OptionParser(option_list=option_list))



## Read in scores (i.e. mappability, gc)
scores = read.csv(opt$score_file, h=F, stringsAsFactors=F, sep='\t')

if (ncol(scores) == 3) {
  scores = GRanges(scores[,1], IRanges(scores[,2], scores[,2]), score=scores[,3])
} else if (ncol(scores) == 4) {
  scores = GRanges(scores[,1], IRanges(scores[,2], scores[,3]), score=scores[,4])
} else {
  stop('Input must be 3-column (chr,pos,score) or 4-column (chr,start,stop,score) format')
}



## Tile if we're doing more than just converting
if (!opt$convert_to_rds) {
  
  ## Read in chr intervals, generate non-overlapping windows
  chr = read.csv(opt$chr_intervals, h=F, stringsAsFactors=F, sep='\t')
  chr = GRanges(chr[,1], IRanges(chr[,2], chr[,3]))
  chr = gr.tile(chr, opt$window_size)
  
  
  ## Guess if scores are out of 1 or 100
  ## Convert to fraction of 1 if needed
  if (!all(sample(scores$score, 1000) <= 1)) {
    scores$score = scores$score / 100
  }
  
  
  ## Aggregate across tiled intervals using the mean, set NA to -1 
  chr = chr %$% scores
  chr$score[is.na(chr$score)] = -1 
  
  
  ## Coordinates may need to be shifted from 0-indexed to 1-indexed
  if (start(chr[1, ]) == 0) {
    chr = chr %+% 1
  }
    
  saveRDS(chr, file=opt$out_file)
  
} else {
  
  saveRDS(scores, file=opt$out_file)
  
}





