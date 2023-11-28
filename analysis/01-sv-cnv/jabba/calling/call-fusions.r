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
## Call fusions on a given junction-balanced genome graph
libs = c('optparse', 'gUtils', 'gGnome', 'Matrix', 'parallel')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


## Get arguments
option_list = list(
  make_option(c("-j", "--jabba_gg"),    type='character', help="JaBbA simplified genome graph"),
  make_option(c("-l", "--gencode_gtf"), type='character', help="Gencode GTF"),
  make_option(c("-o", "--out_file"),    type='character', help="Output RDS file"))
opt = parse_args(OptionParser(option_list=option_list))


## Read jabba result
gg = readRDS(opt$jabba_gg)

## Call fusions
fus = fusions(graph=gg, 
              gencode=opt$gencode_gtf, 
              verbose=T)

## Only filter/write to disk if we have at least one fusion event
if(length(fus) > 0) {
  ## "fusions will output many "near duplicates" which just represent various combinations
  ## of near equivalent transcripts, we can filter these down using gWalk operations"
  fus = fus[!duplicated(genes)]
  
  if (length(fus) > 0) {
    saveRDS(fus, opt$out_file)
  }
}


