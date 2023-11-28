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
## Export filter-pass walks with at least 80% reciprocal overlap 
libs = c('optparse', 'GenomicRanges', 'gGnome', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


RO_CUTOFF = 0.8


## Get arguments
option_list = list(
  make_option(c("-i", "--in_file_txt"),        type='character', help="Txt output of init-classify-ecdna-calls.r"),
  make_option(c("-I", "--in_file_rds"),        type='character', help="RDS output of init-classify-ecdna-calls.r"),
  make_option(c("-o", "--out_dir"),        type='character', help="Output directory"))
opt = parse_args(OptionParser(option_list=option_list))


## Read data
dta = read.csv(opt$in_file_txt, h=T, stringsAsFactors=F, sep='\t')
walks = readRDS(opt$in_file_rds)

### Select high confidence reconstructions
idx.sel = dta$jabba_walk_ro >= RO_CUTOFF
dta = dta[idx.sel, ]
walks = walks[idx.sel]


## For each reconstruction
out.files = c()

for (i in 1:length(walks)) {

  ## Get nodes, preserving orientations
  ## For now, collapse adjacent nodes if their orientation
  ## is the same 
  # res = granges(walks[[i]]$nodes$gr)
  res = reduce(granges(walks[[i]]$nodes$gr), ignore.strand=FALSE)
  res = data.frame(chr=as.character(seqnames(res)), 
                   start=start(res),
                   end=end(res),
                   strand=as.character(strand(res)),
                   connected_to_previous='True')

  ## Write result
  out.file = paste0(opt$out_dir, '/', walks[[i]]$dt$aa_events[1],'.bed')
  out.files = c(out.file, out.files)

  if(any(duplicated(out.files))) {
    stop('Name collision! A single AA event overlaps multiple jabba decompositions')
  }

  write.table(res, out.file, row.names=F, col.names=F, quote=F, sep='\t')
  message(out.file)

}
