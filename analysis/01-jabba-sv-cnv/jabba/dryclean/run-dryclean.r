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
## Run dryclean on a given tumor sample
libs = c('optparse', 'GenomicRanges', 'parallel', 'dryclean', 'data.table', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


## Collect arguments
option_list = list(
  make_option(c("-c", "--coverage"),      type='character',                   help="2 column (sample,path), headerless, tab-delimited list of normal samples and paths to their fragCounter output"),
  make_option(c("-d", "--detergent"),     type='character',                   help="Path to detergent.rds"),
  make_option(c("-o", "--out_file"),      type='character',                   help="Output file"),
  make_option(c("-f", "--germline_file"), type='character',                   help="Germline file used for filtering when germline_mode=FALSE"),
  make_option(c("-g", "--germline_mode"), action='store_true', default=FALSE, help="Run in germline event detection mode"))
opt = parse_args(OptionParser(option_list=option_list))



cov = readRDS(opt$coverage)



## Run differently depending on mode
if (opt$germline_mode) {
  
  if (all(is.na(cov$reads.corrected))) {
    warning('It appears fragcount multiscale coverage correction failed. Using original correction results...')
    field = 'reads.corrected.og'
  } else {
    field = 'reads.corrected'
  }
  
  message(paste0('Running germline event detection for ',basename(dirname(opt$coverage))))
  res = start_wash_cycle(cov=cov,
                         detergent.pon.path=opt$detergent,
                         whole_genome=T,
                         chr=NA,
                         field=field,
                         germline.filter=FALSE)
} else {
  message(paste0('Running dryclean in somatic mode for ',basename(dirname(opt$coverage))))
  
  ## Use multiscale-corrected coverage if possible, otherwise just use original correction
  if (all(is.na(cov$reads.corrected))) {
    warning('It appears fragcount multiscale coverage correction failed. Using original correction results...')
    field = 'reads.corrected.og'
  } else {
    field = 'reads.corrected'
  }
  
  ## Omit the "germline" filter for now -- it's agressive
  ## and removes some true somatic events
  if (!is.null(opt$germline_file) && file.exists(opt$germline_file)) {
  
    res = start_wash_cycle(cov=cov,
                           detergent.pon.path=opt$detergent,
                           whole_genome=T,
                           chr=NA,
                           germline.filter=TRUE,
                           germline.file=opt$germline_file)
  
  } else {
    
    res = start_wash_cycle(cov=cov,
                           detergent.pon.path=opt$detergent,
                           whole_genome=T,
                           chr=NA,
                           field=field,
                           germline.filter=FALSE)
    
  }
}

seqlevelsStyle(res) = 'UCSC'
saveRDS(res, opt$out_file)
