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
## Merge jabba junction and event count summaries
libs = c('optparse')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Get arguments
option_list = list(
  make_option(c("-j", "--jabba_dir"),  type='character', help="Jabba results directory"),
  make_option(c("-t", "--tn_file"),    type='character', help="Tumor-normal pairing file"),
  make_option(c("-p", "--prefix"),     type='character', help="File prefix (everything before '.summary.rds'", default='jabba.events'),
  make_option(c("-o", "--out_dir"),    type='character', help="Output directory"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor-normal pairs
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', comment.char='#')
tn$pair_name = paste0(tn$V1,'--',tn$V2)

## Merge junction summary files
summary.files = paste0(opt$jabba_dir,'/',tn$pair_name,'/',opt$prefix,'.junction.summary.rds')
missing.files = which(!file.exists(summary.files))

if (length(missing.files) > 0) {
  warning('The following files are missing: ', paste(summary.files[missing.files], collapse=', '))
  summary.files = summary.files[-missing.files]
}

res = do.call(rbind, lapply(summary.files, readRDS))
res = as.data.frame(res, stringsAsFactors=F)

if (length(missing.files) > 0){ 
  res$sample = tn$pair_name[setdiff(tn$pair_name, missing.files)]
} else {
  res$sample = tn$pair_name
}

res = res[,c(ncol(res), 1:(ncol(res)-1))]

out.file = paste0(opt$out_dir,'/',opt$prefix,'-junction-summary.csv')
write.csv(res, out.file, row.names=F, quote=F)



## Merge event summary.files
summary.files = paste0(opt$jabba_dir,'/',tn$pair_name,'/',opt$prefix,'.summary.rds')
missing.files = which(!file.exists(summary.files))

if (length(missing.files) > 0) {
  warning('The following files are missing: ', paste(summary.files[missing.files], collapse=', '))
  summary.files = summary.files[-missing.files]
}

res = do.call(rbind, lapply(summary.files, readRDS))
res = as.data.frame(res, stringsAsFactors=F)

if (length(missing.files) > 0){ 
  res$sample = tn$pair_name[setdiff(tn$pair_name, missing.files)]
} else {
  res$sample = tn$pair_name
}

res = res[,c(ncol(res), 1:(ncol(res)-1))]

out.file = paste0(opt$out_dir,'/',opt$prefix,'-event-summary.csv')
write.csv(res, out.file, row.names=F, quote=F)
