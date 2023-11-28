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
## Merge JaBbA junctions into a single file, deduplicating within patients to 
## avoid double-counting 
libs = c('optparse', 'gGnome', 'GenomicRanges', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



loadJunctions = function(f) {
  
  gg = readRDS(f)
  junc = gg$junctions[type == 'ALT'] 
  
  return(junc)
  
}



## Read a list of gGraphs and merge junctions
mergeJunctions = function(x, pad=0) {

  juncs = lapply(x, loadJunctions)
  juncs$pad = pad
  junc.union = do.call(merge, juncs)

  return(junc.union)

}



## Convert merged junction object to bed file
junc2bed = function(x) {

  x = as.data.frame(x$dt)
  res = NULL
  seen.by = grep('seen.by.', colnames(x), fixed=T)

  for (i in 1:nrow(x)) {

    ## If a junction isn't observed in a sample, columns
    ## are set to NA 
    suffix = which(unlist(x[i, seen.by]))[1]
    if (suffix == 1) {
      suffix = ''
    } else {
      suffix = paste0('.', suffix - 1)
    }

    res.i.1 = unname(unlist(x[i, paste0(c('chr1','start1','end1','str1', 'class'), suffix)]))
    res.i.2 = unname(unlist(x[i, paste0(c('chr2','start2','end2','str2', 'class'), suffix)]))

    res.i.1['id'] = i
    res.i.2['id'] = i

    res = rbind(res, res.i.1, res.i.2)

  }

  ## Reformat result
  res = as.data.frame(res)
  colnames(res) = c('chr', 'start', 'end', 'strand', 'name', 'id')

  res$start = as.numeric(res$start)
  res$start = as.numeric(res$end)

  return(res)

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir"),   type='character', help="Directory holding JaBbA results"),
  make_option(c("-t", "--tn"),       type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-o", "--out_file"), type='character', help="Figure output file"))
opt = parse_args(OptionParser(option_list=option_list))


## Read tumor normal pairs, get files 
tn = read.csv(opt$tn, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)
tn$in.files = paste0(opt$in_dir,'/',tn$pair_name,'/jabba.events.rds')

if (any(!file.exists(tn$in.files))) {
  print(tn$in.files[!file.exists(tn$in.files)])
  stop('Missing files!')
  
}


## Merge junctions within patients so that we don't double-count the same breakpoints twice
## Zero padding ensures only identical breakpoints are merged
juncs = tapply(tn$in.files, tn$patient, mergeJunctions)
juncs = lapply(juncs, junc2bed)
juncs = do.call(rbind, juncs)
juncs$patient = gsub('\\..*' ,'', rownames(juncs))


## Write result
write.table(juncs, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)
