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
## Compute a breakpoint dispersion score for FishHook hits
libs = c('optparse', 'GenomicRanges', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Read a text file with header and columns chr/start/end into GRanges
txt2gr = function(f) {

  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  x = makeGRangesFromDataFrame(x, keep.extra.columns=T)
  return(x)

}



## Compute dispersion score for a GRanges object 
dispersion = function(gr) {

  ## Adhere to one breakpoint per-patient per-locus
  gr = gr[!duplicated(gr$patient)]

  ## Nothing to compare if only 1 or 2 breakpoints observed in a bin
  if (length(gr) %in% 1:2) {
    return(NA)
  }

  dist = c()

  for (i in 1:(length(gr)-1)) {
    
    ## Once we've compared the ith junction to the rest, we 
    ## don't need to include it in future comparisons
    dist = c(dist, distance(gr[i], gr[-(1:i)], all=T, ignore.strand=T))

  }
  
  ## Calculate dispersion
  res = mad(dist) / median(dist)
  return(res)

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),               type='character', help="Output of annotate-fishhook-results.r"),
  make_option(c("-j", "--junctions"),             type='character', help="Merged jabba junctions"),
  make_option(c("-x", "--one_patient_per_locus"), type='logical',   help="Restrict to one breakpoint per-patient per-locus?", default=TRUE),
  make_option(c("-o", "--out_file"),              type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



## Read data
dta = txt2gr(opt$in_file)
dta$tile.id = paste0('tile.',dta$tile.id)
junc = txt2gr(opt$junctions)


## Split junctions by bin
junc = gr.val(junc, dta, 'tile.id')
junc = split(junc, junc$tile.id)


## Get dispersion
res.disp = sapply(junc, dispersion)
dta$dispersion = res.disp[dta$tile.id]


## Get junction count 
if (opt$one_patient_per_locus){
  message('Deduplicating')
  res.njunc = sapply(junc, function(x) length(unique(x$patient)))
} else {
  res.njunc = sapply(junc, length)
}

dta$junctions = res.njunc[dta$tile.id]
dta$junctions[is.na(dta$junctions)] = 0


## Write result
dta = as.data.frame(dta)
colnames(dta)[1] = 'chr'
write.table(dta, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)
