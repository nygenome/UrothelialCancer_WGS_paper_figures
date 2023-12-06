#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: Timothy R. Chu 

################################################################# /COPYRIGHT ###
################################################################################

libs = c('optparse', 'deconstructSigs')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)

option_list = list(
  make_option(c("-i", "--input"),   type="character", default=NULL,   help="[REQUIRED] ID83 input file",               metavar="character"),
  make_option(c("-o", "--output"),  type="character", default=NULL,   help="[REQUIRED] output",                        metavar="character"),
  make_option(c("-v", "--version"), type="character", default="v3.4", help="COSMIC version [v3.4,v3.2] default: v3.4", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


if(opt$version == "v3.4"){
  load("signatures.ID.cosmic.v3.4.oct2023.grch37.rda")
  sig_ref = signatures.ID.cosmic.v3.4.oct2023.grch37
}
else if(opt$version=="v3.2"){
  load("signatures.id.cosmic.v3.2.march2021.grch37.rda")
  sig_ref = signatures.id.cosmic.v3.2.march2021.grch37
}else{
  print("error. version must be 'v3.2' or 'v3.4")
}


id83.input = read.table(opt$input, header = T, check.names=F, stringsAsFactors = F)
rownames(id83.input) = id83.input$MutationType
id83.input = id83.input[,2:ncol(id83.input)]
id83.input = t(id83.input)
id83.input = as.data.frame(id83.input)

id83.output = data.frame()
for(samplename in rownames(id83.input)){
  idsigs = whichSignatures(tumor.ref=id83.input, signatures.ref=sig_ref, sample.id=samplename, contexts.needed=T, signature.cutoff = 0)
  id83.output = rbind(id83.output, idsigs$weights)
}

id83.output = cbind(Sample = rownames(id83.output), id83.output)
write.table(id83.output, file=opt$output, sep="\t", quote=FALSE, row.names=FALSE)
