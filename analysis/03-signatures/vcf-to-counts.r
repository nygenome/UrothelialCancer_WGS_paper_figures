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
## Tally counts of mutations per signature/timing combination
libs = c('optparse', 'VariantAnnotation','reshape2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



readMtSigVcf = function(f) { 
  
  x = readVcf(f)
  x = info(x)[, c('CLS','assigned_signature')]
  x$assigned_signature[is.na(x$assigned_signature)] = 'Unassigned'
  x = as.data.frame(x[!is.na(x$CLS), c('CLS','assigned_signature')])
  
  return(x)
    
}



sigs2FoldChange = function(x, sbs, cls=c('subclonal','clonal [late]','clonal [NA]','clonal [early]')) {
  
  x$assigned_signature = factor(x$assigned_signature, levels=sbs)
  x$CLS = factor(x$CLS, levels=cls)
  
  ## Reformat
  x = as.data.frame.matrix(table(x$assigned_signature, x$CLS)) 
  x$signature = rownames(x)
  
  x = x[c('signature', cls)]
  rownames(x) = NULL
  
  return(x)
  
}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="Input VCF"),
  make_option(c("-r", "--sbs_ref"),  type='character', help="COSMIC SBS reference"),
  make_option(c("-o", "--out_file"), type='character', help="Output CSV"))
opt = parse_args(OptionParser(option_list=option_list))



## Read SBS reference
sbs = colnames(read.csv(opt$sbs_ref, h=T, stringsAsFactors=F, sep='\t'))[-1]

## Read data, tabulate mutations per signature/clonality
res = sigs2FoldChange(readMtSigVcf(opt$in_file), sbs=sbs)

## Write result
write.csv(res, opt$out_file, row.names=F, quote=F)
message(opt$out_file)
