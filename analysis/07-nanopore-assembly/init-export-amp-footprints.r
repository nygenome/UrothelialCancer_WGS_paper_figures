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
## Export footprints of interest to BED format
libs = c('optparse', 'GenomicRanges', 'reshape2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))


## Convert a jabba footprint record into GRanges
fp2gr = function(x, id, type) {
  
  x = gsub('\\+(?=,)|\\-(?=,)','',x, perl=T)
  x = gsub('\\+(?=;)|\\-(?=;)','',x, perl=T)
  x = gsub('\\+$|\\-$','',x, perl=T)
  x = strsplit(x,';|,')
  
  res = list()
  
  for (i in 1:length(x)) {
    
    seqnames = gsub(':.*','',x[[i]])
    start = as.numeric(sapply(x[[i]], function(y) unlist(strsplit(gsub('.*:','',y),'-'))[1]))
    end = as.numeric(sapply(x[[i]], function(y) unlist(strsplit(gsub('.*:','',y),'-'))[2]))
    
    res[[i]] = GRanges(seqnames, IRanges(start=start, end=end))
  }
  
  res = do.call(c, res)
  res$id = id
  res$type = type
  
  return(res) 
  
}



## Convert a footrpint table into a single GRanges
fp.table.to.gr = function(x) {

  x$id = make.unique(x$type, sep='_')

  res = mapply(FUN=fp2gr, 
               x=x$footprint, 
               id=x$id, 
               type=x$type, 
               SIMPLIFY=F)

  res = Reduce(c, res)

  return(res)

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="Path to jabba.events.genes.tab file"),
  make_option(c("-s", "--sample"),   type='character', help="Sample ID"),
  make_option(c("-p", "--padding"),  type='numeric',   help="Padding to add to intervals (default=0)", default=0),
  make_option(c("-o", "--out_dir"),  type='character', help="Output directory"))
opt = parse_args(OptionParser(option_list=option_list))


AMPS = c('dm', 'bfb', 'tyfonas')


## Read footprints 
fp = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')
fp = fp[fp$type %in% AMPS, ]
fp = fp.table.to.gr(fp) + opt$padding


## Export each one separately 
out.files = c()
for (i in unique(fp$id)) {

  fp.i = as.data.frame(fp[fp$id == i])[, 1:3]
  fp.i.out.file = paste0(opt$out_dir,'/',opt$sample,'.',i,'.bed')

  write.table(fp.i, fp.i.out.file, row.names=F, col.names=F, quote=F, sep='\t')
  out.files = c(out.files, fp.i.out.file)

}

out.file.list = paste0(opt$out_dir,'/',opt$sample,'.amplification_fp_list.txt')
writeLines(out.files, out.file.list)
cat(out.file.list)