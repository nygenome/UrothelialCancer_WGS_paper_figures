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
## Intersect event footprints witha gene list
libs = c('optparse', 'reshape2', 'stringr', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



footprintsToGRanges = function(x, id) {
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
    
  return(res) 
}


collapseGenes = function(x) {
  
  x = x[!is.na(x)]
  
  if (length(x) == 0) {
    genes = ''
  } else {
    genes = unique(unlist(sapply(x, function(y) unlist(strsplit(y,',')))))
    genes = paste(trimws(genes), collapse=',')
  }    
  
  return(genes)
    
}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),   type='character', help="Input footprint RDS (jabba.events.footrpints.rds)"),
  make_option(c("-g", "--gene_bed"),  type='character', help="BED file with gene coordinates"),
  make_option(c("-o", "--out_file"),  type='character', help="Output BED"))
opt = parse_args(OptionParser(option_list=option_list))



## Read events, convert to granges
## columns: type, footprint, id
events = readRDS(opt$in_file)
events$type = as.character(events$type)
events$id = make.unique(events$type)

## If there aren't any events just touch file and quit
if (nrow(events) == 0) {
  
  system(paste('touch',opt$out_file))  
  warning('No events found in ', opt$in_file)
  
} else {
  
  events.gr = do.call(c, apply(events, 1, function(i) footprintsToGRanges(x=i[2], id=i[3])))
  
  ## Read genes
  genes = read.csv(opt$gene_bed, h=F, stringsAsFactors=F, sep='\t')
  colnames(genes) = c('chr','start','end','gene')
  genes = makeGRangesFromDataFrame(genes, keep.extra.columns=T)
  
  
  ## Aggregate genes across ranges in vcf 
  events.gr = events.gr %$% genes
  
  ## Collapes genes across complex events
  events.gene.list = tapply(events.gr$gene, events.gr$id, collapseGenes)
  events$gene = events.gene.list[events$id]
  events = events[,c('type','footprint','gene')]
  
  ## Write to disk
  write.table(events, opt$out_file, col.names=T, row.names=F, sep='\t', quote=T)
  message(opt$out_file)

}
