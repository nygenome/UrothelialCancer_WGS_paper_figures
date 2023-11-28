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
## Annotate non-clustered SNVs with any associated jabba events
libs = c('optparse', 'GenomicRanges', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)

JBA_AMP = c('dm', 'cpxdm', 'bfb', 'tyfonas')

## Convert a footprint record into GRanges
fp2gr = function(x, id) {
  
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



## Read SigProfilerClusters annotated output into
## GRangesList with an entry for each event
read.sigprofiler = function(f) {

  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  x$end = x$pos
  x = makeGRangesFromDataFrame(x, 
                               seqnames.field='chrom',
                               start.field='pos',
                               end.field='end',
                               keep.extra.columns=T)

  x$jabba_event = 'none'

  return(x)

}



## Write sigprofiler results from a GRangesList
write.sigprofiler = function(x, f) {
  print(head(x))
  x = as.data.frame(x)
  print(head(x))

  x = x[!colnames(x) %in% c('end','width','strand')]
  colnames(x)[1:2] = c('chrom', 'pos')

  write.table(x, f, row.names=F, col.names=T, quote=F, sep='\t')

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file_nc"),   type='character', help="Annotated sigprofilerclusters kataegis results"),
  make_option(c("-I", "--in_file_jba"),  type='character', help="Jabba event footprints, annotated with kataegis and ecDNA status"),
  make_option(c("-o", "--out_file_nc"),  type='character', help="Output kataegis calls with associated jabba events"))
opt = parse_args(OptionParser(option_list=option_list))



## Read input files
if (!file.exists(opt$in_file_nc)) {

  message('No nonclustered events detected')
  quit(save='no', status=0)

} else {

  vcf = read.sigprofiler(opt$in_file_nc)

}



## Case where no jabba events were called for this sample
if (!file.exists(opt$in_file_jba)) {
  
  message('No jabba events detected')

  ## Write kataegis results if we have them
  if (!is.null(vcf)) {
    write.sigprofiler(vcf, opt$out_file_nc)
    message(opt$out_file_nc)
  }
  
  quit(save='no', status=0)

} else {

  jba = read.csv(opt$in_file_jba, h=T, stringsAsFactors=F, sep='\t')

}


## Read footprints into granges
fp = mapply(fp2gr, x=jba$footprint, id=jba$final_type, SIMPLIFY=F)
fp = Reduce(c, fp)

## Annotate each mutation with the jabba event it falls within
vcf$jabba_event = gr.val(vcf, fp, val='id', sep=',')$id
vcf$jabba_event[vcf$jabba_event == ''] = 'none'

## Write annotated kataegis results
write.sigprofiler(vcf, opt$out_file_nc)
message(opt$out_file_nc)
