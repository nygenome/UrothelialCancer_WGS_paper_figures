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
## Report summary stats for ecDNA sizes
libs = c('optparse', 'GenomicRanges')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))



alpha = function(col, alpha) {
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}



## Read jabba footprints table
read.jabba.fp = function(f, sample_id) {

  if (!file.exists(f)) {
    return(NULL)
  }
  
  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  x$sample = sample_id

  return(x)

}



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



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir_jba"),   type='character', help="Jabba output directory"),
  make_option(c("-f", "--ktg_flag"),     type='character', help="Flag to set SigProfilerClusters run type"),
  make_option(c("-t", "--tn_file"),      type='character', help="Tab delimited tumor-normal pairing file"))
opt = parse_args(OptionParser(option_list=option_list))



## Read TN pairs and get input files
tn = read.csv(opt$tn, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)
tn$in_files_jba_fp = paste0(opt$in_dir_jba,'/',tn$pair_name,'/jabba.events.footprints.kataegis.',opt$ktg_flag,'.txt')

## Read input 
jba = do.call(rbind, mapply(read.jabba.fp, f=tn$in_files_jba_fp, sample_id=tn$tumor, SIMPLIFY=F, USE.NAMES=F))

## Select ecDNAs
jba = jba[jba$final_type == 'ecDNA', ]


## Get ecDNA sizes
sizes = sapply(jba$footprint, function(x) sum(width(fp2gr(x=x, id='ecdna'))))

message('Minimum size:', min(sizes)/1E6, ' MB')
message('Median size:', median(sizes)/1E6, ' MB')
message('Mean size:', mean(sizes)/1E6, ' MB')
message('Maximum size:', max(sizes)/1E6, ' MB')

quantile(sizes, c(0.25, 0.5, 0.75))