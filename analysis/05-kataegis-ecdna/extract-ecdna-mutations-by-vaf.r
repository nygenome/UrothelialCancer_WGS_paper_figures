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
## Extract ecDNA mutations into files, splitting by VAF
libs = c('optparse', 'GenomicRanges', 'data.table')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


## Read SigProfilerClusters annotated output
read.sigprofiler = function(f) {

  if (!file.exists(f)) {
    return(NULL)
  }

  x = fread(f, header=T, stringsAsFactors=F, sep='\t', data.table=F)

  x$end = x$pos
  x = makeGRangesFromDataFrame(x, 
                               seqnames.field='chrom',
                               start.field='pos',
                               end.field='end',
                               keep.extra.columns=T)

  return(x)

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir_ktg"),  type='character', help="Directory holding SigProfilerClusters VCFs"),
  make_option(c("-n", "--in_dir_nc"),   type='character', help="Directory holding SigProfilerClusters VCFs (unclustered)"),
  make_option(c("-t", "--tn_file"),     type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-o", "--out_dir"),     type='character', help="Output directory"))
opt = parse_args(OptionParser(option_list=option_list))


## Divide calls into three evenly-sized VAF bins
vaf.bins = c(0, 1/3, 2/3, 1)
vaf.bin.names = c('vaf_0_333', 'vaf_333_667', 'vaf_667_1')
mutation.count.minimum = 50

## Read TN pairs and get input files
tn = read.csv(opt$tn, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))

tn$pair_name = paste0(tn$tumor,'--',tn$normal)
tn$in_files_ktg = paste0(opt$in_dir_ktg,'/',tn$pair_name,'-mutationtimer-snv.annotated.jabba.txt')
tn$in_files_nc = paste0(opt$in_dir_nc,'/',tn$pair_name,'-mutationtimer-snv.annotated.jabba.txt')


## Read input 
ktg = lapply(tn$in_files_ktg, read.sigprofiler)
nc = lapply(tn$in_files_nc, read.sigprofiler)


## Select ecDNA mutations
ktg = lapply(ktg, function(x) x[x$ecdna & x$tricontext != 'DBS'])
nc = lapply(nc, function(x) x[x$jabba_event == 'ecDNA' & x$tricontext != 'DBS'])



## For each tumor 
for (i in 1:nrow(tn)) {


  ## Select columns and combine kataegis/non-clustered
  columns = c('seqnames','start', 'ref', 'alt', 'vaf')
  
  if (is.null(ktg[[i]])) {

    res = as.data.frame(nc[[i]])[, columns]

  } else {

    res = rbind(as.data.frame(ktg[[i]])[, columns], as.data.frame(nc[[i]])[, columns])

  }


  ## Split by VAF
  res$vaf_bin = cut(res$vaf, vaf.bins)
  res = split(res, res$vaf_bin)


  ## Write each bin separately
  for (j in 1:length(vaf.bin.names)) {

    out.file = paste0(opt$out_dir,'/',tn$pair_name[i],'.',vaf.bin.names[j],'.snv.txt')

    ## Skip if there's nothing to write out 
    if (is.null(res[[j]]) || nrow(res[[j]]) <= mutation.count.minimum) {
      next
    }

    ## Select columns and rename chr field
    res.vaf = res[[j]][, columns]
    colnames(res.vaf)[1] = 'chr'

    ## Write result
    write.table(res.vaf, out.file, row.names=F, col.names=T, quote=F, sep='\t')
    message(out.file)

  }

}
