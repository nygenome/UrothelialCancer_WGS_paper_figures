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
## Annotate jabba events with presence/absence of kataegis, and vice versa
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

  if (!file.exists(f)) {
    return(NULL)
  }

  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  x$end = x$pos
  x = makeGRangesFromDataFrame(x, 
                               seqnames.field='chrom',
                               start.field='pos',
                               end.field='end',
                               keep.extra.columns=T)

  x$jabba_event = 'none'
  x$ecdna = FALSE
  x$distance_to_nearest_ecdna = NA

  x = split(x, x$kataegis_id)
  return(x)

}



## Write sigprofiler results from a GRangesList
write.sigprofiler = function(x, f) {

  x = as.data.frame(Reduce(c, x))
  x = x[!colnames(x) %in% c('end','width','strand')]
  colnames(x)[1:2] = c('chrom', 'pos')

  write.table(x, f, row.names=F, col.names=T, quote=F, sep='\t')

}



## Read amplicon architect results. If there weren't any 
## amplicons, just initialize an empty GRanges since we 
## don't have to write these results out 
read.aa = function(f) {

  if (!file.exists(f)) {

    x = GRanges()
    mcols(x)[, c('aa_amplicon', 'aa_id', 'aa_class')] = character()

  } else {

    x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
    x = makeGRangesFromDataFrame(x, keep.extra.columns=T)
    x = x[x$aa_class == 'Cyclic']

  }

  return(x)

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file_ktg"),  type='character', help="Annotated sigprofilerclusters kataegis results"),
  make_option(c("-I", "--in_file_jba"),  type='character', help="Jabba event footprints"),
  make_option(c("-a", "--in_file_aa"),   type='character', help="AmpliconArchitect event bed files"),
  make_option(c("-o", "--out_file_ktg"), type='character', help="Output kataegis calls with associated jabba events"),
  make_option(c("-O", "--out_file_jba"), type='character', help="Output jabba footprints with associated kataegis calls"))
opt = parse_args(OptionParser(option_list=option_list))



## Read input files
ktg = read.sigprofiler(opt$in_file_ktg)
jba = read.csv(opt$in_file_jba, h=T, stringsAsFactors=F, sep='\t')
aa = read.aa(opt$in_file_aa) ## read.aa filters for cyclic amplicons


## Case where no jabba events were called for this sample
if (nrow(jba) == 0) {
  
  message('No jabba events detected')

  ## Write kataegis results if we have them
  if (!is.null(ktg)) {
    write.sigprofiler(ktg, opt$out_file_ktg)
    message(opt$out_file_ktg)
  }
  
  quit(save='no', status=0)

}



## Read footprints into granges
jba$id = make.unique(jba$type)
jba[, c('has_kataegis', 'kataegis_id', 'kataegis_mutations', 'apobec_associated', 
        'aa_cyclic', 'aa_intersection_fp' , 'aa_intersection_fp_bp', 'final_type')] = NA
fp = mapply(fp2gr, x=jba$footprint, id=jba$id, SIMPLIFY=F)
names(fp) = jba$id

## Case where no kataegis events were called for this sample
if (is.null(ktg)) {

  message('No kataegis events detected')
  jba$final_type = jba$type
  jba$has_kataegis = FALSE
  
  write.table(jba, opt$out_file_jba, row.names=F, col.names=T, quote=F, sep='\t')
  message(opt$out_file_jba)
  quit(save='no', status=0)
  
}


## For each event
for (i in 1:length(fp)) {

  ## Check to see if this jabba event overlaps 
  ## a cyclic AmpliconArchitect call
  jba$aa_cyclic[i] = any(fp[[i]] %^% aa)


  ## If this is a AA cyclic AND jabba bfb/tyfonas/dm,
  ## call it ecdna
  if (jba$aa_cyclic[i] && jba$type[i] %in% JBA_AMP) {
    jba$final_type[i] = 'ecDNA'
  } else {
    jba$final_type[i] = jba$type[i]
  }


  ## Annotate with intersection between AA and jabba
  if(jba$final_type[i] == 'ecDNA') {

    ## Retain intersection of AA and jabba intervals 
    gr.union = disjoin(grbind(fp[[i]], aa))
    gr.union = gr.union[(gr.union %^% fp[[i]]) & (gr.union %^% aa)]

    jba$aa_intersection_fp[i] = paste(gr.string(gr.union), collapse=',')
    jba$aa_intersection_fp_bp[i] = sum(width(gr.union))

  }


  ## Check to see if any of the kataegis events fall completely
  ## within this event
  ktg.in.fp = grl.in(ktg, fp[[i]], some=FALSE, only=FALSE, logical=TRUE, ignore.strand=TRUE)
  jba$has_kataegis[i] = any(ktg.in.fp)

  ## If so,
  if (any(ktg.in.fp)) {
    
    ## Annotate jabba event table
    jba$kataegis_id[i] = paste(names(ktg)[ktg.in.fp], collapse=',')
    jba$kataegis_mutations[i] = paste(sapply(ktg[ktg.in.fp], length), collapse=',')
    jba$apobec_associated[i] = paste(sapply(ktg[ktg.in.fp], function(x) any(x$apobec_associated)), collapse=',')

    ## Annotate kataegis event
    for (j in which(ktg.in.fp)) {
      ktg[[j]]$jabba_event = jba$id[i]
      ktg[[j]]$ecdna = jba$final_type[i] == 'ecDNA'
    }

  } 

}



## Compute distance to nearest ecDNA event 
for (i in 1:length(fp)) {
  fp[[i]]$type = jba$final_type[i]
}

fp = Reduce(c, fp)
fp = fp[fp$type == 'ecDNA']

for (i in 1:length(ktg)) {

  dist.to.ecdna = distanceToNearest(ktg[[i]], fp)
  ktg[[i]]$distance_to_nearest_ecdna[queryHits(dist.to.ecdna)] = mcols(dist.to.ecdna)$distance

}


## Write annotated kataegis results
write.sigprofiler(ktg, opt$out_file_ktg)
message(opt$out_file_ktg)

# ## Write annotated jabba results
write.table(jba, opt$out_file_jba, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file_jba)
