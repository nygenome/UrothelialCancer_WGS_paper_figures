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
## Call events on a junction-balanced genome graph
libs = c('optparse', 'gGnome', 'parallel')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


## Get arguments
option_list = list(
  make_option(c("-g", "--jabba_gg"), type='character', help="JaBbA simplified genome graph"),
  make_option(c("-o", "--out_file"), type='character', help="Output RDS"))
opt = parse_args(OptionParser(option_list=option_list))


## Generate some output filenames based on the main output file (a marked genome graph)
summary.file.events = gsub('rds$','summary.rds',opt$out_file)
summary.file.junctions = gsub('rds$','junction.summary.rds',opt$out_file)
summary.file.mh = gsub('rds$','microhomology.junction.summary.rds',opt$out_file)
summary.file.footprints = gsub('rds$','footprints.rds',opt$out_file)
summary.file.footprints.txt = gsub('rds$','footprints.txt',opt$out_file)


## Read jabba result
dta = readRDS(opt$jabba_gg)
dta = dta[seqnames %in% paste0('chr',c(1:22,'X')), ]


## Call events
res = gGnome::events(dta, QRP=TRUE)



## Get event count, per-event junction count, and footprints
sv.classes = c('tra','inv','invdup','del','dup','pyrgo','rigma','chromoplexy','chromothripsis','tic','qrp','bfb','dm','tyfonas','cpxdm')
counts = c()
junctions = c()
footprints = c()

for (s in sv.classes) {

  if (s %in% c('dm','fbi','bfb','tyfonas','cpxdm')) {

    fp = res$meta$amp$footprint[res$meta$amp$type == s]

  } else if (s %in% c('tra','inv','invdup')) {

    fp = res$meta$simple$footprint[res$meta$simple$type == s]

  } else if (s == 'qrp' ) {

    ## QRP is split into qrp-positive, qrp-minus, qrp-mixed
    qrp.mta = res$meta[c('qrpmin','qrpmix','qrppos')]
    if(all(!sapply(qrp.mta, is.null))) {
      fp = do.call(function(...) rbind(..., use.names=F), res$meta[c('qrpmin','qrpmix','qrppos')])$footprint
    } else {
      fp = NULL
    }

  } else {

    fp = res$meta[[s]]$footprint
    
  }

  ## Event count, junction count (by event)
  count = length(fp)
  counts = c(counts, count)
  
  ## Individual events called by simple() aren't annotated directly,
  ## but we can use the event count because TRA are 1 junction, INV/INVDUP are 2
  if (s %in% c('tra')) {

    js = count

  } else if (s %in% c('inv','invdup')) {

    js = count * 2

  } else if (s == 'qrp') {
    
    js = sum(apply(res$edges$dt[, c('qrpmin','qrpmix','qrppos')], 1, function(x) any(!is.na(x))))

  } else {

    js = sum(!is.na(unlist(res$edges$dt[[s]])))
  }
  junctions = c(junctions, js)

  ## Footprints
  fp = data.frame(type=rep(s, tail(counts, 1)), footprint=fp)
  footprints = rbind(footprints, fp)

}

names(junctions) = sv.classes
names(counts) = sv.classes



## Get count of unclassified junctions and total junction count 
sv.classes = c('simple','del','dup','pyrgo','rigma','chromoplexy','chromothripsis','tic','bfb','qrpmix','qrpmin','qrppos','dm','tyfonas','cpxdm')
edge.df = as.data.frame(res$edges[type == 'ALT']$dt, stringsAsFactors=F)


if (nrow(edge.df) == 0) {
  junctions['unclassified'] = 0
} else {
  junctions['unclassified'] = sum(apply(edge.df[,sv.classes], 1, function(y) all(is.na(y))))
}

junctions['njunc'] = nrow(edge.df)


## Write results
saveRDS(res, opt$out_file)
saveRDS(counts, summary.file.events)
saveRDS(junctions, summary.file.junctions)
saveRDS(footprints, summary.file.footprints)
write.table(footprints, summary.file.footprints.txt, col.names=T, row.names=F, sep='\t', quote=F)
