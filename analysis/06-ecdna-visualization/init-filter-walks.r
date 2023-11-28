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
## Filter decomposed walks for circular walks across high-copy nodes
libs = c('optparse', 'reshape2', 'GenomicRanges', 'gGnome')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


HIGHCOPY = 8
RO_CUTOFF = 0 ## More stringent filtering during walk export

## Get arguments
option_list = list(
  make_option(c("-s", "--sample"),             type='character', help="Tumor-normal pair"),
  make_option(c("-i", "--in_file_aa"),         type='character', help="Output of init-summarize-aa-results.r"),
  make_option(c("-w", "--in_file_walks"),      type='character', help="RDS with walks from init-decompose-amplicons.r"),
  make_option(c("-o", "--out_file"),           type='character', help="Output RDS"))
opt = parse_args(OptionParser(option_list=option_list))


## Read walks
walk.list = readRDS(opt$in_file_walks)

## Read AA intervals 
aa = makeGRangesFromDataFrame(read.csv(opt$in_file_aa, h=T, stringsAsFactors=F, sep='\t'), keep.extra=T)
aa$sample = gsub('amplicon.*_(?=PM)', '', aa$id, perl=T)
aa = aa[aa$sample %in% opt$sample & aa$class == 'ecDNA']


## We have a set of walks for each amplicon that 
## overlapped an AA ecDNA call
for (i in 1:length(walk.list)) {

  walks = walk.list[[i]]

  ## Filter for circular high-copy walks 
  walks = walks[circular]
  if (length(walks) == 0) {
    walk.list[[i]] = walks
    next
  }

  ## Check to make sure all nodes are amplified 
  walks = walks[walks$eval(all(cn >= HIGHCOPY))]
    if (length(walks) == 0) {
    walk.list[[i]] = walks
    next
  }

  ## For walks that passed initial filters, 
  ## check reciprocal overlap with associated AA events
  ro = c()
  aa.events = c()
  for (j in 1:length(walks)) {

    walk = walks[j]
    walk.aa.events = unique(aa$id[aa %^% walk$footprint])

    if (length(walk.aa.events) > 0) {
      
      bp.overlap = tapply(aa[aa$id %in% walk.aa.events] %o% walk$footprint, aa$id[aa$id %in% walk.aa.events], sum)
      aa.width = tapply(width(aa[aa$id %in% walk.aa.events]), aa$id[aa$id %in% walk.aa.events], sum)

      walk.ro = (2 * bp.overlap) / (sum(width(walk$footprint)) + aa.width)

      ## We might still be overlapping more than one AA call
      ## Take the one with the best reciprocal overlap
      walk.ro = walk.ro[which.max(walk.ro)]
      walk.aa.events = paste(names(walk.ro), collapse=',')

    } else {

      walk.ro = 0
      walk.aa.events = NA

    }

    ro = c(ro, walk.ro)
    aa.events = c(aa.events, walk.aa.events)
  
  }
  
  walks$set(reciprocal_overlap=ro)
  walks$set(aa_events=aa.events)

  ## Filter on reciprocal overlap
  walks = walks[reciprocal_overlap > RO_CUTOFF]

  ## If multiple walks overlap the same AA call, choose the one with the highest reciprocal overlap
  walks = walks[order(walks$dt$aa_events, walks$dt$reciprocal_overlap, decreasing=T)]
  walks = walks[!duplicated(walks$dt$aa_events)]
  walk.list[[i]] = walks

}


## Save result
saveRDS(walk.list, opt$out_file)
