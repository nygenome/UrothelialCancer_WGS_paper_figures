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
## Decompose amplicons overlapping an AA event into constituent walks
libs = c('optparse', 'reshape2', 'GenomicRanges', 'gGnome')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



JABBA_AMPS = c('bfb','dm','tyfonas')



## Get arguments
option_list = list(
  make_option(c("-s", "--sample"),             type='character', help="Tumor-normal pair"),
  make_option(c("-i", "--in_file_comparison"), type='character', help="Output of init-compare-aa-to-jabba.r"),
  make_option(c("-I", "--in_file_aa"),         type='character', help="Output of init-summarize-aa-results.r"),
  make_option(c("-j", "--in_file_jabba"),      type='character', help="Jabba events RDS"),
  make_option(c("-o", "--out_file"),           type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



## Read input file 
dta = read.csv(opt$in_file_comparison, h=T, stringsAsFactors=F, sep='\t')

## Read jabba graph
gg = readRDS(opt$in_file_jabba)

## Read AA intervals 
aa = makeGRangesFromDataFrame(read.csv(opt$in_file_aa, h=T, stringsAsFactors=F, sep='\t'), keep.extra=T)
aa$sample = gsub('amplicon.*_(?=PM)', '', aa$id, perl=T)
aa = aa[aa$sample %in% opt$sample & aa$class == 'ecDNA']

## Subset to ecDNA calls overlapping with jabba amps, sort on reciprocal overlap
dta = dta[!is.na(dta$overlapping_jabba_amp) & dta$class == 'ecDNA', ]
dta = dta[order(dta$overlapping_jabba_amp_ro, decreasing=T), ]

## Add reciprocal overlap to AA intervals
aa$ro = round(dta$overlapping_jabba_amp_ro[match(aa$id, dta$aa_id)], 3)
aa$plt_id = paste(aa$plt_id,'\nRO=',aa$ro)

## Get list of jabba events that we have to plot 
jba.events= dta[,c('sample', 'overlapping_jabba_amp')]
jba.events = jba.events[!duplicated(apply(jba.events, 1, paste, collapse=',')), ]
jba.events = jba.events$overlapping_jabba_amp[jba.events$sample == opt$sample]




if (length(jba.events) > 0) {

  res = list()

  ## For each jabba event that we need to plot 
  for (i in 1:length(jba.events)) {

    ## Get amplification IDs
    amp.type = jba.events[i]
    amp.df = as.data.frame(gg$edges[!is.na(get(amp.type))]$dt)
    amp.ids = unique(amp.df[, amp.type])

    ## For each amp event of this type     
    for (id in amp.ids) {
      
      ## Get associated subgraph
      amp.subg.gr = gg$nodes[!is.na(get(amp.type)) & get(amp.type) == id]$gr
      amp.subg = gg$copy$subgraph(amp.subg.gr, d=1E5)

      ## Make sure it overlaps any % with ecDNA before proceeding 
      if (!any(amp.subg$footprint %^% aa)) {
        next
      }

      ## Peel off CN, preferring flows that maximize JCN
      res[[length(res) + 1]] = peel(amp.subg, field='cn', verbose=T)

    }
    
  }

  saveRDS(res, opt$out_file)
  message(opt$out_file)

}
