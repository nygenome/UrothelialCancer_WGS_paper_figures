#!/nfs/sw/R/R-4.0.0/bin/Rscript
#!/bin/bash
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
## Plot alignment lengths of split reads supporting an ecDNA
libs = c('optparse', 'ggplot2', 'patchwork')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="TSV with alignment lengths of junction-supporting reads"),
  make_option(c("-o", "--out_file"), type='character', help="Figure output file (SVG)"))
opt = parse_args(OptionParser(option_list=option_list))

MIN_SUPPORT=1E3
POINT_COL='#7F0E7F'

dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')



##########
## Plot ##
##########

svg(opt$out_file, width=6.5, height=6.5)

for (i in unique(dta$event_id)) {

  dta.i = dta[dta$event_id == i, ]

  message('Total number of reads: ', nrow(dta.i))
  message('Total number of reads with both split alignments > 1KB: ', sum(dta.i$span_left > 1E3 & dta.i$span_right > 1E3))
  message('Total number of reads with both split alignments > 5KB: ', sum(dta.i$span_left > 5E3 & dta.i$span_right > 5E3))
  message('Total number of reads with both split alignments > 10KB: ', sum(dta.i$span_left > 1E4 & dta.i$span_right > 1E4))

  ## Main scatter plot 
  axis.limits = c(0, 1.05 * max(c(dta.i$span_left, dta.i$span_right)))
  scatter.breaks = c(1E3, 5E3, 1E4, 1.5E4, 2E4, 2.5E4)
  scatter.labels = c('1KB', '5KB', '10KB', '15KB', '20KB', '25KB')
  hist.binwidth=1E3


  plt.main = ggplot(dta.i, aes(x=span_right, y=span_left, col=both_above_cutoff)) +
    geom_point(alpha=0.75, col=POINT_COL) +
    geom_hline(yintercept=MIN_SUPPORT, linetype='dashed', color='#000000') +
    geom_vline(xintercept=MIN_SUPPORT, linetype='dashed', color='#000000') +
    scale_x_continuous(limits=axis.limits, name='ecDNA right flank', breaks=scatter.breaks, labels=scatter.labels) +
    scale_y_continuous(limits=axis.limits, name='ecDNA left flank', breaks=scatter.breaks, labels=scatter.labels) +
    theme_bw() + 
    theme(text=element_text(size=15),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          legend.position='none')

  x = ggplot(dta.i, aes(x=span_right)) + 
        geom_histogram(binwidth=hist.binwidth, center=hist.binwidth/2, col='#000000', lwd=0.2) +
        geom_vline(xintercept=MIN_SUPPORT, linetype='dashed', color='#000000') +
        scale_x_continuous(limits=axis.limits) +
        scale_y_continuous(expand=c(0,0), name='ONT\nReads') +
        theme_bw() +
        theme(text=element_text(size=15),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          legend.position='none') 

  y = ggplot(dta.i, aes(x=span_left)) + 
        geom_histogram(binwidth=hist.binwidth, center=hist.binwidth/2, col='#000000', lwd=0.2) +
        geom_vline(xintercept=MIN_SUPPORT, linetype='dashed', color='#000000') +
        scale_x_continuous(limits=axis.limits) +
        scale_y_continuous(expand=c(0,0), name='ONT\nReads') +
        theme_bw() +
        theme(text=element_text(size=15),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          legend.position='none') +
        coord_flip()
  
  
  plot(
    wrap_plots(
      x, 
      plot_spacer(), 
      plt.main,
      y,
      nrow=2,
      byrow=T,
      widths=c(1, 0.2), 
      heights=c(0.2, 1)
    ))

}

dev.off()
message(opt$out_file)
    