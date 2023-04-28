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
## Manhattan plot of significant fishhook loci, with nearest CGC genes annotated
libs = c('optparse', 'ggplot2', 'viridis', 'ggrepel')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="Output of annotate-fishhook-results.r"),
  make_option(c("-f", "--fdr"),      type='numeric',   help="FDR cutoff"),
  make_option(c("-o", "--out_file"), type='character', help="Output SVG"))
opt = parse_args(OptionParser(option_list=option_list))



## Read data
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')
dta$chr = factor(dta$chr, levels=paste0('chr',c(1:22,'X')))

## Don't use facets
dta$abs_pos = dta$start / tapply(dta$start, dta$chr, max)[dta$chr]
dta$abs_pos = dta$abs_pos + (as.numeric(dta$chr) - 1)

## log transform Y axis, update chr labels
dta$y = -log10(dta$fdr)
levels(dta$chr) = gsub('chr', '', levels(dta$chr))

## Handle nearest gene formatting
dta$nearest_distance[!dta$significant] = NA
dta$nearest[!dta$significant] = NA

## If multiple bins share the same nearest gene, mark the bin with the best FDR
dta = dta[order(dta$y, decreasing=T), ]
dta$nearest[duplicated(dta$nearest)] = NA
dta$nearest_distance = log10(dta$nearest_distance + 1)



##########
## Plot ##
##########

## Plot formatting 
vlines = seq(0, length(levels(dta$chr))-1)
ticks = seq(0, length(levels(dta$chr))-1) + 0.5

legend.breaks = 0:6
legend.labels = c('1bp', '10bp', '100bp', '1KB', '10KB', '100KB', '>=1MB')


svg(opt$out_file, width=14, height=5)

ggplot(dta, aes(x=abs_pos, y=y, col=nearest_distance)) + 
  geom_point() +
  geom_hline(yintercept=-log10(opt$fdr), alpha=0.15, linetype='solid') +
  geom_vline(xintercept=vlines, alpha=0.4, linetype='dashed') +
  scale_color_viridis(breaks=legend.breaks, labels=legend.labels, limits=c(0,6.01), name='Nearest CGC Gene') +
  xlab('Chromosome') +
  ylab(expression(-log[10](q-value))) +
  scale_x_continuous(limits=c(0,length(levels(dta$chr))), expand=c(0,0), breaks=ticks, labels=levels(dta$chr)) +
  scale_y_continuous(limits=c(0,2), expand=c(0,0)) +
  geom_text_repel(aes(label = nearest),
                      col='black',
                      nudge_x = .1,
                      box.padding = 0.5,
                      nudge_y = 0.1,
                      segment.alpha=1,
                      size=2.5) +
  ggtitle('Nearest CGC genes to loci with significantly recurrent breakpoints as detected by FishHook') + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.spacing=unit(0, 'mm'))

dev.off()
message(opt$out_file)
