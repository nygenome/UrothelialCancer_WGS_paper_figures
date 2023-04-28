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
## Plot empirical cumulative distributions of kataegis distance to the nearest junction

libs = c('optparse', 'data.table', 'ggplot2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Read SigProfilerClusters annotated output
read.sigprofiler = function(f) {

  if (!file.exists(f)) {
    return(NULL)
  }

  x = fread(f, header=T, stringsAsFactors=F, sep='\t', data.table=F)

  return(x)

}



## Get arguments
option_list = list(
  make_option(c("-I", "--in_dir_ktg"),  type='character', help="Directory holding SigProfilerClusters kataegis VCFs"),
  make_option(c("-n", "--in_dir_nc"),   type='character', help="Directory holding SigProfilerClusters VCFs (unclustered)"),
  make_option(c("-t", "--tn_file"),     type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-o", "--out_file"),    type='character', help="SVG output"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor-normal pairs and get input files
tn = read.csv(opt$tn, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)

tn$in_files_ktg = paste0(opt$in_dir_ktg,'/',tn$pair_name,'-mutationtimer-snv.annotated.jabba.txt')
tn$in_files_nc = paste0(opt$in_dir_nc,'/',tn$pair_name,'-mutationtimer-snv.annotated.jabba.txt')


## Read input 
ktg = Reduce(rbind, lapply(tn$in_files_ktg, read.sigprofiler))

ktg$jabba_event[ktg$ecdna] = 'ecDNA'

col.sel = c('assigned_signature','jabba_event', 'distance_to_nearest_junction', 'kataegis_id')
ktg$assigned_signature = 'kataegis'
ktg$jabba_event = gsub('\\..*', '', ktg$jabba_event)


## Remove cases where kataegis falls on a chromosome with no breakpoints
ktg = ktg[!is.na(ktg$distance_to_nearest_junction), ]


## Collapse each kataegis event, using mean distance to nearest breakpoint
## This way each event is only represented once 
mean.dist = tapply(ktg$distance_to_nearest_junction, ktg$kataegis_id, median)
ktg = ktg[!duplicated(ktg$kataegis_id), ]
ktg$distance_to_nearest_junction = mean.dist[ktg$kataegis_id]


## Simplify event types
ktg$jabba_event[ktg$jabba_event == 'none'] = 'No SV'
ktg$jabba_event[!ktg$jabba_event %in% c('ecDNA', 'No SV')] = 'JaBbA SV (Non-ecDNA)'
ktg$jabba_event = factor(ktg$jabba_event, levels=rev(c('No SV', 'JaBbA SV (Non-ecDNA)', 'ecDNA')))



###############
## K-S tests ##
###############

comparisons = combn(levels(ktg$jabba_event), m=2)

plt.test.coords = data.frame(x0=rep(NA, ncol(comparisons)),
                             y0=rep(NA, ncol(comparisons)),
                             y1=rep(NA, ncol(comparisons)),
                             p=rep(NA, ncol(comparisons)))


for (i in 1:ncol(comparisons)) {

  dta.compr.1 = ktg$distance_to_nearest_junction[ktg$jabba_event == comparisons[1, i]]
  dta.compr.2 = ktg$distance_to_nearest_junction[ktg$jabba_event == comparisons[2, i]]

  message(comparisons[1, i], ' (n=',length(dta.compr.1),') vs. ', 
          comparisons[2, i], ' (n=',length(dta.compr.2),')\n')

  ## Run test
  ks = unlist(ks.test(x=dta.compr.1, y=dta.compr.2, alternative='two.sided'))

  print(wilcox.test(x=dta.compr.1, y=dta.compr.2, alternative='two.sided'))

  ## Get data for plotting (adapted from https://rpubs.com/mharris/KSplot)
  cdf1 = ecdf(dta.compr.1) 
  cdf2 = ecdf(dta.compr.2)

  ## Find min and max statistics to draw line between points of greatest distance
  minMax = seq(min(dta.compr.1, dta.compr.2), max(dta.compr.1, dta.compr.2), length.out=length(dta.compr.1)) 
  x0 = min(minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )])

  plt.test.coords$x0[i] = x0
  plt.test.coords$y0[i] = cdf1(plt.test.coords$x0[i]) 
  plt.test.coords$y1[i] = cdf2(plt.test.coords$x0[i]) 
  plt.test.coords$p[i] = as.numeric(ks['p.value'])
  plt.test.coords$py[i] = (plt.test.coords$y0[i] + plt.test.coords$y1[i]) / 2

  ## Format p-value
  plt.test.coords$p.original[i] = plt.test.coords$p[i]
  if (plt.test.coords$p[i] <= 0.001) {
    plt.test.coords$p[i] = '***'
  } else if(plt.test.coords$p[i] <= 0.01) {
    plt.test.coords$p[i] = '**'
  } else if(plt.test.coords$p[i] <= 0.05) {
    plt.test.coords$p[i] = '*'
  } else {
    plt.test.coords$p[i] = 'ns'
  }
  
}



##########
## Plot ##
##########

COLORS = c(`No SV`='#bebebe', `JaBbA SV (Non-ecDNA)`='#000000', `ecDNA`='#26567e')

svg(opt$out_file, width=6, height=4)


ggplot() +
  stat_ecdf(data=ktg, aes(x=distance_to_nearest_junction, color=jabba_event), size=1) +
  scale_fill_manual(values=COLORS) +
  scale_color_manual(values=COLORS) +
  scale_x_continuous(trans='log10', 
                     limits=c(min(ktg$distance_to_nearest_junction),1E8),
                     breaks=c(10, 1E2, 1E3, 1E4, 1E5, 1E6, 1E7, 1E8), 
                     labels=c('10bp', '100bp', '1KB', '10KB', '100KB', '1MB', '10MB', '100MB')) +
  
  xlab('Distance to nearest junction') +
  ylab('Proportion of kataegis events') +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        text=element_text(size=17),
        legend.position=c(0.275, 0.875),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title=element_blank())


dev.off()
message(opt$out_file)



###########################
## Write out source data ##
###########################

## Just pull quantiles 0, 0.01, 0.02, ..., 1
res = do.call(cbind, tapply(ktg$distance_to_nearest_junction, ktg$jabba_event, quantile, probs=seq(0,1,by=0.01), na.rm=T))

out.file.txt = gsub('\\.svg$', '.txt', opt$out_file)
write.table(res, out.file.txt, row.names=T, col.names=T, quote=F, sep='\t')
message(out.file.txt)
