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
## Plot VAF distributions for mutations on ecDNAs

libs = c('optparse', 'GenomicRanges', 'data.table', 'ggplot2', 'ggpubr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))



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
ktg = Reduce(c, lapply(tn$in_files_ktg, read.sigprofiler))
nc = Reduce(c, lapply(tn$in_files_nc, read.sigprofiler))


## Select ecDNA mutations
ktg = ktg[ktg$ecdna]


## Collapse each kataegis event, using mean VAF
## This way each event is only represented once 
mean.vaf = tapply(ktg$vcf_hc_vaf, ktg$kataegis_id, mean)
ktg = ktg[!duplicated(ktg$kataegis_id)]
ktg$vaf = unname(mean.vaf[ktg$kataegis_id])

nc = nc[nc$jabba_event == 'ecDNA']


## Pool chemo, exclude everything else
nc$assigned_signature[nc$assigned_signature %in% paste0('SBS', c(31,35))] = 'Platinum chemotherapy\n(SBS31 + SBS35)'
nc = nc[nc$assigned_signature %in% c('Platinum chemotherapy\n(SBS31 + SBS35)')]


## Combine kataegic/non-clustered mutations
ktg$assigned_signature = 'Kyklonas'
nc = as.data.frame(mcols(nc))[, c('assigned_signature', 'vaf')]
ktg = as.data.frame(mcols(ktg))[, c('assigned_signature', 'vaf')]
res = rbind(nc, ktg)


## Assign colors and enforce ordering
col = c(`Kyklonas`='#26567e', `Platinum chemotherapy\n(SBS31 + SBS35)`='#e41a1c')
res$assigned_signature = factor(res$assigned_signature, levels=names(col))


## Generate list of comparisons to test
compr = combn(names(col), m=2)
compr = lapply(1:ncol(compr), function(i) compr[,i])



##########
## Plot ##
##########

svg(opt$out_file, width=4.25, height=4)

ggplot(res, aes(x=assigned_signature, y=vaf, fill=assigned_signature)) +
  geom_violin() +
  geom_boxplot(width=0.075, fill='white', outlier.shape=NA) +
  scale_fill_manual(values=col) +
  scale_y_continuous(limits=c(0,1.05), breaks=c(0, 0.5, 1)) +
  stat_compare_means(comparisons=compr, method='wilcox') +
  theme_bw() +
  xlab('Mutation Type') +
  ylab('VAF') +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position='none',
        text=element_text(size=17))

dev.off()
message(opt$out_file)



###########################
## Write out source data ##
###########################

ct = table(res$assigned_signature)
res = as.data.frame(do.call(rbind, tapply(res$vaf, res$assigned_signature, summary)))
res$n = as.numeric(ct)

out.file.txt = gsub('\\.svg$', '.txt', opt$out_file)

write.table(res, out.file.txt, row.names=T, col.names=T, quote=F, sep='\t')
message(out.file.txt)
