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
## Plot proportion of clonal/subclonal signature associated mutations

libs = c('optparse', 'reshape2', 'ggplot2','ggpubr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


SIG_COLS = c('#4daf4a','#377eb8','#4daf4a','#ff7f00','#984ea3','#ffff33','#a65628','#f781bf','#999999','#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3','#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#4e79a7','#f28e2c','#e15759','#76b7b2','#664040','#59a14f','#edc949','#af7aa1','#ff9da7','#b37654','#bab0ab','#1f77b4','#ff7f0e','#2ca02c','#a81d1d','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#664f40','#17becf','#666540','#446640','#40664e','#406660','#405966','#404a66','#464066','#574066','#664065','#664050','#664049','#732222','#737322','#407322','#22733d','#227368','#734822','#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666','#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')
SIGS = c("SBS1","SBS2","SBS3","SBS4","SBS5","SBS6","SBS7a","SBS7b","SBS7c","SBS7d","SBS8","SBS9","SBS10a","SBS10b","SBS10c","SBS10d","SBS11","SBS12","SBS13","SBS14","SBS15","SBS16","SBS17a","SBS17b","SBS18","SBS19","SBS20","SBS21","SBS22","SBS23","SBS24","SBS25","SBS26","SBS27","SBS28","SBS29","SBS30","SBS31","SBS32","SBS33","SBS34","SBS35","SBS36","SBS37","SBS38","SBS39","SBS40","SBS41","SBS42","SBS43","SBS44","SBS45","SBS46","SBS47","SBS48","SBS49","SBS50","SBS51","SBS52","SBS53","SBS54","SBS55","SBS56","SBS57","SBS58","SBS59","SBS60","SBS84","SBS85","SBS86","SBS87","SBS88","SBS89","SBS90","SBS91","SBS92","SBS93","SBS94")

SIG_COLS = SIG_COLS[1:length(SIGS)]
names(SIG_COLS) = SIGS
SIG_COLS = c(SIG_COLS, `Aging\n(SBS1)`='#4daf4a', `Chemotherapy\n(SBS31+SBS35)`='#e41a1c', `APOBEC\n(SBS2+SBS13)`='#377eb8', SBS92='#984ea3')



sigs2foldchange = function(x, timing, min.counts=500) {
  
  cols = c('subclonal', 'clonal [late]', 'clonal [NA]', 'clonal [early]', 'clonal')
  
  ## Reformat
  rownames(x) = x$signature
  x = x[, !colnames(x) %in% c('signature')]
  
  ## Pool chemo/apobec
  x['Chemo', ] = x['SBS31', ] + x['SBS35', ]
  x['APOBEC', ] = x['SBS2', ] + x['SBS13', ]
  x['Tobacco', ] = x['SBS4', ] + x['SBS92', ]
  
  ## Check for SBS with low counts  
  x$low_counts = rowSums(x) < min.counts
  
  ## Pool clonal 
  x$clonal = x$`clonal [early]` + x$`clonal [late]` + x$`clonal [NA]`
  
  ## Need to normalize by the total mutations in each timing category
  timing = structure(timing$count, names=timing$timing)
  timing['clonal'] = sum(timing[!names(timing) %in% 'subclonal'])
  timing = timing[cols]
  
  for (i in 1:nrow(x)) {
    x[i, cols] = x[i, cols] / timing    
  }
  
  ## Compute fold changes
  x$clonal_subclonal = x$clonal / x$subclonal
  x$early_late =  x$`clonal [early]` / x$`clonal [late]`
    
  ## Reformat   
  x$sbs = rownames(x)
  x = x[,c('sbs','clonal_subclonal', 'early_late', 'low_counts')]
  x = melt(x, id.vars=c('sbs','low_counts'), variable.name='comparison', value.name='fold_change')
  
  return(x)
  
}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir"),          type='character', help="Input directory with MutationTimer results"),
  make_option(c("-t", "--tn_file"),         type='character', help="Tumor normal pairing file"),
  make_option(c("-m", "--metadata"),        type='character', help="Sample metadata"),
  make_option(c("-u", "--timing_summary"),  type='character', help="MutationTimer summary with number of mutations per timing category"),
  make_option(c("-s", "--suffix"),          type='character', help="Which signatures are we showing?", default='-signature-clonality-counts.csv'),
  make_option(c("-o", "--out_file"),        type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))

sig.keep = c('APOBEC', 'SBS1', 'Chemo')

## Read metadata
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')

## Read mutationtimer summary
timing = read.csv(opt$timing_summary, h=T, stringsAsFactors=F, sep='\t')

## Read TN file, select post-chemo samples
tn = read.csv(opt$tn_file, h=F, sep='\t', stringsAsFactors=F, col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor, '--', tn$normal)
tn$in_file = paste0(opt$in_dir, '/', tn$pair_name, opt$suffix)

tn = tn[tn$tumor %in% mta$tumor, ]

if(!all(file.exists(tn$in_file))) {
  stop('Missing files!')
}


## Read files, tabulate mutations per signature/clonality
dta = lapply(tn$in_file, read.csv, h=T, stringsAsFactors=F, check.names=F)

res = lapply(1:length(dta), function(i) sigs2foldchange(
                              x=dta[[i]], 
                              timing=timing[timing$sample == tn$tumor[i], ]))
res = do.call(rbind, res)

res$sample = rep(tn$tumor, each=(nrow(dta[[1]])+3)*2)



## Tabulate mutations per signature/sample
for(i in 1:length(dta)) {
  
  dta[[i]] = data.frame(signature=dta[[i]]$signature, 
                        sample=tn$tumor[i],
                        count=as.character(rowSums(dta[[i]][,-1])))
  
  
  dta[[i]] = rbind(dta[[i]], c('Chemo',   tn$tumor[i], sum(as.numeric(dta[[i]]$count[dta[[i]]$signature %in% c('SBS31', 'SBS35')]))))
  dta[[i]] = rbind(dta[[i]], c('APOBEC',  tn$tumor[i], sum(as.numeric(dta[[i]]$count[dta[[i]]$signature %in% c('SBS2', 'SBS13')]))))
  dta[[i]] = rbind(dta[[i]], c('Tobacco', tn$tumor[i], sum(as.numeric(dta[[i]]$count[dta[[i]]$signature %in% c('SBS4', 'SBS92')]))))
    
}

dta = do.call(rbind, dta)
dta$count = as.numeric(dta$count)


## Add some metadata
res[, c('Chemotherapy', 'Tobacco')] = mta[match(res$sample, mta$tumor), c('platinum_chemotherapy', 'smoking_status')]
res$Tobacco = factor(res$Tobacco, levels=c('Never', 'Former', 'Current', 'N/A'))

## Select signatures
res = res[!is.nan(res$fold_change) & !is.infinite(res$fold_change), ]
res$comparison = factor(gsub('_', '/', res$comparison), levels=c('early/late','clonal/subclonal'))

res = merge(x=res, y=dta, by.x=c('sbs','sample'), by.y=c('signature','sample'), all.x=T)

res$sbs[res$sbs == 'Chemo'] = 'Chemotherapy\n(SBS31+SBS35)'
res$sbs[res$sbs == 'APOBEC'] = 'APOBEC\n(SBS2+SBS13)'
res$sbs[res$sbs == 'SBS1'] = 'Aging\n(SBS1)'



################################
## Chemo-treated samples only ## 
################################

## Select chemo, manually-defined ordering
sig.order = c('APOBEC\n(SBS2+SBS13)', 'Aging\n(SBS1)', 'Chemotherapy\n(SBS31+SBS35)')
res.chemo = res[res$Chemotherapy == 'Post-chemo' & res$sbs %in% sig.order, ]
res.chemo$sbs = factor(res.chemo$sbs, levels=sig.order)

compr = combn(rev(unique(as.character(res.chemo$sbs))), m=2)
compr = lapply(1:ncol(compr), function(i) compr[,i])
compr


## Early/late
ano = data.frame(x=levels(res.chemo$sbs),
                 lbl=tapply(res.chemo$count[res.chemo$comparison=='early/late'], res.chemo$sbs[res.chemo$comparison=='early/late'], function(x) round(mean(x, na.rm=T))))
ano$lbl[1] = paste0('mean count:\n', ano$lbl[1])
ano$lbl[2:nrow(ano)] = paste0('\n', ano$lbl[2:nrow(ano)])

chemo.early.late = ggplot(res.chemo[res.chemo$comparison=='early/late', ], aes(x=sbs, y=fold_change, fill=sbs)) +
  geom_boxplot(outlier.shape=NA) +
  geom_hline(yintercept=1) +
  geom_jitter(alpha=0.75, height=0, width=0.1) +
  # geom_text(data=ano, aes(x=factor(x, levels=x), y=100, label=lbl), inherit.aes=F) +
  stat_compare_means(comparisons=compr, method='wilcox', step.increase=1.5, tip.length=0) +
  scale_fill_manual(values=SIG_COLS) +
  scale_y_continuous(limits=c(0.01,300), breaks=c(0.01,0.1,1,10,100)) +
  coord_trans(y='log10') +
  xlab('COSMIC SBS') + 
  ylab('Fold Change [Early Clonal / Late Clonal]') +
  theme_bw() + 
  theme(legend.position='none',
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())


## Clonal/subclonal
chemo.clonal.subclonal = ggplot(res.chemo[res.chemo$comparison=='clonal/subclonal', ], aes(x=sbs, y=fold_change, fill=sbs)) +
  geom_boxplot(outlier.shape=NA) +
  geom_hline(yintercept=1) +
  geom_jitter(alpha=0.75, height=0, width=0.1) +
  # geom_text(data=ano, aes(x=factor(x, levels=x), y=100, label=lbl), inherit.aes=F) +
  stat_compare_means(comparisons=compr, method='wilcox', step.increase=1.5, tip.length=0) +
  scale_fill_manual(values=SIG_COLS) +
  scale_y_continuous(limits=c(0.01,300), breaks=c(0.01,0.1,1,10,100)) +
  coord_trans(y='log10') +
  xlab('COSMIC SBS') +
  ylab('Fold change [Clonal / Subclonal]') +
  theme_bw() + 
  theme(legend.position='none',
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())



############################
## All samples/signatures ## 
############################

## Select non-chemo
res = res[res$Chemotherapy == 'Pre-chemo', ]

## Order by median
m = tapply(res$fold_change, list(res$comparison, res$sbs), median, na.rm=T)

## Early/late
res$sbs = factor(res$sbs, levels=colnames(m)[order(m['early/late',])])
pre.early.late = ggplot(res[res$comparison=='early/late', ], aes(x=sbs, y=fold_change, fill=sbs)) +
  geom_boxplot(outlier.shape=NA) +
  geom_hline(yintercept=1) +
  geom_jitter(alpha=0.75, height=0, width=0.1) +
  scale_fill_manual(values=SIG_COLS) +
  scale_y_continuous(trans='log10', limits=c(0.01,250), breaks=c(0.1,1,10,100)) +
  xlab('COSMIC SBS') +
  ylab('Fold Change [Early Clonal / Late Clonal]') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position='none')

## Clonal/subclonal
res$sbs = factor(res$sbs, levels=colnames(m)[order(m['clonal/subclonal',])])
pre.clonal.subclonal = ggplot(res[res$comparison=='clonal/subclonal', ], aes(x=sbs, y=fold_change, fill=sbs)) +
  geom_boxplot(outlier.shape=NA) +
  geom_hline(yintercept=1) +
  geom_jitter(alpha=0.75, height=0, width=0.1) +
  scale_fill_manual(values=SIG_COLS) +
  scale_y_continuous(trans='log10', limits=c(0.01,250), breaks=c(0.1,1,10,100)) +
  xlab('COSMIC SBS') +
  ylab('Fold change [Clonal / Subclonal]') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position='none')


## Plot PDF of all plots
pdf(opt$out_file, width=8, height=8)

chemo.early.late
chemo.clonal.subclonal
pre.early.late
pre.clonal.subclonal

dev.off()
message(opt$out_file)


## Plot individual SVGs for inclusion in manuscript figures 
svg(gsub('\\.pdf$', '.chemo.early_late.svg', opt$out_file), width=5, height=6)
chemo.early.late
dev.off()

svg(gsub('\\.pdf$', '.chemo.clonal_subclonal.svg', opt$out_file), width=5, height=6)
chemo.clonal.subclonal
dev.off()

svg(gsub('\\.pdf$', '.prechemo.early_late.svg', opt$out_file), width=10, height=6)
pre.early.late
dev.off()

svg(gsub('\\.pdf$', '.prechemo.clonal_subclonal.svg', opt$out_file), width=10, height=6)
pre.clonal.subclonal
dev.off()


res.chemo = res.chemo[, c('sample', 'comparison', 'Chemotherapy', 'sbs', 'fold_change')]
res.chemo$sample = mta$sample_id[match(res.chemo$sample, mta$tumor)]

res.chemo$sbs = gsub('\n', ' ', res.chemo$sbs, fixed=T)
write.table(res.chemo, gsub('\\.pdf', '.txt', opt$out_file), row.names=F, col.names=T, quote=F, sep='\t')
message(gsub('\\.pdf', '.txt', opt$out_file))


## Summary stats for source data
med = as.data.frame.matrix(tapply(res.chemo$fold_change, list(res.chemo$comparison, res.chemo$sbs), function(x) round(median(x), 3)))
mn = as.data.frame.matrix(tapply(res.chemo$fold_change, list(res.chemo$comparison, res.chemo$sbs), function(x) round(min(x), 3)))
mx = as.data.frame.matrix(tapply(res.chemo$fold_change, list(res.chemo$comparison, res.chemo$sbs), function(x) round(max(x), 3)))

for (i in 1:nrow(med)) {
  for (j in 1:ncol(med)) {

    med[i, j] = paste0(med[i, j], ' (',mn[i,j],'-',mx[i,j], ')')

  }
}

write.table(t(med), sep='\t', row.names=T, quote=F, col.names=T)
