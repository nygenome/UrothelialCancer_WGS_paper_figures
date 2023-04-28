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
## Plot comparison of JCN between pre/post-chemo samples

libs = c('optparse', 'ggplot2', 'reshape2', 'ggpubr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="Summary of ecDNA features, annotated with sample metadata"),
  make_option(c("-m", "--metadata"), type='character', help="Sample metadata"),
  make_option(c("-o", "--out_file"), type='character', help="SVG output"))
opt = parse_args(OptionParser(option_list=option_list))


plt.features = c('Mean JCN', 'Max JCN')


## Read data
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')
dta$num_cgc_genes = sapply(strsplit(dta$cgc_genes, ',', fixed=T), function(x) sum(x != ''))

## Read metadata 
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')

## Filter for ecDNAs less than 3MB in length
dta = dta[dta$footprint_mb < 3, ]

dta = dta[, !colnames(dta) %in% c('type', 'final_type', 'has_kataegis', 'cgc_genes', 'footprint')]
dta.orig = dta

dta = melt(dta, id.var=c('sample', 'patient', 'platinum_chemotherapy'))

dta$platinum_chemotherapy = factor(dta$platinum_chemotherapy, levels=c('Pre-chemo','Post-chemo'))



## Reformat names 
name.reformat = c(`max_cn`='Max interval CN', `num_junctions`='Junction count', `max_jcn`='Max JCN', 
                  `jcn_shannon_entropy`='JCN Shannon entropy', `genomic_mass_mb`='Genomic mass (MB*CN)', 
                  `footprint_mb`='Genomic footprint (MB)', `mean_segment_size_mb`='Mean segment size (MB)', 
                  `sd_segment_size_mb`='Segment size std dev (MB)', `num_cgc_genes`='CGC gene count',
                  `width_weighted_cn`='Width-weighted CN', `mean_jcn`='Mean JCN')


levels(dta$variable) = name.reformat[levels(dta$variable)]
dta = dta[as.character(dta$variable) %in% plt.features, ]

dta$variable = factor(dta$variable, levels=plt.features)


## Summary stats for the text
for (i in plt.features) {

  message(i)
  print(tapply(dta$value[dta$variable == i], dta$platinum_chemotherapy[dta$variable == i], summary))

}



##########
## Plot ##
##########

svg(opt$out_file, width=8, height=5)

col = c(`Pre-chemo`='#BC8AB6',`Post-chemo`='#8E305B')

ggplot(dta, aes(x=platinum_chemotherapy, y=value, fill=platinum_chemotherapy)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.15) +
  scale_fill_manual(values=col) +
  stat_compare_means(method='wilcox.test', 
                    label='p.format',
                    label.x=1.4) +
  facet_wrap(. ~ variable, 
             scale='free', 
             strip.position='left') +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        text=element_text(size=15),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_blank(),
        strip.placement='outside',
        legend.position='none') 

dev.off()
message(opt$out_file)


## Sort and write out source data
out.file.txt = gsub('\\.svg$', '.txt', opt$out_file)
dta = dta[order(dta$variable, dta$platinum_chemotherapy), ]

dta$sample = mta$sample_id[match(dta$sample, mta$tumor)]

write.table(dta, out.file.txt, row.names=F, col.names=T, quote=F, sep='\t')
message(out.file.txt)
