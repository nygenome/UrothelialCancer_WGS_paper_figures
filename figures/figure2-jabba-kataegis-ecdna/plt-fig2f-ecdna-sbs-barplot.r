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
## Plot mutational signatures for ecDNA mutations

libs = c('optparse', 'reshape2', 'ggplot2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="deconstructSigs COSMIC SBS results for ecDNA mutations binned by VAF"),
  make_option(c("-m", "--metadata"), type='character', help="Sample metadata"),
  make_option(c("-o", "--out_file"), type='character', help="SVG output"))
opt = parse_args(OptionParser(option_list=option_list))



## Read metadata and align with tn pairs
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')

sbs.highlight = c(`SBS2`='#377eb8', `SBS13`='#2166ac',
                  `SBS31`='#e41a1c', `SBS35`='#b2182b')


## Read data and reshape 
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')
dta$Sample = factor(dta$Sample, levels=dta$Sample[order(dta$SBS31 + dta$SBS35, dta$SBS2 + dta$SBS13, decreasing=T)])
dta$VAF = colsplit(dta$Sample,'\\.', c('Sample', 'VAF'))$VAF
dta = melt(dta, id.vars=c('Sample', 'VAF'), variable.name='SBS', value.name='contribution')


## Only highlight specific signatures
levels(dta$SBS)[!levels(dta$SBS) %in% names(sbs.highlight)] = 'Other'


## Update VAF labels 
vaf.labels = c(`vaf_0_333`='VAF <= 0.333', `vaf_333_667`='0.333 < VAF <= 0.667', `vaf_667_1`='VAF > 0.667')
dta$VAF = factor(vaf.labels[dta$VAF], levels=vaf.labels)



##########
## Plot ##
##########

## Assign color for non-highlighted signatures
sbs.highlight['Other'] = 'grey90'

svg(opt$out_file, width=10, height=4)

ggplot(dta, aes(x=Sample, y=contribution, fill=SBS)) +
  geom_bar(stat='identity') +
  facet_grid(. ~ VAF, space='free', scale='free', switch='both') +
  scale_fill_manual(values=sbs.highlight) +
  scale_y_continuous(expand=c(0,0), name='Signature contribution') +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.spacing=unit(3, 'mm'),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        text=element_text(size=15),
        legend.position='right',
        legend.title=element_blank(),
        strip.background = element_rect(color="#FFFFFF", fill="#FFFFFF"))

dev.off()
message(opt$out_file)



## Reformat table for source data
dta = dcast(formula=Sample*VAF ~ SBS, data=dta, fun.aggregate=sum)
dta$Sample = gsub('--.*', '', dta$Sample)
dta$Sample = mta$sample_id[match(dta$Sample, mta$tumor)]

## Write source data
out.file.txt = gsub('\\.svg', '.txt', opt$out_file)
write.table(dta, out.file.txt, row.names=F, col.names=T, quote=F, sep='\t')
message(out.file.txt)