#!/nfs/sw/R/R-4.0.0/bin/Rscript
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: Rahul R. Singh

################################################################# /COPYRIGHT ###
################################################################################
## Plot barplot showing the prevalence of ecDNA+ amplification of CCND1 

libs = c('optparse', 'ggplot2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="Input file with amplicon counts for each cancer type"),
  make_option(c("-o", "--out_file"), type='character', help="Figure output file (SVG)"))
opt = parse_args(OptionParser(option_list=option_list))


dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')


##########
## Plot ##
##########

svg(opt$out_file, width=6, height=4)

ggplot(dta, aes(x=typenum, y=number, fill=Eventnum)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("#DADAEB", "#9E9AC8", "#6A51A3"))

dev.off()
message(opt$out_file)
