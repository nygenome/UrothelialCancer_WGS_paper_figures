#!/nfs/sw/R/R-4.0.0/bin/Rscript
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: Timothy R. Chu

################################################################# /COPYRIGHT ###
################################################################################
## Plot normal urothelium -> tumor barplots

libs = c('optparse')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="Tab delimited file with mutation counts (urothelium_counts.txt)"),
  make_option(c("-o", "--out_file"), type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



## Read data
counts = read.table(opt$in_file, header=T)



##########
## Plot ##
##########
svg(opt$out_file)

barplot(counts, 
          main="Fig 1D 1E",
          xlab="Number of SNV", 
          col=c("darkblue","blue"),
          legend=rownames(counts))

dev.off()
message(opt$out_file)