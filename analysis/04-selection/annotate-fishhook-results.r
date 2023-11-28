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
## Annotate FishHook output 
libs = c('optparse', 'GenomicRanges', 'gUtils', 'fishHook')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),         type='character', help="RDS with fishook result object"),
  make_option(c("-f", "--fdr"),             type='numeric',   help="FDR cutoff", default=0.1),
  make_option(c("-d", "--distance_cutoff"), type='numeric',   help="Distance to nearest gene cutoff", default=1E6),
  make_option(c("-g", "--genes"),           type='character', help="BED file with eligible regions"),
  make_option(c("-c", "--cytoband"),        type='character', help="UCSC cytoband file"),
  make_option(c("-o", "--out_file"),        type='character', help="Output TSV"))
opt = parse_args(OptionParser(option_list=option_list))


## Read data, filter by FDR
dta = readRDS(opt$in_file)
dta = dta$res
dta$significant = dta$fdr < opt$fdr

## Read CGC genes
genes = read.csv(opt$genes, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','start','end', 'gene'))
genes = makeGRangesFromDataFrame(genes, keep.extra.columns=T)
genes$gene = gsub('\\|.*', '', genes$gene)

## Read cytoband
cyto = read.csv(opt$cytoband, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','start','end','cytoband','stain'))
cyto = makeGRangesFromDataFrame(cyto, keep.extra.columns=T)

## Annotate intersecting and nearest genes
dta = dta %$% genes
dta = gr.val(dta, cyto, 'cytoband')


## Nearest gene is subject to distance cutoff
dta$nearest = ''
dta$nearest_distance = NA

nearest.genes = distanceToNearest(x=dta, subject=genes)
mcols(nearest.genes)$distance[mcols(nearest.genes)$distance > opt$distance_cutoff] = opt$distance_cutoff
dta$nearest[queryHits(nearest.genes)] = genes$gene[subjectHits(nearest.genes)]
dta$nearest_distance[queryHits(nearest.genes)] = mcols(nearest.genes)$distance


## Convert to data frame
dta = as.data.frame(dta)
colnames(dta)[1] = 'chr'


## Write result
write.table(dta, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)
