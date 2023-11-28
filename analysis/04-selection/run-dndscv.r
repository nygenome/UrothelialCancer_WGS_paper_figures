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
## Run dNdScv to estimate positive and negative selection
libs = c('optparse', 'data.table', 'dndscv')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


QVAL_CUTOFF = 0.1


## Read RDA file with a single object into 
## a variable, similar to readRDS
read.cov.rda = function(f) {

  x = get(load(f))
  return(x)

}


## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="Input tab-delimited file with the columns sampleID, patient, chr, pos, ref, mut"),
  make_option(c("-m", "--metadata"), type='character', help="Metadata file"),
  make_option(c("-r", "--ref"),      type='character', help="dNdScv reference rda file"),
  make_option(c("-c", "--cov"),      type='character', help="dNdScv covariate rda file"),
  make_option(c("-o", "--out_file"), type='character', help="Output RDS"))
opt = parse_args(OptionParser(option_list=option_list))



## Read covariates
cov = read.cov.rda(opt$cov)

## Read metadata
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')

## Read input file, add metadata and key for deduplication
dta = data.table::fread(opt$in_file, h=T, sep='\t', data.table=F, stringsAsFactors=F)
dta$chr = gsub('^chr', '', dta$chr)
dta$key = apply(dta[,c('patient','chr','pos','ref','mut')], 1, paste, collapse='|')

## Init result object
res = list()

## Run dNdScv with all samples
dta = dta[!duplicated(dta$key), c("sampleID", "chr", "pos", "ref", "mut")]
res[['all']] = dndscv::dndscv(dta, cv=cov, refdb=opt$ref)
  


###################
## Write results ##
###################

## Save result object
saveRDS(res, opt$out_file)
message(opt$out_file)

## Export global rates of selection
global.dnds = do.call(rbind, lapply(names(res), function(x) {res[[x]]$globaldnds$cohort=x; res[[x]]$globaldnds}))
out.file.global = gsub('\\.rds$','.globaldnds.csv', opt$out_file)
write.csv(global.dnds, out.file.global, row.names=T, quote=F)

## Export gene-level seleciton metrics
sig.genes = do.call(rbind, lapply(names(res), function(x) {res[[x]]$sel_cv$cohort=x; res[[x]]$sel_cv[res[[x]]$sel_cv$qglobal_cv < QVAL_CUTOFF, ]}))
out.file.sig = gsub('\\.rds$','.siggenes.csv', opt$out_file)
write.csv(sig.genes, out.file.sig, row.names=T, quote=F)
