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
## Run fishhook to check for regions of recurrent structural variation
libs = c('optparse', 'GenomicRanges', 'gUtils', 'fishHook')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),               type='character', help="Merged junction BED file"),
  make_option(c("-c", "--chr_len"),               type='character', help="Chromosome lengths"),
  make_option(c("-e", "--eligible"),              type='character', help="BED file with eligible regions"),
  make_option(c("-v", "--covariate_dir"),         type='character', help="Directory with covariate GRanges"),
  make_option(c("-b", "--bin_size"),              type='numeric',   help="Bin size", default=1E5),
  make_option(c("-s", "--step_size"),             type='numeric',   help="Bin step size", default=5E4), #5E4
  make_option(c("-x", "--one_patient_per_locus"), type='logical',   help="Restrict to one breakpoint per-patient per-locus?", default=TRUE),
  make_option(c("-o", "--out_file"),              type='character', help="Output RDS"))
opt = parse_args(OptionParser(option_list=option_list))


## Read data
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')
dta = makeGRangesFromDataFrame(dta, keep.extra.columns=T)
dta$patient.junction.id = paste0(dta$patient,'.', 1:length(dta))


## Get bins
chr.len = read.table(opt$chr_len, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','end'))
chr.len$start = 1
chr.len = chr.len[!chr.len$chr %in% c('chrM', 'chrY'), ]
chr.len = makeGRangesFromDataFrame(chr.len)
bins = gr.tile(chr.len, opt$bin_size)

if (opt$step_size < opt$bin_size) {

  message('Using step size: ', opt$step_size)

  ## Original bins, plus those shifted right by the step size
  ## Make sure to exclude bins that fall off the edge of the chromosome
  bins = c(bins, gUtils::`%+%`(bins ,opt$step_size))
  bins = bins[bins %O% chr.len == 1]

  ## Have to update IDs to avoid duplicates
  bins$tile.id = 1:length(bins)

}


## Filter for one breakpont per patient per bin
bins$tile.id = as.character(bins$tile.id)
dta = gr.val(dta, bins, 'tile.id', sep=',')

if (opt$one_patient_per_locus) {
  message('Filtering to one patient per locus')
  dta = dta[!duplicated(paste0(dta$patient,'--',dta$tile.id))]
}


## Set eligible regions for calling 
eligible = data.table::fread(opt$eligible, header=F, stringsAsFactors=F, sep='\t', col.names=c('chr','start','end'))
eligible = makeGRangesFromDataFrame(eligible)


## Get covariates
cov.files = Sys.glob(paste0(opt$covariate_dir,'/covariate_*',opt$bin_size,'.rds'))
cov.names = gsub('_[0-9]+\\.rds', '', gsub('covariate_', '', basename(cov.files)))

covariates = c()

for (i in 1:length(cov.files)) {

  cov = readRDS(cov.files[i])

  ## Fishhook doesn't like special characters
  cov.names[i] = gsub('-','_',cov.names[i])
  colnames(mcols(cov)) = gsub('-','_',colnames(mcols(cov)))

  type = class(mcols(cov)[[1]])
  if (type == 'character') {
    type = 'interval'
  }

  ## Build covariate object and add to list 
  cov = Cov(cov, field=colnames(mcols(cov)), name=cov.names[i], type=type)

  print(cov)

  if (length(covariates) == 0) {
    covariates = cov
  } else {
    covariates = c(covariates, cov)
  }

}


## Build fishhook object and run
fish = Fish(hypotheses=bins, 
            events=dta, 
            eligible=eligible,
            covariates=covariates,
            mc.cores = 2)

fish = fish[fish$data$frac.eligible > 0.75, ]
fish$score()
fish$qqp(plotly = FALSE)

saveRDS(fish, opt$out_file)
