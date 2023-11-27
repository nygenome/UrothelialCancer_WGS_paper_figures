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
## Generate dryclean PON following the README tutorial: https://github.com/mskilab/dryclean
libs = c('optparse', 'parallel', 'magrittr', 'dryclean', 'data.table', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


PARFILE='/gpfs/commons/resources/GRCh38_full_analysis_set_plus_decoy_hla/internal/dryclean/par.rds'
warning('This script is only configured to work on GRCh38 samples!')


## Collect arguments
option_list = list(
  make_option(c("-n", "--normals"),        type='character',                   help="2 column (sample,path), headerless, tab-delimited list of normal samples and paths to their fragCounter output"),
  make_option(c("-o", "--outdir"),         type='character',                   help="Directory to save results"),
  make_option(c("-s", "--cluster_subset"), action='store_true', default=FALSE, help="Use clustering-based approach to select the most informative samples"),
  make_option(c("-c", "--cores"),          type='numeric',                     help="Number of cores to use"))
opt = parse_args(OptionParser(option_list=option_list))


## Normal table needs to be an RDS data.table
dta.normal = read.csv(opt$normals, h=F, stringsAsFactors=F, sep='\t')
colnames(dta.normal) =  c('sample','normal_cov')
dta.normal = as.data.table(dta.normal)

normal.file = paste0(dirname(opt$normals),'/normals.rds')
saveRDS(dta.normal, normal.file)


## Run PON generation 
if (opt$cluster_subset) {
  dryclean::prepare_detergent(normal.table.path=normal.file,
                              path.to.save=opt$outdir,
                              num.cores=opt$cores,
                              use.all=T,
                              save.pon=T,
                              build='hg38',
                              PAR.file=PARFILE,
                              use.all=F,
                              choose.by.clustering=T)
} else {
  dryclean::prepare_detergent(normal.table.path=normal.file,
                              path.to.save=opt$outdir,
                              num.cores=opt$cores,
                              use.all=T,
                              save.pon=T,
                              build='hg38',
                              PAR.file=PARFILE)  
}

