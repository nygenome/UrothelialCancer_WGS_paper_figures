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
## Run mutation timing with MutationTimeR
libs = c('optparse', 'VariantAnnotation', 'GenomicRanges','MutationTimeR')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)

## Set seed for reproducibility
set.seed(8162021)


## Get arguments
option_list = list(
  make_option(c("-c", "--cn"),       type='character', help="Allele-specific copy number output RDS of init-mutation-timer-input.r"),
  make_option(c("-p", "--purity"),   type='numeric',   help="Sample purity"),
  make_option(c("-t", "--tumor"),    type='character', help="Tumor ID"),
  make_option(c("-v", "--vcf"),      type='character', help="Sample VCF"),
  make_option(c("-g", "--gender"),   type='character', help="Sample gender"),
  make_option(c("-o", "--out_file"), type='character', help="Output RDS"))
opt = parse_args(OptionParser(option_list=option_list))


## Read data
cnv = readRDS(opt$cn)
vcf = VariantAnnotation::readVcf(opt$vcf)


## Filter for high confidence variants, or high confidence that were
## added as part of the union pileup
vcf = vcf[info(vcf)$HighConfidence | info(vcf)$TYPE == 'ADDED']

info(vcf)$t_ref_count = sapply(geno(vcf)$AD[,opt$tumor], function(x) x[1])
info(vcf)$t_alt_count = sapply(geno(vcf)$AD[,opt$tumor], function(x) x[2])


## Exclude mutations marked 'PartOfMNV'
vcf = vcf[rowRanges(vcf)$FILTER == 'PASS']


## Update seqlevels style/ From a look at the MutationTimeR code it 
## expects no chr prefix
seqlevels(cnv) = gsub('^chr','',seqlevels(cnv))
seqlevels(vcf) = gsub('^chr','',seqlevels(vcf))


## Run MutationTimeR
mt = MutationTimeR::mutationTime(vcf=vcf,
                                 cn=cnv,
                                 purity=opt$purity,
                                 gender=opt$gender,
                                 n.boot=200)


## Write result
saveRDS(mt, opt$out_file)
message(opt$out_file)
