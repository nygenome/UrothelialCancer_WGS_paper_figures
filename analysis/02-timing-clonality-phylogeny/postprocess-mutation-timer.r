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
## Postprocess mutation timing results
libs = c('optparse', 'VariantAnnotation', 'GenomicRanges', 'MutationTimeR', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Workaround for GenomeInfoDb issue
to.ucsc = function(x) {
  
  missing.chr = grep('^chr', x, fixed=T, invert=T)
  x[missing.chr] = paste0('chr', x[missing.chr])
  
  return(x)
  
}



## Compute a CCF using mutation timer multiplicity
add.mt.cff = function(vcf, purity, gender) {
  
  ncn = structure(rep(2, 23), names=c(1:22, 'X'))
  
  ## Check gender 
  if (!gender %in% c('male','female')) {
    stop('Gender should be male or female. Given: ', gender)
  }
    
  if (gender == 'male') {
    ncn['chrX'] = 1
  }
  
  ## f = p*m*phi / psi
  locus.total.cn = info(vcf)$MajCN + info(vcf)$MinCN
  avg.ploidy = (purity * locus.total.cn) + ((1 - purity) * ncn[as.character(seqnames(rowRanges(vcf)))]) 
  vaf = info(vcf)$t_alt_count / (info(vcf)$t_alt_count + info(vcf)$t_ref_count)
 
  info(vcf)$MT_CCF = (vaf * avg.ploidy) / (purity * info(vcf)$MutCN)
  
  ## Cap at 1 (clonal mutations will be distributed about 1)
  info(vcf)$MT_CCF[info(vcf)$MT_CCF > 1] = 1 
  
  return(vcf)
    
}



## Get arguments
option_list = list(
  make_option(c("-m", "--mt"),               type='character', help="MutationTimeR result RDS"),
  make_option(c("-c", "--cn"),               type='character', help="Allele-specific copy number output RDS of init-mutation-timer-input.r"),
  make_option(c("-u", "--purity"),           type='numeric',   help="Sample purity"),
  make_option(c("-g", "--gender"),           type='character', help="Sample gender"),
  make_option(c("-v", "--vcf"),              type='character', help="Sample VCF"),
  make_option(c("-t", "--tumor"),            type='character', help="Tumor ID"),
  make_option(c("-b", "--out_file_bed"),     type='character', help="Output BED annotated with MutationTimeR results"),
  make_option(c("-o", "--out_file_vcf"),     type='character', help="Output VCF annotated with MutationTimeR results"),
  make_option(c("-s", "--out_file_summary"), type='character', help="Summary txt file"),
  make_option(c("-p", "--out_file_pdf"),     type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



## Read data
mt = readRDS(opt$mt)
vcf = VariantAnnotation::readVcf(opt$vcf)


## Filter for high confidence variants, or high confidence that were
## added as part of the union pileup
vcf = vcf[info(vcf)$HighConfidence | info(vcf)$TYPE == 'ADDED']



info(vcf)$t_ref_count = sapply(geno(vcf)$AD[,opt$tumor], function(x) x[1])
info(vcf)$t_alt_count = sapply(geno(vcf)$AD[,opt$tumor], function(x) x[2])

cnv = readRDS(opt$cn)


## Make sure seqinfos are consistent
seqlengths(cnv) = seqlengths(vcf)[names(seqlengths(cnv))]
seqlevels(cnv) = gsub('^chr','',seqlevels(cnv))
seqlevels(vcf) = gsub('^chr','',seqlevels(vcf))


## Annotate bed/vcf
mcols(cnv) = cbind(mcols(cnv), mt$T)

vcf = MutationTimeR::addMutTime(vcf, mt$V)

vcf = add.mt.cff(vcf, purity=opt$purity, gender=opt$gender)



## Write annotated results
cnv.out = cnv
seqlevels(cnv.out) = to.ucsc(seqlevels(cnv.out))

cnv.out = as.data.frame(cnv.out)
colnames(cnv.out)[colnames(cnv.out) == 'seqnames'] = '#chr'
cnv.out = cnv.out[, colnames(cnv.out) != 'timing_param'] 
write.table(cnv.out, opt$out_file_bed, row.names=F, col.names=T, sep='\t', quote=F)

vcf.out = vcf
seqlevels(vcf.out) = to.ucsc(seqlevels(vcf.out))


## Add header documenting CCF
vcf.out.header = DataFrame(Number=c('1'), 
                       Type=c('Float'),
                       Description=c(paste0('CCF computed using MutationTimeR multiplicity, assuming standard ', opt$gender,' normal copy configuration')), 
                       row.names=c('MT_CCF'))
info(header(vcf.out)) = rbind(info(header(vcf.out)), vcf.out.header)

## Write VCF
writeVcf(vcf.out, opt$out_file_vcf)


## Write summary
result.summary = as.data.frame(table(info(vcf)$CLS), stringsAsFactors=F)
colnames(result.summary) = c('clonality','count')
result.summary$sample = opt$tumor
result.summary = result.summary[, c('sample', 'clonality', 'count')]
write.table(result.summary, opt$out_file_summary, row.names=F, col.names=T, quote=F, sep='\t')


## Plot 
png(opt$out_file_pdf, width=1200, height=1200)
MutationTimeR::plotSample(vcf=vcf, cn=cnv, regions=si2gr(seqinfo(cnv)))
dev.off()
