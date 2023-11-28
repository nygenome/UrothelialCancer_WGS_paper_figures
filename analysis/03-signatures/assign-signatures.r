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
## For each SNV, assign posterior probabilities of belonging to each signature
libs = c('optparse', 'VariantAnnotation', 'BSgenome.Hsapiens.UCSC.hg38', 'data.table')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Convenience functions for accessing ref and alt alleles
## vcf: VariantAnnotaton VCF object
refc = function(vcf) as.character(VariantAnnotation::ref(vcf))
altc = function(vcf) as.character(unlist(VariantAnnotation::alt(vcf)))



## Convert to trinucleotide contexts in the same manner as deconstructSigs
## Adapted from https://github.com/raerose01/deconstructSigs/blob/master/R/mut.to.sigs.input.R
## vcf: VariantAnnotaton VCF object
add.trinuc.context = function(vcf) {
  
  bsg = BSgenome.Hsapiens.UCSC.hg38
  info(vcf)$context = BSgenome::getSeq(bsg, as.character(seqnames(rowRanges(vcf))), start(rowRanges(vcf))-1, end(rowRanges(vcf))+1, as.character=T) 
  info(vcf)$mutcat = paste0(info(vcf)$ref, ">", info(vcf)$alt)
  
  ## Check to make sure REF alleles match what we got from bsgenome
  if(any(info(vcf)$ref != substr(info(vcf)$context, 2, 2))){
    stop('Check ref bases -- not all match context')
  }
  
  # Reverse complement the G's and A's
  gind = grep("G",substr(info(vcf)$mutcat,1,1))
  tind = grep("A",substr(info(vcf)$mutcat,1,1))
  
  info(vcf)$std.mutcat = info(vcf)$mutcat
  info(vcf)$std.mutcat[c(gind, tind)] = gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", info(vcf)$std.mutcat[c(gind, tind)])))) # to lowercase
  info(vcf)$std.mutcat[c(gind, tind)] = gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", info(vcf)$std.mutcat[c(gind, tind)])))) # complement
  
  info(vcf)$std.context = info(vcf)$context
  info(vcf)$std.context[c(gind, tind)] = gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", info(vcf)$std.context[c(gind, tind)])))) # to lowercase
  info(vcf)$std.context[c(gind, tind)] = gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", info(vcf)$std.context[c(gind, tind)])))) # complement
  info(vcf)$std.context[c(gind, tind)] = sapply(strsplit(info(vcf)$std.context[c(gind, tind)], split = ""), function(str) {paste(rev(str), collapse = "")}) # reverse
  
  # Make the tricontext
  info(vcf)$tricontext = paste(substr(info(vcf)$std.context, 1, 1), "[", info(vcf)$std.mutcat, "]", substr(info(vcf)$std.context, 3, 3), sep = "")
  
  return(vcf)  
}



## Annotate each mutation with a probability of belonging to each signature
## Adapted from code originally written by Marcin Imielinski's group
## vcf: VariantAnnotaton VCF object with trinucleotide context annotated in info(vcf)$tricontext
## sig.ref: data.frame with signature reference data. Should match the reference used to generate sig.cont. Rows are trinucleotide contexts,
##          i.e. C[C>T]G and columns are signature IDs
## sig.cont: Signature contributions for this sample. Should be a named vector with names corresponding to signature IDs
## cutoff: Minimum probability for inclusion in output VCF
## precision: Rounding precision for output probabilities
add.posteriors = function(vcf, sig.ref, sig.cont, cutoff, precision=3) {
  
  ## Numerator in bayes rule
  p_intermat = t(apply(sig.ref, 1, function(x) x*sig.cont))
  t_dt = as.data.frame(p_intermat)
  t_dt[, "mut_ind"] = rownames(t_dt)
  t_dt = as.data.table(t_dt)
  setkey(t_dt, mut_ind)

    
  ## 1/p_mi is the denominator in bayes rule
  p_mi = as.matrix(sig.ref) %*% sig.cont
  c = rownames(p_mi)
  p_mi = as.data.table(p_mi)
  p_mi[, mut_ind2 := c]
  
  setnames(p_mi, "V1", "p_mi")
  p_mi[, inv_p_mi := 1/p_mi]
  setkey(p_mi, mut_ind2)
  
  
  ## Per-mutation posteriors
  info(vcf)$inv_p = p_mi[info(vcf)$tricontext, inv_p_mi]
  p_mi_sig = as.data.frame(info(vcf)$inv_p * t_dt[info(vcf)$tricontext, .SD, .SDcols = 1:length(sig.cont)])
  
  
  ## Build info field, only including signatures with probability above cutoff
  info(vcf)$signature_posterior = ''
  for (i in 1:nrow(p_mi_sig)) {
    
    keep = which(p_mi_sig[i, , drop=T] > cutoff)
    
    if (length(keep) > 0) {
      sig.name.keep = colnames(p_mi_sig)[keep]
      sig.prob.keep = round(unname(unlist(p_mi_sig[i, keep])), precision)
      
      ord = order(sig.prob.keep,decreasing=T)
      
      info(vcf)$signature_posterior[i] = paste(paste0(sig.name.keep[ord],':',sig.prob.keep[ord]), collapse=',')  
    } 
    
  }
    
  return(vcf)
  
}



## Get arguments
option_list = list(
  make_option(c("-i", "--vcf"),       type='character', help="Input VCF (should have already run through MutationTimer"),
  make_option(c("-t", "--sample_id"), type='character', help="Sample ID corresponding to the entry in the signature results Sample column"),
  make_option(c("-s", "--sbs"),       type='character', help="COSMIC SBS contributions"),
  make_option(c("-r", "--sbs_ref"),   type='character', help="COSMIC SBS reference, should be the same file used to generate sample-level results"),
  make_option(c("-c", "--cutoff"),    type='numeric',   help="Minimum posterior probability for including a signature in the output VCF", default=0.25),
  make_option(c("-o", "--out_file"),  type='character', help="Output VCF"))
opt = parse_args(OptionParser(option_list=option_list))


message('WARNING: THIS SCRIPT IS ONLY COMPATIBLE WITH GRCH38 FOR NOW')


## Read data
vcf = readVcf(opt$vcf)

## Filter for high confidence variants, or high confidence that were
## added as part of the union pileup
vcf = vcf[info(vcf)$HighConfidence | info(vcf)$TYPE == 'ADDED']

## Read SBS results and select sample of interest
sbs = read.csv(opt$sbs, h=T, stringsAsFactors=F, sep='\t')
sbs = sbs[sbs$Sample == opt$sample_id, setdiff(colnames(sbs), 'Sample')]

## Read SBS reference data
sbs.ref = read.csv(opt$sbs_ref, h=T, stringsAsFactors=F, sep='\t')
rownames(sbs.ref) = sbs.ref$Type
sbs.ref = sbs.ref[, setdiff(colnames(sbs.ref), 'Type')]

## Make sure columns are aligned
sbs = unlist(sbs[,colnames(sbs.ref)])


## Subset to SNVs
vcf = vcf[sapply(alt(vcf), length) == 1]
info(vcf)$ref = refc(vcf)
info(vcf)$alt = altc(vcf)
vcf = vcf[nchar(info(vcf)$ref) == 1 & nchar(info(vcf)$alt) == 1]


## Add trinucleotide context, followed by posterior probabilities
vcf = add.trinuc.context(vcf)
vcf = add.posteriors(vcf=vcf, sig.ref=sbs.ref, sig.cont=sbs, cutoff=opt$cutoff)


## Remove some intermediate INFO fields, add header for the ones we're keeping
info(vcf) = info(vcf)[, setdiff(colnames(info(vcf)), c('context','mutcat','std.context','std.mutcat', 'ref', 'alt', 'inv_p'))]

sig.header = DataFrame(Number=c('.', '.'), 
                       Type=c('String', 'String'),
                       Description=c('Normalized trinucleotide context',
                                      'Posterior probability that the mutation came from a specific COSMIC mutational signature'),
                       row.names=c('tricontext', 'signature_posterior'))

info(header(vcf)) = rbind(info(header(vcf)), sig.header)


## Write result
writeVcf(vcf, opt$out_file)
