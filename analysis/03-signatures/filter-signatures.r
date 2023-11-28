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
## Sanity check signature assignment
libs = c('optparse', 'VariantAnnotation', 'stringr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Parse signature probability assignments
## x: string of the form SBS1:numeric,SBS2:numeric,...
## min.diff: If there are multiple signature probabilities, what is the minimum
##           difference allowed between the top 2 signatures? This should help
##           control for cases where we can't unambiguously assign a signature
parse.signature.field = function(x, min.diff=0.2, min.prob) {
  
  ## Case where no signature was assigned
  if (length(x) == 0) {
    return(NA)
  }
  
  ## Parse string into named vector
  x = unlist(strsplit(x,':'))
  x = structure(as.numeric(x[1:length(x) %% 2 == 0]), 
                names=x[!1:length(x) %% 2 == 0])
  
  ## Filter for minimum probability
  x = x[x > min.prob]
  
  ## Case where none meet the minimum posterior probability
  if (length(x) == 0) {
    return(NA)
  }
  
  ## Case where only one signature was assigned
  if (length(x) == 1) {
    return(names(x))
  }  
  
  ## Case where multiple signatures assigned
  x = x[1:2]
  
  if (x[1] - x[2] > min.diff) {
    return(names(x)[1])
  } else {
    return(NA)
  }
  
  
}



vcf2signatures = function(f, sigs) {
  
  x = readInfo(file=f, x='chosen_signature')
  x = sapply(x, parse.signature.field)  
  
  res = table(factor(x, levels=sigs)) / length(x)
  return(res)
  
}



## Convert SBS reference to binary matrix by setting top X peaks to TRUE
sbs2bin = function(sbs, x) {
  
  ## For each signature
  for (i in 1:ncol(sbs)) {
    
    ## Find top x peaks and set them to TRUE
    j = order(sbs[, i], decreasing=T)[1:x]
    
    sbs[, i] = FALSE
    sbs[j, i] = TRUE
    
  }
  
  return(sbs)
  
}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),          type='character', help="Input VCF with info fields signature_posterior and tricontext"),
  make_option(c("-s", "--sbs_ref"),          type='character', help="COSMIC SBS reference"),
  make_option(c("-f", "--sbs_31_35_filter"), type='character', help="COSMIC SBS reference"),
  make_option(c("-c", "--context_cutoff"),   type='numeric',   help="Number of top SBS reference peaks to retain", default=5),
  make_option(c("-p", "--min_probability"),  type='numeric',   help="Minimum posterior probability", default=0.5),
  make_option(c("-m", "--metadata"),         type='character', help="Metadata file with platinum treatment and sample collection info"),
  make_option(c("-t", "--tumor_id"),         type='character', help="Sample tumor ID", default=0.5),
  make_option(c("-o", "--out_file"),         type='character', help="Output VCF"))
opt = parse_args(OptionParser(option_list=option_list))



## Read input VCF
vcf = readVcf(opt$in_file)

## Read metadata, select data for this tumor sample
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
mta = mta[!is.na(mta$tumor), ]
mta$platinum_chemotherapy = mta$platinum_chemotherapy == 'Post-chemo'
mta = mta[mta$tumor == opt$tumor_id, ]


## Read cosmic reference 
sbs.ref = read.csv(opt$sbs_ref, h=T, stringsAsFactors=F, sep='\t')
rownames(sbs.ref) = sbs.ref$Type
sbs.ref = sbs.ref[, -1]

## Convert to binary matrix that tells us whether to keep a 
## particular context/signature assignment
sbs.ref.bin = sbs2bin(sbs.ref, opt$context_cutoff)


## Chose a signature if possible
info(vcf)$assigned_signature = sapply(info(vcf)$signature_posterior, parse.signature.field, min.prob=opt$min_probability)


## Filter for mutations characteristic of each signature (i.e. don't 
## keep assignments for low peaks in a given signatures mutational spectra)
sig.keep = unlist(mapply(function(i, j) ifelse(is.na(j), F, sbs.ref.bin[i, j]), 
                                                         i=info(vcf)$tricontext, 
                                                         j=info(vcf)$assigned_signature))

## Special rules for SBS2/SBS13 (APOBEC) and SBS31/SBS35 (platinum chemotherapy)
allowed.substitutions = list(SBS2=c('C>T'), 
                             SBS13=c('C>A', 'C>G'), 
                             SBS31=c('C>T', 'T>A'), 
                             SBS35=c('C>A', 'C>G', 'C>T', 'T>A'))

substitutions = str_extract(unlist(info(vcf)$tricontext), pattern='(?<=\\[).*(?=\\])')
sig.keep.special = unlist(mapply(function(i, j) ifelse(is.na(j), F, i %in% allowed.substitutions[[j]]), 
                                   i=substitutions, 
                                   j=info(vcf)$assigned_signature))

## Take sig.keep, unless the assigned signature should be manually curated. In those cases, use the manual curation (sig.keep.special)
sig.keep = sig.keep.special | (sig.keep & !is.na(info(vcf)$assigned_signature) & !info(vcf)$assigned_signature %in% names(allowed.substitutions))
info(vcf)$assigned_signature[!sig.keep] = NA


## Additionally filter a list of SBS31/SBS35 mutations
if (!is.null(opt$sbs_31_35_filter)) {
  
  bl = read.csv(opt$sbs_31_35_filter, F, stringsAsFactors=T, sep='\t', col.names=c('chr','pos','ref','alt'))
  bl$key = paste0(bl$chr, ':', bl$pos, '_', bl$ref, '/', bl$alt)
  info(vcf)$in_sbs_31_35_filter_list = !is.na(info(vcf)$assigned_signature) & 
                                        info(vcf)$assigned_signature %in% c('SBS31','SBS35') & 
                                        names(rowRanges(vcf)) %in% bl$key
  
  info(vcf)$assigned_signature[info(vcf)$in_sbs_31_35_filter_list] = NA
  
}


if (nrow(mta) > 0 && !mta$platinum_chemotherapy) {
  n.filtered = sum(info(vcf)$assigned_signature %in% c('SBS31','SBS35'))
  message('Filtering out ', n.filtered, ' SBS31/SBS35 mutations based on treatment info')
  info(vcf)$assigned_signature[info(vcf)$assigned_signature %in% c('SBS31','SBS35')] = NA
}


## Update header
sig.header = DataFrame(Number=c('1','0'), 
                       Type=c('String', 'Flag'),
                       Description=c('Assigned signature','SBS31/SBS35-assigned variant found in filter list'), 
                       row.names=c('assigned_signature', 'in_sbs_31_35_filter_list'))
info(header(vcf)) = rbind(info(header(vcf)), sig.header)


## Write result
writeVcf(vcf, opt$out_file)
message(opt$out_file)
