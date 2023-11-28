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
## Add an additional layer of filtering for patients with chemotherapy-treated samples
libs = c('optparse', 'VariantAnnotation')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir"),     type='character', help="Input directory with VCF input"),
  make_option(c("-p", "--patient"),    type='character', help="Patient"),
  make_option(c("-t", "--tn_file"),    type='character', help="Tumor-normal pairing file"),
  make_option(c("-o", "--out_dir"),    type='character', help="Output dir"),
  make_option(c("-s", "--in_suffix"),  type='character', help="Input suffix"),
  make_option(c("-x", "--out_suffix"), type='character', help="Output suffix"), 
  make_option(c("-m", "--metadata"),   type='character', help="Sample metadata"))
opt = parse_args(OptionParser(option_list=option_list))


## Read TN file, select patient of interest 
tn = read.csv(opt$tn_file, h=F, sep='\t', stringsAsFactors=F, col.names=c('tumor','normal','gender','patient'))
tn = tn[tn$patient == opt$patient, ]
tn$pair_name = paste0(tn$tumor, '--', tn$normal)
tn$in_file = paste0(opt$in_dir, '/', tn$pair_name, opt$in_suffix)
tn$out_file = paste0(opt$out_dir, '/', tn$pair_name, opt$out_suffix)

if (!all(file.exists(tn$in_file))) {
  stop('Missing files!')
}

## Read metadata
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
mta = mta[!is.na(mta$tumor), ]
mta$platinum_chemotherapy = mta$platinum_chemotherapy == 'Post-chemo'
tn$platinum = mta$platinum_chemotherapy[match(tn$tumor, mta$tumor)]


## If none of the samples received platinum, we don't need to 
## do any additional filtering
if (all(!tn$platinum)) {
  
  message('No platinum-treated samples. Skipping filtering...')
  cmd = paste('cp', tn$in_file, tn$out_file)
  print(sapply(cmd, system, intern=T))
  
}


## If all of the samples received platinum, a mutation assigned to 
## platinum in one sample is enough to "rescue" the same mutation in other samples
if (all(tn$platinum)) {
  
  message('All samples received platinum therapy')
  
  ## Read all variants
  vcf = lapply(tn$in_file, VariantAnnotation::readVcf)  
  
  ## Select all platinum-assigned mutations  
  mut.31 = lapply(vcf, function(x) names(rowRanges(x)[!is.na(info(x)$assigned_signature) & 
                                              info(x)$assigned_signature %in% c('SBS31')]))
  mut.31 = unique(unlist(mut.31))
  
  mut.35 = lapply(vcf, function(x) names(rowRanges(x)[!is.na(info(x)$assigned_signature) & 
                                                        info(x)$assigned_signature %in% c('SBS35')]))
  mut.35 = unique(unlist(mut.31))
  
  
  ## Go back and "rescue" mutations that weren't previously assigned to platinum
  for (i in 1:length(vcf)) {
    
    rescue.31 = (is.na(info(vcf[[i]])$assigned_signature) | 
                 !info(vcf[[i]])$assigned_signature %in% c('SBS31', 'SBS35')) &
                 names(rowRanges(vcf[[i]])) %in% mut.31
    
    rescue.35 = (is.na(info(vcf[[i]])$assigned_signature) | 
                   !info(vcf[[i]])$assigned_signature %in% c('SBS31', 'SBS35')) &
                    names(rowRanges(vcf[[i]])) %in% mut.35
    
    info(vcf[[i]])$assigned_signature[rescue.31] = 'SBS31'
    info(vcf[[i]])$assigned_signature[rescue.35] = 'SBS35'
    
    message('Rescued ', sum(rescue.31 | rescue.35),' mutations for ', tn$tumor[i])
    
  }
  
  
  ## Write result
  for (i in 1:length(vcf)) {
     
    VariantAnnotation::writeVcf(vcf[[i]], tn$out_file[i])
    
  }
  
}


## If only some of the samples received chemo, compare untreated vs treated variants
if (!all(tn$platinum) && any(tn$platinum)) {
  
  message('Only some samples received platinum therapy')
  
  vcf = lapply(tn$in_file, VariantAnnotation::readVcf)
  
  ## Select all non-chemo mutations
  mut = lapply(vcf[!tn$platinum], function(x) names(rowRanges(x)))
  mut = unique(unlist(mut))
  
  for (i in which(tn$platinum)) {
    
    filter.mut = !is.na(info(vcf[[i]])$assigned_signature) & 
                  info(vcf[[i]])$assigned_signature %in% c('SBS31', 'SBS35') & 
                  names(rowRanges(vcf[[i]])) %in% mut
    
    
    info(vcf[[i]])$assigned_signature[filter.mut] = NA
    
    message('Filtered ', sum(filter.mut),' mutations for ', tn$tumor[i])
    
  }
  
  
  ## Write result
  for (i in 1:length(vcf)) {
    
    VariantAnnotation::writeVcf(vcf[[i]], tn$out_file[i])
    
  }
  
}
