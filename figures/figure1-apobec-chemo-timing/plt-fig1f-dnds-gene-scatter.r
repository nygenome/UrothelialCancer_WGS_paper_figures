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
## Plot global dN/dS values, overlay apobec and count information

libs = c('optparse', 'VariantAnnotation', 'reshape2', 'ggplot2',  'ggrepel', 'patchwork')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Parse VEP annotations from a vcf object (default field is CSQ)
parse.vep = function(vcf, field='CSQ') { 
  
  ## Get format 
  frm = info(header(vcf))[field, 'Description']
  frm = gsub('^.*Format: ', '', frm)  
  frm = unlist(strsplit(frm, '|', fixed=T))  
  
  ## Read field
  colnames(info(vcf))
  x = unlist(info(vcf)[, field])
  
  ## Split into data frame
  x = colsplit(x, pattern='\\|', names=frm)
  
  x$IMPACT = trimws(x$IMPACT)
  
  return(x)
  
}



## Read VCF and extract mutations/signatures in genes of interest
read.vcf = function(f, genes) {

  ## Filter for high confidence
  x = VariantAnnotation::readVcf(f)
  x = x[info(x)$HighConfidence | info(x)$TYPE == 'ADDED']

  ## Fill in dummy signature column if needed
  if (!'assigned_signature' %in% colnames(info(x))) {
    info(x)$assigned_signature = 'Unassigned'
  }

  ## Filter for SNVs
  x = x[nchar(as.character(rowRanges(x)$REF)) == 1 & 
       nchar(as.character(unlist(rowRanges(x)$ALT))) == 1] 

  info(x)$assigned_signature[is.na(info(x)$assigned_signature)] = 'Unassigned'

  ## Parse VEP annotations
  vep = parse.vep(x)

  ## Whitelist mis-annotated HRAS based on COSMIC IDs in Existing variation
  ## checked by hand
  hras.pos = c(532749, 533874, 533875, 534289)
  hras.chr = 'chr11'
  vep$SYMBOL[as.character(seqnames(rowRanges(x))) == hras.chr & start(rowRanges(x)) %in% hras.pos] = 'HRAS'
  vep$IMPACT[as.character(seqnames(rowRanges(x))) == hras.chr & start(rowRanges(x)) %in% hras.pos] = 'MODERATE'

  ## Filter for genes of interest
  calls.keep = vep$SYMBOL %in% genes & vep$IMPACT %in% c('MODERATE', 'HIGH')
  x = x[calls.keep]
  vep = vep[calls.keep, ]

  ## Simplify output
  res = data.frame(gene=vep$SYMBOL, signature=info(x)$assigned_signature)

  return(res)

}



## Get arguments
option_list = list(
  make_option(c("-t", "--tn_file"),          type='character', help="Tumor normal pairs file"),
  make_option(c("-S", "--sig_genes_file"),   type='character', help="CSV with dndscv significant genes"),
  make_option(c("-i", "--in_dir_sig"),       type='character', help="Directory with signature-assigned VCFs"),
  make_option(c("-u", "--in_dir_union"),     type='character', help="Directory with patient union VCFs"),
  make_option(c("-s", "--in_dir_singleton"), type='character', help="Directory with singleton VCFs (e.g singletons with low purity)"),
  make_option(c("-o", "--out_file"),         type='character', help="Output SVG"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tn pairs and significant genes
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$tumor = gsub('^#', '', tn$tumor)
tn$pair_name = paste0(tn$tumor,'--',tn$normal)

sig.genes = read.csv(opt$sig_genes_file, h=T, stringsAsFactors=F)



## Construct file names, prefering signature-assigned VCFs if they exist, 
## patient-union if those don't, and finally FFPE-filtered singleton VCFs
## if neither of the first two exist 
tn$in_file = paste0(opt$in_dir_sig,'/',tn$pair_name,'.mutationtimer-snv.withSigProbs.filtered.vcf')
for (i in 1:nrow(tn)) {

  ## Try union
  if (!file.exists(tn$in_file[i])) {
    tn$in_file[i] = paste0(opt$in_dir_union,'/',tn$patient[i],'/',tn$pair_name[i],'.highconf.ffpefilter.union.vcf')
  }

  ## Try single-sample VCF
  if (!file.exists(tn$in_file[i])) {
    tn$in_file[i] = paste0(opt$in_dir_singleton,'/',tn$pair_name[i],'/',tn$pair_name[i],'.snv.indel.final.v6.annotated.ffpePON_filtered.vcf')
  }

}

## We shouldn't have any missing files by this point 
if (any(!file.exists(tn$in_file))) {
  stop('Missing files!')
}



## Read input files
vcf = lapply(tn$in_file, function(f) read.vcf(f, genes=sig.genes$gene_name))
vcf = do.call(rbind, vcf)


## Tabulate
vcf$gene = factor(vcf$gene, levels=sig.genes$gene_name)
res = as.data.frame.matrix(table(vcf$gene, vcf$signature))
res$gene = rownames(res)


## Reshape for plotting
res = melt(res, id.var='gene', variable.name='COSMIC SBS Signature', value.name='SNV Count')


## Format signatures
colors = c(`Unassigned`='grey90',
           `Other SBS`='grey50',
           `SBS31`='#e41a1c', 
           `SBS35`='#b2182b',
           `SBS2`='#377eb8', 
           `SBS13`='#2166ac')
           

res$`COSMIC SBS Signature` = as.character(res$`COSMIC SBS Signature`)
res$`COSMIC SBS Signature`[!res$`COSMIC SBS Signature` %in% names(colors)] = 'Other SBS'
res$`COSMIC SBS Signature` = factor(res$`COSMIC SBS Signature`, levels=names(colors))


## Sort genes by mutation count 
ct = tapply(res$`SNV Count`, res$gene, sum)
res$gene = factor(res$gene, levels=names(ct)[order(ct, decreasing=T)])


## Summarize by proportion SBS2/SBS13 and total count 
sig.genes$total_snv = tapply(res$`SNV Count`, res$gene, sum)[sig.genes$gene_name]
sig.genes$total_apobec = tapply(res$`SNV Count`[res$`COSMIC SBS Signature` %in% c('SBS2', 'SBS13')], res$gene[res$`COSMIC SBS Signature` %in% c('SBS2', 'SBS13')], sum)[sig.genes$gene_name]
sig.genes$prop_apobec = sig.genes$total_apobec / sig.genes$total_snv

## Temporary, while not using full TN pairs list 
sig.genes$prop_apobec[is.na(sig.genes$prop_apobec)] = 0

## Bin apobec proportion
sig.genes$`APOBEC SNV Proportion` = cut(sig.genes$prop_apobec, c(0,0.25, 0.5, 0.75,1), include.lowest=T)
sig.genes.col = structure(c('#87898a', '#9abbd6', '#569cd6', '#2166ac'), names=levels(sig.genes$`APOBEC SNV Proportion`))



##########
## Plot ##
##########

svg(opt$out_file, width=9, height=6)

## Base plot object
sig.genes$`Total Coding SNV Count` = sig.genes$total_snv
plt = ggplot(sig.genes, aes(x=wmis_cv, y=wnon_cv, col=`APOBEC SNV Proportion`, size=`Total Coding SNV Count`)) +
        geom_point() +
        scale_color_manual(values=sig.genes.col) +
        geom_text_repel(aes(label=gene_name),
                            col='black',
                            nudge_x=0.15,
                            box.padding=0.5,
                            nudge_y=1,
                            segment.alpha=0.3,
                            segment.color = NA,
                            size=4) +
        xlab('Missense dN/dS') +
        ylab('Nonsense dN/dS') +
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              text=element_text(size=15)) 


## Larger, main plot
p1 = plt + 
       xlim(0,40) + 
       ylim(0,60) +
       guides(size='legend')


## Use smaller inset for TP53
tp53.idx = which(sig.genes$gene_name == 'TP53')
pad = 5

p2 = plt + 
       scale_x_continuous(breaks=round(sig.genes$wmis_cv[tp53.idx]), 
                          label=round(sig.genes$wmis_cv[tp53.idx]), 
                          limits=c(round(sig.genes$wmis_cv[tp53.idx]) - pad, round(sig.genes$wmis_cv[tp53.idx]) + pad)) +
       scale_y_continuous(breaks=round(sig.genes$wnon_cv[tp53.idx]), 
                          label=round(sig.genes$wnon_cv[tp53.idx]), 
                          limits=c(round(sig.genes$wnon_cv[tp53.idx]) - pad, round(sig.genes$wnon_cv[tp53.idx]) + pad)) +                          
       theme(legend.position='none',
             plot.background=element_rect(fill='transparent', color=NA),
             plot.margin = unit(c(0, 0, 0, 0), "cm"),
             axis.title.x=element_blank(),
             axis.title.y=element_blank())


## Draw plot
p1 + inset_element(p2, 0.8, 0.8, 1, 1)

dev.off()
message(opt$out_file)


## Write source data 
out.file.txt = gsub('\\.svg$', '.txt', opt$out_file)
write.table(sig.genes[, -1], out.file.txt, row.names=F, col.names=T, quote=F, sep='\t')
message(out.file.txt)