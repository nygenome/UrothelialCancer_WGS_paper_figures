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
## Report summary stats for ecDNAs
libs = c('optparse', 'gGnome', 'gUtils', 'GenomicRanges')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


alpha = function(col, alpha) {

  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)

}



## Read jabba footprints table
read.jabba.fp = function(f, sample_id) {

  if (!file.exists(f)) {
    return(NULL)
  }
  
  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  x$sample = sample_id

  return(x)

}



## Convert a footprint record into GRanges
fp2gr = function(x, id) {
  
  x = gsub('\\+(?=,)|\\-(?=,)','',x, perl=T)
  x = gsub('\\+(?=;)|\\-(?=;)','',x, perl=T)
  x = gsub('\\+$|\\-$','',x, perl=T)
  x = strsplit(x,';|,')
  
  res = list()
  
  for (i in 1:length(x)) {
    
    seqnames = gsub(':.*','',x[[i]])
    start = as.numeric(sapply(x[[i]], function(y) unlist(strsplit(gsub('.*:','',y),'-'))[1]))
    end = as.numeric(sapply(x[[i]], function(y) unlist(strsplit(gsub('.*:','',y),'-'))[2]))
    
    res[[i]] = GRanges(seqnames, IRanges(start=start, end=end))
  }
  
  res = do.call(c, res)
  res$id = id
  
  return(res) 
  
}



## Compute shannon entropy
## Normalize by dividing by maximum possible entropy 
shannon.entropy = function(p) {
  if (min(p) < 0 || sum(p) <= 0) return(NA)
  
  p.norm = p[p>0]/sum(p)
  h = -sum(log2(p.norm)*p.norm)
  h.max = log2(length(p.norm))

  return(h / h.max)

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir_jba"), type='character', help="Jabba output directory"),
  make_option(c("-m", "--metadata"),   type='character', help="Sample metadata"),
  make_option(c("-c", "--cgc_bed"),   type='character', help="BED file with CGC genes"),
  make_option(c("-f", "--ktg_flag"),   type='character', help="Flag to set SigProfilerClusters run type"),
  make_option(c("-t", "--tn_file"),    type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-o", "--out_file"),   type='character', help="Output file"))
opt = parse_args(OptionParser(option_list=option_list))



## Read CGC genes
cgc = makeGRangesFromDataFrame(read.csv(opt$cgc_bed, 
                                        h=F, 
                                        stringsAsFactors=F, 
                                        sep='\t', 
                                        col.names=c('chr','start','end','gene')), 
                                        keep.extra.columns=T)


## Read metadata
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
mta = mta[!is.na(mta$tumor), ]


## Read TN pairs and get input files
tn = read.csv(opt$tn, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)

tn$in_files_jba_fp = paste0(opt$in_dir_jba,'/',tn$pair_name,'/jabba.events.footprints.kataegis.',opt$ktg_flag,'.txt')
tn$in_files_jba_gg = paste0(opt$in_dir_jba,'/',tn$pair_name,'/jabba.events.rds')


## Read footprints and select ecDNA
jba = do.call(rbind, mapply(read.jabba.fp, f=tn$in_files_jba_fp, sample_id=tn$tumor, SIMPLIFY=F, USE.NAMES=F))
jba = jba[jba$final_type == 'ecDNA', ]


## Check to see if pre or post-chemo samples are enriched in ecDNA presence 
mta$has_ecdna = ifelse(mta$tumor %in% jba$sample, 'ecDNA+', 'ecDNA-')
mta$has_jabba = mta$tumor %in% tn$tumor

ecdna.chemo.summary = table(mta$has_ecdna[mta$has_jabba], mta$platinum_chemotherapy[mta$has_jabba])
ecdna.chemo.summary
fisher.test(ecdna.chemo.summary)



#####################
## Annotate ecDNAs ##
#####################

mta.cols = c('platinum_chemotherapy', 'patient')
jba[, mta.cols] = mta[match(jba$sample, mta$tumor), mta.cols]


for (i in 1:nrow(jba)) {

  ## Get footprint and kataegis event
  id = jba$sample[i]
  fp.i = fp2gr(x=jba$footprint[i], id=jba$id[i])

  ## Read genome graph
  gg = readRDS(tn$in_files_jba_gg[tn$tumor == id])

  ## Extract ecDNA footprint  
  gg = gg$trim(fp.i)

  ## Compute summary statistics
  jba$max_cn[i] = max(gg$nodes$dt$cn)
  jba$genomic_mass_mb[i] = sum(gg$nodes$gr$cn * width(gg$nodes$gr)) / 1E6
  jba$cgc_genes[i] = paste(cgc$gene[cgc %^% gg$nodes$gr], collapse=',')
  jba$footprint_mb[i] = sum(width(gg$nodes$gr)) / 1E6
  jba$mean_segment_size_mb[i] = mean(width(gg$nodes$gr)) / 1E6
  jba$sd_segment_size_mb[i] = sd(width(gg$nodes$gr)) / 1E6
  jba$num_junctions[i] = nrow(gg$edges[type == 'ALT']$dt)
  jba$max_jcn[i] = max(gg$edges[type == 'ALT']$dt$cn)
  jba$mean_jcn[i] = mean(gg$edges[type == 'ALT']$dt$cn)
  jba$jcn_shannon_entropy[i] = shannon.entropy(table(gg$edges[type == 'ALT']$dt$cn))
  jba$width_weighted_cn[i] = sum(gg$nodes$gr$cn * (width(gg$nodes$gr) / sum(width(gg$nodes$gr))))

}

col.sel = c('sample','patient','platinum_chemotherapy', 'type', 'final_type', 'has_kataegis', 
            'max_cn', 'num_junctions', 'max_jcn', 'mean_jcn', 'jcn_shannon_entropy', 'genomic_mass_mb', 
            'footprint_mb', 'mean_segment_size_mb', 'sd_segment_size_mb', 'width_weighted_cn',
            'cgc_genes')
jba = jba[, col.sel]


## Write result
write.table(jba, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)
