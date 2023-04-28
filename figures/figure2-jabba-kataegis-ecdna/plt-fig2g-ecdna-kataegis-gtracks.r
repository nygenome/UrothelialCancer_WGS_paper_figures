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
## Plot gTracks for ecDNA events with kataegis

libs = c('optparse', 'GenomicRanges', 'gGnome', 'gTrack', 'data.table')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))



## Read SigProfilerClusters annotated output
read.sigprofiler = function(f) {

  if (!file.exists(f)) {
    return(NULL)
  }

  x = fread(f, header=T, stringsAsFactors=F, sep='\t', data.table=F)
  x$end = x$pos
  x = makeGRangesFromDataFrame(x, 
                               seqnames.field='chrom',
                               start.field='pos',
                               end.field='end',
                               keep.extra.columns=T)

  return(x)

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



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir_jba"),   type='character', help="Jabba output directory"),
  make_option(c("-I", "--in_dir_ktg"),   type='character', help="Directory holding SigProfilerClusters kataegis VCFs"),
  make_option(c("-f", "--ktg_flag"),     type='character', help="Flag to set SigProfilerClusters run type"),
  make_option(c("-n", "--in_dir_nc"),    type='character', help="Directory holding SigProfilerClusters VCFs (unclustered)"),
  make_option(c("-t", "--tn_file"),      type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-o", "--out_file"),     type='character', help="PDF output"))
opt = parse_args(OptionParser(option_list=option_list))



## Read TN pairs and get input files
tn = read.csv(opt$tn, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)

tn$in_files_ktg = paste0(opt$in_dir_ktg,'/',tn$pair_name,'-mutationtimer-snv.annotated.jabba.txt')
tn$in_files_nc = paste0(opt$in_dir_nc,'/',tn$pair_name,'-mutationtimer-snv.annotated.jabba.txt')
tn$in_files_jba_fp = paste0(opt$in_dir_jba,'/',tn$pair_name,'/jabba.events.footprints.kataegis.',opt$ktg_flag,'.txt')
tn$in_files_jba_gg = paste0(opt$in_dir_jba,'/',tn$pair_name,'/jabba.events.rds')


## Read input 
ktg = Reduce(c, lapply(tn$in_files_ktg, read.sigprofiler))
nc = Reduce(c, lapply(tn$in_files_nc, read.sigprofiler))
jba = do.call(rbind, mapply(read.jabba.fp, f=tn$in_files_jba_fp, sample_id=tn$tumor, SIMPLIFY=F, USE.NAMES=F))


## Only plot ecDNA events with kataegis 
jba = jba[jba$final_type == 'ecDNA' & jba$has_kataegis, ]
jba$ecdna_id = make.unique(jba$sample)


ktg = ktg[ktg$ecdna]
ktg$assigned_signature = 'Kyklonas'
nc = nc[nc$jabba_event == 'ecDNA' & nc$assigned_signature %in% paste0('SBS', c(2,13,31,35))]

nc$assigned_signature[nc$assigned_signature %in% c('SBS2', 'SBS13')] = 'APOBEC (SBS2 + SBS13)'
nc$assigned_signature[nc$assigned_signature %in% c('SBS31', 'SBS35')] = 'Platinum chemotherapy (SBS31 + SBS35)'


colormap = c(`Kyklonas`='#2E3192', 
             `APOBEC (SBS2 + SBS13)`='#005BAA', 
             `Platinum chemotherapy (SBS31 + SBS35)`='#E21F26')

nc = nc[order(factor(nc$assigned_signature, levels=names(colormap)))]



###########################
## Plot each ecDNA event ##
###########################

pdf(opt$out_file, width=13, height=10)

for (i in 1:nrow(jba)) {

  ## Get footprint and kataegis event
  id = jba$sample[i]
  fp.i = fp2gr(x=jba$footprint[i], id=jba$id[i])
  ktg.i =  ktg[ktg$sample == id, ]
  nc.i =  nc[nc$sample == id, ]
  
  ## Read genome graph
  gg = readRDS(tn$in_files_jba_gg[tn$tumor == id])

  ## Color ecDNA 
  gg$edges[type=='ALT' & col == 'purple']$mark(col='#2EC3E7')
  gg$nodes[col == 'purple']$mark(col='#2EC3E7')
  gg$edges[type=='ALT' & is.na(col)]$mark(col='gray50')


  ## Expand footprint a little bit for some context 
  gg.subgraph = gg$copy$subgraph(fp.i, d=5E5)
  fp.i = gg.subgraph$footprint


  ## Build gTracks and plot
  gg.gt = gg$gtrack(name=id, ylab='CN', col.loose='gray50')
  ktg.gt = gTrack(ktg.i, 
                  circle=T, 
                  stack.gap=0, 
                  height=3,
                  colormaps=colormap['Kyklonas'],
                  y.field='vcf_hc_vaf',  
                  y0=0, y1=1,
                  yaxis.pretty=4,
                  gr.colorfield='assigned_signature',
                  name="Kataegis", 
                  ylab='VAF')

  nc.gt = gTrack(nc.i[nc.i %^% fp.i], 
                  circle=T, 
                  stack.gap=0, 
                  height=3,
                  colormaps=colormap,
                  y.field='vaf', 
                  y0=0, y1=1, 
                  yaxis.pretty=4,
                  gr.colorfield='assigned_signature',
                  name="Non-clustered", 
                  ylab='VAF',
                  xaxis.unit=1e6, 
                  xaxis.suffix='MB',
                  xaxis.chronly=TRUE,
                  xaxis.width=FALSE,
                  chr.sub=FALSE)

  plot(c(nc.gt, ktg.gt, gg.gt), fp.i)


  ## Write source data ####################################

  nodes = as.data.frame(gg.subgraph$nodes$dt)
  edges = as.data.frame(gg.subgraph$edges$dt)

  nodes.out.txt = gsub('\\.pdf$', paste0('.', jba$ecdna_id[i],'.nodes.txt'), opt$out_file)
  edges.out.txt = gsub('\\.pdf$', paste0('.', jba$ecdna_id[i],'.edges.txt'), opt$out_file)

  write.table(nodes, nodes.out.txt, row.names=F, col.names=T, quote=F, sep='\t')
  write.table(edges, edges.out.txt, row.names=F, col.names=T, quote=F, sep='\t')
  
  mut.columns = c('seqnames', 'start', 'assigned_signature', 'vaf')

  mutations = rbind(as.data.frame(ktg.i)[, mut.columns], as.data.frame(nc.i[nc.i %^% fp.i])[, mut.columns])
  mutations = mutations[!is.na(mutations$vaf), ]

  mutations.out.txt = gsub('\\.pdf$', paste0('.', jba$ecdna_id[i],'.mutations.txt'), opt$out_file)
  write.table(mutations, mutations.out.txt, row.names=F, col.names=T, quote=F, sep='\t')
  
  ########################################################


}

dev.off()
message(opt$out_file)
