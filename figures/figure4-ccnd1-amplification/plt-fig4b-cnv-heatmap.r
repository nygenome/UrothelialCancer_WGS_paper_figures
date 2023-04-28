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
## Plot CNV heatmap given a bed file with intervals of interest 

libs = c('optparse', 'GenomicRanges', 'gUtils', 'ComplexHeatmap', 'circlize')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Read CNVs 
read.cnv = function(f, id) {
  
  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t', col.names=c('chr','start','end', 'cn', 'gene'))
  x[,id] = x$cn
  x = x[,c('chr', 'start', 'end', id)]
    
  x = makeGRangesFromDataFrame(x, keep.extra.columns=T)
    
}



## Build color map from metadata fields
build.colormap = function(var, col) {
  
  var.col = structure(col, names=var)
  var.col = na.omit(var.col[!duplicated(var)])
  return(var.col)
  
}



## https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
numeric.mode = function(x, na.rm) {
  
  if (na.rm == TRUE) {

    x = x[!is.na(x)]

  }
  
  ux = unique(x)
  ux = ux[which.max(tabulate(match(x, ux)))]

  return(ux)

}



MUT_COL = c(`None`='white', `Moderate Impact`='#0000005c', `High Impact`='#000000B3')

ecdna.colors = c(`ecDNA-`='#FFFFFF', `ecDNA+`='aquamarine3')

events.long = c('Simple (TRA/INV)', 'Deletion', 'Duplication', 'Chromoplexy', 'Rigma', 'Pyrgo', 'Chromothripsis', 'TIC', 'QRP', 'BFB', 'Double Minute', 'Complex DM', 'Tyfonas','None')
events.short = c('simple', 'del', 'dup', 'chromoplexy', 'rigma', 'pyrgo', 'chromothripsis', 'tic', 'qrp', 'bfb', 'dm', 'cpxdm','tyfonas', 'None')
names(events.long) = events.short



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir"),      type='character', help="Directory holding JaBbA results"),
  make_option(c("-t", "--tn_file"),          type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-p", "--pp"),          type='character', help="Purity/ploidy file"),
  make_option(c("-m", "--metadata"),    type='character', help="Metadata file"),
  make_option(c("-r", "--bed"),         type='character', help="Regions to plot"),
  make_option(c("-g", "--gene_status"), type='character', help="Somatic mutation status of RB1"),
  make_option(c("-a", "--amp_status"), type='character', help="Amplification status of CCND1"),
  make_option(c("-b", "--binsize"),     type='numeric',   help="Size of bins to use. Large == more coarse grained", default=1E6),
  make_option(c("-o", "--out_file"),    type='character', help="Figure output file"))
opt = parse_args(OptionParser(option_list=option_list))


## Read TN pairs
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)
tn$name = gsub('-(Case-WGS|WGS)','', tn$tumor)

## Read purity/ploidy, align with TN pairs
pp = read.csv(opt$pp, h=F, stringsAsFactors=F, col.names=c('tumor','purity','ploidy'))
tn[,c('purity','ploidy')] = pp[match(tn$tumor, pp$tumor), c('purity','ploidy')]

## Treat XY ploidy differently for males
tn$ploidy_xy = tn$ploidy * ifelse(tn$gender=='male', 0.5, 1)

## Read in metadata, align with TN pairs
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
mta = mta[match(tn$tumor, mta$tumor), ]
mta$sample = paste0(mta$tumor,'--',mta$normal)

mta$Patient = mta$patient
mta$patient[mta$patient_color == '#FFFFFF'] = 'Singleton'
mta$platinum_chemotherapy = factor(mta$platinum_chemotherapy, levels=c('Pre-chemo', 'Post-chemo'))
mta$localization = factor(mta$localization, levels=c('Primary', 'Metastatic'))


## Make sure all expected files are present
tn$in.files = paste0(opt$in_dir,'/',tn$pair_name,'/jabba.simple.cnv.bed')
if (any(!file.exists(tn$in.files))) {
  print(tn$in.files[!file.exists(tn$in.files)])
  stop('Missing files!')
  
}


## Read regions bed
bed = read.csv(opt$bed, h=T, stringsAsFactors=F, sep='\t', check.names=F)
bed = makeGRangesFromDataFrame(bed, seqnames.field='#chr', keep.extra.columns=T)


## Read gene status
gene = read.csv(opt$gene_status, h=T, stringsAsFactors=F, sep='\t')
gene = gene[gene$sample %in% tn$pair_name, ]
gene = gene[match(gene$sample, tn$pair_name), ]

gene$RB1_status = gsub('HIGH','High Impact', gene$RB1_status)
gene$RB1_status = gsub('MODERATE','Moderate Impact', gene$RB1_status)
gene$RB1_status = gsub('NONE','None', gene$RB1_status)


## Read amp status
amp = read.csv(opt$amp_status, h=T, stringsAsFactors=F, sep='\t')
amp$jabba_class[amp$jabba_class == ''] = 'None'
amp$jabba_class = gsub(',.*','',amp$jabba_class)
amp$ecdna = ifelse(amp$ecdna, 'ecDNA+', 'ecDNA-')

cols = c('aa_class', 'jabba_class', 'ecdna')
mta[, cols] = amp[match(mta$sample, amp$sample), cols]
mta$jabba_class = events.long[mta$jabba_class]

mta$jabba_class[is.na(mta$jabba_class)] = 'None'
mta$aa_class[is.na(mta$aa_class)] = 'None'
mta$ecdna[is.na(mta$ecdna)] = 'ecDNA-'

mta$ccnd1_amp = ifelse(mta$sample %in% amp$sample, 'CCND1 Amplification', ' ')


## Read CNV files 
cnv.all = mapply(FUN=read.cnv, f=tn$in.files, id=mta$sample_id)

## Get CN for each sample in 1MB bins
tiles = gr.tile(bed, opt$binsize)
cytoband = bed$cytoband[tiles$query.id]

## Enforce ordering
cytoband = factor(cytoband, levels=c('9p21.3', '11q13.3', 'CDK4', 'CDK6', 'RB1'))
tiles = granges(tiles)

for (i in 1:length(cnv.all)) {
  tiles = tiles %$% cnv.all[[i]]
}

seqlevelsStyle(tiles) = 'NCBI'
chr = as.factor(seqnames(tiles))
res = t(as.matrix(mcols(tiles)))


## Ploidy-normalize
res.norm = res
res.norm[, !chr %in% c('X','Y')] = res.norm[, !chr %in% c('X','Y')] / tn$ploidy
res.norm[, chr %in% c('X','Y')] = res.norm[, chr %in% c('X','Y')] / tn$ploidy_xy

res.norm[res.norm == 0] = 1E-9
res.norm = log2(res.norm)



#############################
## Some summary statistics ##
############################# 

# res.norm.amp.summary = 2^res.norm[rownames(res.norm) %in% mta$sample_id[mta$ccnd1_amp == 'CCND1 Amplification'], ]
res.norm.amp.summary = 2^res.norm[rownames(res.norm) %in% mta$sample_id[mta$ecdna == 'ecDNA+'], ]

mta$ecdna

modal.amp.cn  = c()
for (i in 1:nrow(res.norm.amp.summary)) {
  modal.amp.cn = rbind(modal.amp.cn, tapply(res.norm.amp.summary[i, ], cytoband, numeric.mode, na.rm=T))
}

message('Mean fold change relative to ploidy:')
apply(modal.amp.cn, 2, mean)



########## 
## Plot ##
##########

svg(opt$out_file, width=17, height=11)
  

## Color scales
col.raw = colorRamp2(c(0, 2, 6), c("blue", "white", "red"))
col.norm = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

## Metadata
left.ano = ComplexHeatmap::rowAnnotation(`CCND1 ecDNA`=mta$ecdna,
                                        `Patient`=mta$patient,
                                          col=list(`CCND1 ecDNA`=ecdna.colors,
                                                   `Patient`=build.colormap(mta$patient, mta$patient_color)),
                                        gap=unit(c(0, 0, 0), "mm"),
                                        gp = gpar(col="black"))

right.ano = ComplexHeatmap::rowAnnotation(`RB1 Mutation`=gene$RB1_status,
                                          col=list(`RB1 Mutation`=MUT_COL),
                                          gp = gpar(col="black"))


mta$ccnd1_amp = factor(mta$ccnd1_amp, levels=c('CCND1 Amplification', ' '))


## Raw CN
ht = Heatmap(res, 
        col=col.raw,
        column_split=cytoband, 
        row_split=mta$ccnd1_amp,
        left_annotation=left.ano,
        right_annotation=right.ano,
        heatmap_legend_param = list(title = "CN"),
        row_names_gp = gpar(fontsize=8),
        border=T,
        cluster_columns=F, 
        cluster_rows=T, 
        show_row_names=T)


## log2(CN/ploidy)
ht = Heatmap(res.norm,
        col=col.norm,
        column_split=cytoband,
        row_split=mta$ccnd1_amp,
        left_annotation=left.ano,
        right_annotation=right.ano,
        heatmap_legend_param = list(title = "log2(CN/ploidy)"),
        row_names_gp = gpar(fontsize=8),
        border=T,
        cluster_columns=F,
        cluster_rows=T,
        show_row_names=T)
        
draw(ht, column_title=paste0('Ploidy-normalized Jabba Copy Number Heatmap, ',opt$binsize/1E6,'MB Bins'))

dev.off()
message(opt$out_file)



##############################
## Output source data table ##
##############################


res.norm = cbind(c(NA, mta$ecdna), 
                 c(NA, mta$Patient),
                 c(NA, gene$RB1_status), 
                 rbind(as.character(cytoband), res.norm))

out.file.txt = gsub('\\.svg$', '.txt', opt$out_file)

write.table(res.norm, out.file.txt, row.names=T, col.names=F, quote=F, sep='\t')
message(out.file.txt)