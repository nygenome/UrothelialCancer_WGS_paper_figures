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
## Plot oncoprint of jabba events that overlap FishHook hits

libs = c('optparse', 'ComplexHeatmap', 'circlize', 'GenomicRanges', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Pulled directly from https://github.com/mskilab/GxG/blob/2e5d3f88075520e90762b409f95eff2694fc422b/R/utils.R
alpha = function(col, alpha) {
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}



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



parseFootprints = function(f, gr) {
  
  events.long = c('Translocation', 'Inversion', 'Inverted Duplication', 'Deletion', 'Duplication', 'Chromoplexy', 
                  'Rigma', 'Pyrgo', 'Chromothripsis', 'TIC', 'QRP', 'BFB','Complex DM', 'Double Minute', 'Tyfonas')
  events.short = c('tra','inv', 'invdup', 'del', 'dup', 'chromoplexy', 'rigma', 'pyrgo', 'chromothripsis', 'tic', 
                   'qrp', 'bfb', 'cpxdm', 'dm', 'tyfonas')
  
  ## Read file 
  if (file.info(f)$size == 0) {

     ## Set up result matrix
    res = as.data.frame(matrix(0, nrow=length(unique(gr$name)), ncol=length(events.short)))
    colnames(res) = events.short
    rownames(res) = unique(gr$name)

    colnames(res) = events.long
    return(res)

  }

  ## Read footprints
  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  fp = mapply(fp2gr, x=x$footprint, id=x$type)
  fp = Reduce(c, fp)
  
  fp$id = factor(fp$id, levels=events.short)
  levels(fp$id) = events.long

  ## Overlap with bed of interest
  fp = gr.val(fp, gr, 'name', sep=',')

  if (any(grepl(',', fp$name, fixed=T))) {
    message('Warning: some footprints overlap multiple BED intervals')
    fp$name = sapply(fp$name, function(z) unlist(strsplit(z, ',', fixed=T))[1])
  }

  ## Factor so we always get the same-sized matrix
  fp$name = factor(fp$name, levels=unique(gr$name))

  ## Handle multiple overlaps for a single interval/name
  ## Just take largest overlap
  fp$overlap = fp %O% gr
  fp = fp[order(fp$overlap, decreasing=T)]
  fp = fp[!duplicated(fp$name)]
  
  res = as.data.frame.matrix(table(fp$name, fp$id))

  return(res)

}



## Need to use a closure to dynamically generate coloring function
colorClosure = function(x, y, w, h, col, alpha=1) {
  
  use.col = alpha(col, alpha)
  
  function(x, y, w, h) {
    grid.rect(x, y, 0.95*w, 0.95*h, gp=gpar(fill=use.col, col=NA))
  }
  
}



event.colors  = c(`Translocation`='#B2DF8A', `Inversion`='#05FF03', `Inverted Duplication`='#FF1293',`Deletion`='#A6CEE3', `Duplication`='#FB9A99', 
                  `Chromoplexy`='#33A02C', `Rigma`='#1F78B4', `Pyrgo`='#E31A1C', `Chromothripsis`='darkgreen', `TIC`='orchid4', 
                  `QRP`='gold2', `BFB`='firebrick', `Tyfonas`='gray25', `Double Minute`='darkorange', `Complex DM`='darkorange3')



## Get arguments
option_list = list(
  make_option(c("-j", "--jba_dir"),  type='character', help="Directory holding JaBbA results"),
  make_option(c("-t", "--tn"),       type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-b", "--bed"),      type='character', help="BED file to overlap with event footprints"),
  make_option(c("-m", "--metadata"), type='character', help="Tab delimited headered metadata file"),
  make_option(c("-o", "--out_file"), type='character', help="Figure output SVG"))
opt = parse_args(OptionParser(option_list=option_list))



## Parse gene info
bed = read.csv(opt$bed, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','start','end','name'))
bed = makeGRangesFromDataFrame(bed, keep.extra.columns=T)

message('Using intervals:')
print(as.data.frame(bed))


## Read tumor normal pairing file
tn = read.csv(opt$tn, h=F, stringsAsFactors=F, sep='\t')
colnames(tn) = c('tumor','normal','gender','patient')
tn$names = gsub('-Case-WGS','',tn$tumor)
tn$pair_name = paste0(tn$tumor,'--',tn$normal)


## Read in the genes associated with event footprints
jba.files = paste0(opt$jba_dir,'/',tn$pair_name,'/jabba.events.genes.tab')
res = lapply(jba.files, parseFootprints, gr=bed)


## Read in metadata, align with tumor-normal pairs
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
mta$sample = paste0(mta$tumor,'--',mta$normal)
mta = mta[match(tn$pair_name, mta$sample), ]
names(res) = mta$sample_id


## Reformat so we have one matrix per event, rather than one per sample
## Also, build color functions
res.oncoprint = list()
color.fxn = list()

for (i in names(event.colors)) {
  
  res.oncoprint[[i]] = do.call(cbind, lapply(res, function(x) x[, i]))
  rownames(res.oncoprint[[i]]) = unique(bed$name)
  color.fxn[[i]] = colorClosure(x, y, w, h, col=event.colors[i]) 
  
}


## Construct last 'dummy' alteration category to work around thin double minute line 
res.oncoprint$DUMMY = res.oncoprint[[1]]
for (i in 1:ncol(res.oncoprint$DUMMY)) {
  res.oncoprint$DUMMY[,i] = 0
}

color.fxn[['DUMMY']] = colorClosure(x, y, w, h, col="#F7F7F7") 
color.fxn[['background']] = colorClosure(x, y, w, h, col="#E8E8E8") 


## Metadata colors
mta$Patient = mta$patient
mta$patient[mta$patient_color == '#FFFFFF'] = 'Singleton'

patient.colors = structure(mta$patient_color, names=mta$patient)
patient.colors = na.omit(patient.colors[!duplicated(mta$patient)])

chemo.colors = structure(mta$platinum_chemotherapy_color, names=mta$platinum_chemotherapy)
chemo.colors = na.omit(chemo.colors[!duplicated(mta$platinum_chemotherapy)])

localization.colors = structure(mta$localization_color, names=mta$localization)
localization.colors = na.omit(localization.colors[!duplicated(mta$localization)])



##########
## Plot ##
##########

svg(opt$out_file, width=14, height=11)
ano = HeatmapAnnotation(`Patient`=mta$patient,
                        `Chemotherapy`=mta$platinum_chemotherapy,
                        `Localization`=mta$localization,
                         col=list(`Patient`=patient.colors,
                                  `Chemotherapy`=chemo.colors,
                                  `Localization`=localization.colors),
                          gp = gpar(col="black"))

plt = oncoPrint(res.oncoprint, 
                alter_fun=color.fxn, 
                col=c(event.colors, `DUMMY`='white'), 
                column_title='Jabba SV events overlapping with FishHook Hits\nCollapsed on cytoband and nearest CGC gene (within 500KB)',
                row_title='FishHook Hits',
                show_column_names=T,
                show_pct=F,
                border=T,
                bottom_annotation=ano,
                top_annotation=NULL,
                # column_split=tn$patient, 
                remove_empty_rows=T, 
                column_names_gp=gpar(fontsize=8))  
  
print(plt)

dev.off()
message(opt$out_file)



###################################
## Format for source data output ##
###################################

for (i in 1:length(res.oncoprint)) {

  res.oncoprint[[i]] = as.matrix(res.oncoprint[[i]])
  res.oncoprint[[i]][res.oncoprint[[i]] == 1] = names(res.oncoprint)[i]
  res.oncoprint[[i]][res.oncoprint[[i]] == 0] = ''

}


## Paste everything together 
res = Reduce(f=function(x, y) matrix(paste(x,y,sep=','), nrow=nrow(x), dimnames=dimnames(x)), res.oncoprint)
res = gsub(',+$|^,+', '', res)
res = gsub(',+', ',', res)
res = t(cbind(t(res),mta[, c('Patient', 'platinum_chemotherapy', 'localization')]))


## Write result
out.file.txt = gsub('\\.svg$', '.txt', opt$out_file)
write.table(res, out.file.txt, row.names=T, col.names=T, quote=F, sep='\t')
message(out.file.txt)
