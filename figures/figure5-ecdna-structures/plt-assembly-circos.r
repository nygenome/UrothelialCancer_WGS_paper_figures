#!/nfs/sw/R/R-4.0.0/bin/Rscript
#!/bin/bash
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
## Generate circos plot of circular contigs 
libs = c('optparse', 'GenomicRanges','gTrack', 'bamUtils', 'gUtils', 'circlize', 'reshape2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Convert a jabba footprint record into GRanges
fp2gr = function(x, id, type) {
  
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
  res$type = type
  
  return(res) 
  
}



## Select intervals in subject overlapping query, and clipping them 
## NOTE: only works when query is a single interval
gr.select = function(subject, query, trim=T) {

  if (length(query) != 1) {
    stop('gr.select only works with a single query interval')
  }

  subject = subject[subject %^% query]

  min.start = min(start(query))
  if (any(start(subject) < min.start)) {
    start(subject)[start(subject) < min.start] = min.start
  }

  max.end = max(end(query))
  if (any(end(subject) > max.end)) {
    end(subject)[end(subject) > max.end] = max.end
  }

  return(subject)

}




## Convert a footrpint table into a single GRanges
fp.table.to.gr = function(x) {

  x$id = make.unique(x$type, sep='_')

  res = mapply(FUN=fp2gr, 
               x=x$footprint, 
               id=x$id, 
               type=x$type, 
               SIMPLIFY=F)

  res = Reduce(c, res)

  return(res)

}



## Read gencode GTF and format 
read.gencode = function(f, genes=NULL) { 

    ## Read parsed GTF
    x = readRDS(f)

    ## Optionally select genes of interest
    if (!is.null(genes)) {
      x = x[genes]
    }

    ## For each gene 
    for (i in 1:length(x)) {

      ## Filter out gene body
      x[[i]] = x[[i]][!is.na(x[[i]]$transcript_name)]

      ## Get widths of each transcript 
      w = sapply(split(x[[i]], x[[i]]$transcript_name), function(y) sum(width(y)))

      ## Choose longest transcript 
      x[[i]] = x[[i]][x[[i]]$transcript_name %in% names(w)[which.max(w)]]

    }

    ## Combine into single granges
    x = Reduce(f=c, x)
    return(x)

}



read.regulatory.elements = function(f) {

  x = read.csv(f, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr', 'start','end','score','color'))
  x = makeGRangesFromDataFrame(x, keep.extra=T)

  ## Parse colors to hex
  x$color = gsub('\\(|\\)', '', x$color)
  x$color = sapply(strsplit(x$color, ','), function(y) rgb(y[1], y[2], y[3], maxColorValue=255))

  return(x)

}



## Initialize circos track from granges
init.circos = function(gr) {

  ctg.info = as.data.frame(gr)[, c('seqnames', 'start', 'end')]
  ctg.info$seqnames = as.character(ctg.info$seqnames)

  circos.genomicInitialize(ctg.info, tickLabelsStartFromZero=F, major.by=1E5)

}



## Plot a single read 
plot.read = function(read) {

  read = as.data.frame(read)
  
  circos.genomicTrack(read, ylim = c(-0.1, 0.1), stack=F,
    panel.fun = function(region, value, ...) {
            current_tx_start = min(region$start)
            current_tx_end = max(region$end)

            circos.lines(c(current_tx_start, current_tx_end), 
                         c(0, 0),  
                         lwd=1.1, 
                         col="#000000")

            circos.genomicRect(region, 
                               ytop=0.05, 
                               ybottom =-0.05, 
                               col="#7F0E7F", 
                               border=NA)

  }, bg.border = NA, track.height = 0.1)

}



## Plot transcripts
plot.transcripts = function(tx) {

  tx = as.data.frame(tx)
  
  circos.genomicTrack(tx, ylim = c(-0.1, 0.1), stack=F,
    panel.fun = function(region, value, ...) {
      all_tx = unique(value$gene_name)
      for (i in all_tx) {
        
        current_tx_start = min(region$start[value$gene_name == i])
        current_tx_end = max(region$end[value$gene_name == i])
        arrow.dir = ifelse(unique(value$strand[value$gene_name == i]) == '+', 'end', 'start')

        circos.genomicRect(region[value$gene_name == i, ], 
                           ytop=0.05, 
                           ybottom =-0.05, 
                           col="grey50", 
                           border=NA)
                           
         circos.arrow(x1=current_tx_start, 
                      x2=current_tx_end, 
                      y=0,
                      width=0.01,
                      arrow.head.length=mm_x(1),
                      arrow.position=arrow.dir,
                      col="#000000")

        midpoint = current_tx_start + ((current_tx_end - current_tx_start) / 2)
        circos.text(x=midpoint, 
                    y=0.1, 
                    cex=0.6,
                    labels=i, 
                    facing='bending.inside', 
                    niceFacing=T)
      }
  }, bg.border = NA, track.height = 0.2)

}



## Plot generic range data (e.g. promoters)
plot.ranges = function(x, use.colors=F) {

  h = 0.05

  if (!use.colors) {
    x$color = '#000000'
  }

  x = as.data.frame(x)

  set_track_gap(mm_h(0.01))
  circos.genomicTrack(x, ylim = c(0, h), stack=F,
    panel.fun = function(region, value, ...) {

            current_tx_start = min(region$start)
            current_tx_end = max(region$end)

            circos.genomicRect(region, 
                               ytop=h, 
                               ybottom=0, 
                               col=value$color, 
                               border=NA)
  }, bg.border = NA, track.height = h)

}



## Get arguments
option_list = list(
  make_option(c("-r", "--regions"),        type='character', help="JaBbA event footprints"),
  make_option(c("-c", "--ont_cov"),        type='character', help="Path to ONT mosdepth coverage"),
  make_option(c("-f", "--flye_dir"),       type='character', help="Top-level Flye output directory"),
  make_option(c("-g", "--gtf"),            type='character', help="gTrack-parsed GTF with gene models"),
  make_option(c("-p", "--promoters"),      type='character', help="gTrack-parsed GTF with gene models"),
  make_option(c("-e", "--enhancers"),      type='character', help="gTrack-parsed GTF with gene models"),
  make_option(c("-x", "--superenhancers"), type='character', help="gTrack-parsed GTF with gene models"),
  make_option(c("-s", "--sample"),         type='character', help="Sample ID"),
  make_option(c("-o", "--out_file"),       type='character', help="Figure output file (PDF)"))
opt = parse_args(OptionParser(option_list=option_list))


AMPS = c('dm', 'bfb', 'tyfonas')
COV_FLAG = '.pcr_free'
GENES = c('CCND1', 'MYC', 'LTO1', 'MYEOV', 'FGF3', 'FGF4', 'FGF19')

## Read footprints 
jba.fp = read.csv(opt$regions, h=T, stringsAsFactors=F, sep='\t')
jba.fp = jba.fp[jba.fp$type %in% AMPS, ]
jba.fp = fp.table.to.gr(jba.fp)

## Read ONT coverage 
cov = rtracklayer::import(opt$ont_cov)
cov$name = as.numeric(cov$name)
cov$name = cov$name / median(cov$name)

## Read gencode GTFs 
gtf = read.gencode(opt$gtf, genes=GENES)

## Read promoters/enhancers/superenhancers
prm = read.regulatory.elements(opt$promoters)
enh = read.regulatory.elements(opt$enhancers)
sup = rtracklayer::import(opt$superenhancers)



##########
## Plot ##
##########

pdf(opt$out_file, width=6, height=6)

for (i in unique(jba.fp$id)) {

  ## Find minimap2-aligned assembly BAM
  message('Flag ', COV_FLAG)
  bam = paste0(opt$flye_dir,'/',opt$sample,'.',i,COV_FLAG,'/minimap2/',opt$sample,'.',i,'-assembly.bam')

    ## Only plot circular contigs 
  message('Selecting circular contigs')
  asm.info = paste0(opt$flye_dir,'/',opt$sample,'.',i,COV_FLAG,'/assembly_info.txt')
  asm.info = read.csv(asm.info, h=T, stringsAsFactors=F, sep='\t')
  circ.contig.ids = asm.info$`X.seq_name`[asm.info$circ. == 'Y']
  
  ## Initialize circos plot 
  init.circos(jba.fp[jba.fp$id == i])
  circos.par('track.height' = 0.1, 'track.margin' = c(0.01, 0.01))


  ## Add assembly to plot 
  asm = read.bam(bam=bam, 
                   intervals=jba.fp[jba.fp$id == i], 
                   ignore.indels=T)

  asm = asm[[which(sapply(asm, function(x) all(x$qname %in% circ.contig.ids)))]]
  
  plot.read(asm)


  ## Add coverage 
  sel.fp = jba.fp[jba.fp$id == i] - 1E3
  circos.genomicTrack(as.data.frame(cov[cov %^% sel.fp]), numeric.column = 6, 
    panel.fun = function(region, value, ...) {
        circos.genomicPoints(region, value, pch=20, cex=0.3, ...)
    })


  ## Add CGC gene models 
  plot.transcripts(gtf[gtf %^% sel.fp])


  ## Add promoters 
  plot.ranges(gr.select(sup, sel.fp))
  plot.ranges(gr.select(enh, sel.fp))
  plot.ranges(gr.select(prm, sel.fp), use.colors=T)

  ## Add size of assembled amplicon
  asm.length = paste(round(asm$qwidth[asm$flag == 0] / 1E3, 2), 'KB')
  text(0, 0, paste0(asm.length, ' ecDNA'))

  circos.clear()

}

dev.off()
message(opt$out_file)
    