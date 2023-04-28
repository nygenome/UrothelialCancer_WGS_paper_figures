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
## Plot jabba track and coverage track for 

libs = c('optparse', 'GenomicRanges', 'gGnome', 'gTrack', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))



alpha = function(col, alpha) {
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}



## If intervals are really close to eachother, gtrack will shift some vertically
## Normally this is desired but here it isn't
collapse.on = function(x, f) {
  
  res = unlist(reduce(split(x, mcols(x)[, f])))
  res = res %$% x
  names(res) = c()
  
  return(res)
  
}


aa.cols  = c(`Complex non-cyclic`='springgreen3',
             `Cyclic`='palevioletred2',
             `Linear amplification`='deepskyblue3')


events.long = c('Simple (TRA/INV)', 'Deletion', 'Duplication', 'Chromoplexy', 'Rigma', 'Pyrgo', 'Chromothripsis', 'TIC', 'QRP', 'BFB', 'Double Minute', 'Complex DM', 'Tyfonas')
events.short = c('simple', 'del', 'dup', 'chromoplexy', 'rigma', 'pyrgo', 'chromothripsis', 'tic', 'qrp', 'bfb', 'dm', 'cpxdm','tyfonas')
event.colors  = c(`Simple (TRA/INV)`='#B2DF8A', `Deletion`='#A6CEE3', `Duplication`='#FB9A99',
                  `Chromoplexy`='#33A02C', `Rigma`='#1F78B4', `Pyrgo`='#E31A1C', `Chromothripsis`='darkgreen',
                  `TIC`='orchid4', `QRP`='gold2', `BFB`='firebrick', `Tyfonas`='gray25', `Double Minute`='darkorange',
                  `Complex DM`='darkorange3')


## Get arguments
option_list = list(
  make_option(c("-j", "--jabba_gg"),        type='character', help="Path to jabba.events.rds file"),
  make_option(c("-c", "--coverage"),        type='character', help="Coverage GRanges object (RDS)"),
  make_option(c("-n", "--normal_coverage"), type='character', help="Coverage GRanges object (RDS)"),
  make_option(c("-f", "--field"),           type='character', help="Field in --coverage file to use when plotting coverage", default='foreground'),
  make_option(c("-a", "--aa"),              type='character', help="Optional BED file with AmpliconArchitect results"),
  make_option(c("-b", "--bed"),             type='character', help="BED file with regions for plotting"),
  make_option(c("-p", "--padding"),         type='numeric',   help="Padding for taking subgraph associated with BED intervals"),
  make_option(c("-o", "--out_file"),        type='character', help="Figure output file (PNG)"),
  make_option(c("-x", "--out_file_fp"),     type='character', help="Optional output of subgraph footprint"))
opt = parse_args(OptionParser(option_list=option_list))



## Read BED, format intervals
## This is the first mandatory track
bed = read.csv(opt$bed, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','start','end','label'))
bed = makeGRangesFromDataFrame(bed, keep.extra.columns=T)

names(bed) = bed$label


bed.gtrack = gTrack(bed, 
                    col=alpha('black', 0.5), name=' ', height=3,  labels.suppress=F, gr.labelfield='label',
                    grl.labelfield='label',  vadj.label=-3, cex.label=3, xaxis.cex.tick=1.5,
                    xaxis.cex.label=1.5, xaxis.unit=1e6, xaxis.suffix='MB')


## Read tumor coverage
## This is the second mandatory track
td.cov = gTrack(readRDS(opt$coverage), y.field=opt$field, col=alpha('black', 0.2), name='Cov')



## Optionally add genome graph colored by event type
if (!is.null(opt$jabba_gg)) {
  ## Read genome graph and coverage
  jabba.gg = readRDS(opt$jabba_gg)
  
  ## Reset node/edge coloring
  jabba.gg$nodes$mark(col='gray')
  jabba.gg$edges[type=='ALT']$mark(col='gray50')
  
  
  ## Mark edges by event type 
  for (e in 1:length(events.short)) {
    
    jabba.gg$edges[!is.na(get(events.short[e]))]$mark(col=event.colors[events.long[e]])
  }
  
  ## Optionally mark "altered" nodes according to cutoff
  if (!is.null(opt$fga_cutoff) && !is.null(opt$ploidy) && !is.null(opt$gender)) {
    
    jabba.gg$nodes[abs(log2(cn/opt$ploidy)) >= opt$fga_cutoff & !seqnames %in% c('chrX','chrY')]$mark(col='black')
    jabba.gg$nodes[abs(log2(cn/(opt$ploidy/2))) >= opt$fga_cutoff & seqnames %in% c('chrX','chrY')]$mark(col='black')
    
  }
  
  ## Build gtrack
  ## Only show colors in legend that are in plot
  jabba.gg = jabba.gg$subgraph(bed, opt$padding)
  event.colors = event.colors[which(event.colors %in% unique(as.data.frame(jabba.gg$edges$dt)$col))]
  print(event.colors)
  jabba.gtrack = jabba.gg$gtrack(colormaps=list(`Jabba event`=event.colors))
  plt.gtrack = c(td.cov, jabba.gtrack)
  
  if (!is.null(jabba.gg$agtrack)) {
    plt.gtrack = c(jabba.gg$agtrack, plt.gtrack)
  }
  
} else {
  
  plt.gtrack = td.cov
  
}



## Optionally add normal coverage track
if (!is.null(opt$normal_coverage)) {
  
  ncov = gTrack(readRDS(opt$normal_coverage), y.field=opt$field, col=alpha('black', 0.2), name='Normal Cov')
  plt.gtrack = c(plt.gtrack, ncov)
  
}



## Optionally add AmpliconArchitect results
if (!is.null(opt$aa) && file.exists(opt$aa)) {
  
  aa = makeGRangesFromDataFrame(read.csv(opt$aa, h=T, stringsAsFactors=F, sep='\t'), keep.extra.columns=T)
  aa = collapse.on(aa[aa$aa_id != ''], 'aa_id')

  aa.gtrack = gTrack(aa, 
                     name='Amplicon\nArchitect', 
                     colormaps=list(aa_class=aa.cols), 
                     height=2)
  

  plt.gtrack = c(plt.gtrack, aa.gtrack)
  
  if (any(aa %^% bed)) {
    cat(paste(unique(aa$aa_class[aa %^% bed]), collapse=','))
  }


  ## Export source data
  aa.out = as.data.frame(aa[aa %^% jabba.gg$footprint])[, c('seqnames', 'start', 'end', 'aa_class')]
  aa.out.txt = gsub('\\.svg$', '.aa.txt', opt$out_file)

  write.table(aa.out, aa.out.txt, row.names=F, col.names=T, quote=F, sep='\t')
  message(aa.out.txt)

}


##########
## Plot ## 
##########

if (substr(opt$out_file, nchar(opt$out_file)-2, nchar(opt$out_file)) == 'png') {
  png(opt$out_file, width=2000, height=1000)  
} else {
  svg(opt$out_file, width=18, height=13)
}


region = bed + opt$padding
  
# plot(c(bed.gtrack, plt.gtrack), region, y.pad=1, cex.ylabel=0.6, cex.xlabel=1.5)
plot(c(bed.gtrack, plt.gtrack), jabba.gg$footprint)
  
x = dev.off()
message(opt$out_file)



########################
## Export source data ##
########################


node.cols = c('snode.id', 'seqnames', 'start', 'end', 'cn', 'loose.left',	'loose.right', 'col')
edge.cols = c('sedge.id', 'chr1', 'start1', 'chr2', 'start2', 'type', 'str1', 'str2', 'n1.side', 'n2.side', 'n1', 'n2','col')

nodes = as.data.frame(jabba.gg$nodes$dt)[, node.cols]
edges = as.data.frame(jabba.gg$edges$dt)[, edge.cols]

nodes.out.txt = gsub('\\.svg$', '.nodes.txt', opt$out_file)
edges.out.txt = gsub('\\.svg$', '.edges.txt', opt$out_file)

write.table(nodes, nodes.out.txt, row.names=F, col.names=T, quote=F, sep='\t')
write.table(edges, edges.out.txt, row.names=F, col.names=T, quote=F, sep='\t')

message(nodes.out.txt)
message(edges.out.txt)




## Write out footprints for downstream use 
if (!is.null(opt$out_file_fp)) {
  
  fp = cbind(as.character(seqnames(jabba.gg$footprint)), 
             start(jabba.gg$footprint),
             end(jabba.gg$footprint))
  
  write.table(fp, opt$out_file_fp, col.names=F, row.names=F, sep='\t', quote=F)
  message(opt$out_file_fp)
  
}

