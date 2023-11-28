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
## Summarize AmpliconArchitect results
libs = c('optparse', 'rtracklayer', 'GenomicRanges', 'ggplot2', 'ggbeeswarm', 'reshape2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


## Import all BED files associated with an amplicon
read.amplicons = function(f, id) {

  ## Read file(s)
  res = lapply(f, rtracklayer::import, format='bed')

  ## Filter empty files
  non.empty.files = which(sapply(res, length) > 0)
  res = res[non.empty.files]
  id = id[non.empty.files]

  ## Add ID
  for (i in 1:length(res)) {
    res[[i]]$id = id[i]
  }
  
  ## Collapse to single GRanges and return
  res = Reduce(c, res)
  return(res)

}



## Read amplicon cycles file 
## Try computing complexity score as desribed in AA/AC 
## Barrett's esophagus paper methods
read.cycles = function(f) { 

  ## Read cycles
  cmd = paste('grep "^Cycle"', f)
  x = system(cmd, intern=T)

  x = colsplit(x, pattern=';', names=c('cycle','copy_count','segments'))
  x$cycle = gsub('Cycle=', '', x$cycle, fixed=T)
  x$copy_count = as.numeric(gsub('Copy_count=', '', x$copy_count, fixed=T))
  x$segments = gsub('Segments=', '', x$segments, fixed=T)

  ## Exclude "0" nodes that mark start/end
  x$num_segments = sapply(strsplit(x$segments, ','), length) - 2 


  ## Read segments into Granges
  cmd = paste('grep "^Segment"', f)
  seg = system(cmd, intern=T)
  seg = colsplit(seg, pattern='\t', names=c('segment', 'id', 'chr','start', 'end'))[,-1]
  seg = makeGRangesFromDataFrame(seg, keep.extra.columns=T)

  ## Compute path lengths
  x$path_length = NA
  for (i in 1:nrow(x)) {

    seg.ids = unlist(strsplit(x$segments[i], ','))
    seg.ids = gsub('\\+|-', '', seg.ids)
    seg.ids = seg.ids[seg.ids != '0']
    x$path_length[i] = sum(width(seg[seg$id %in% seg.ids]))

  }

  
  ## Compute T: length weighted copy number of all decomposed paths
  T = sum(x$copy_count * x$path_length)

  ## Compute D: fraction of total CN T contributed by each path
  x$D = (x$copy_count * x$path_length) / T

  ## Sort on D in ascending order
  x = x[order(x$D, decreasing=F), ]

  ## Mark paths that are part of the residual by
  ## finding index j, index corresponding the the sum of (1:j) that is less
  ## than the 80th percentile 
  quantile.cutoff = 0.8
  j = max(which(cumsum(x$D) < quantile(x$D, quantile.cutoff)))

  ## Mark residual as indices > j + 1
  x$residual = (1:nrow(x)) %in% ((j+2):nrow(x))

  ## Compute residual (epsilon)
  eps = sum(x$D[x$residual])

  ## k is just the number of segments in the amplicon
  k = length(seg)


  ## Compute complexity 
  res = (-eps * log(eps)) - (sum(x$D[!x$residual] * log(x$D[!x$residual]))) - log(1/k)
  
  if (is.na(res)) {
    print(x)
  }

  # return(res)
  return(nrow(x))

}



## Get arguments
option_list = list(
  make_option(c("-t", "--tn_file"),                          type='character', help="Tumor-normal pairs file"),
  make_option(c("-A", "--amplicon_classification_profiles"), type='character', help="AmpliconArchitect classification table"),
  make_option(c("-f", "--aa_run_flag"),                      type='character', help="AmpliconArchitect run flag"),
  make_option(c("-x", "--pm636_z2_profiles"),                type='character', help="PM636-Z2 AmpliconArchitect classification table (4CN 10X run failed)"),
  make_option(c("-X", "--pm636_z2_aa_run_flag"),             type='character', help="PM636-Z2 AmpliconArchitect run flag (4CN 10X run failed)"),
  make_option(c("-a", "--aa_dir"),                           type='character', help="Top-level AmpliconArchitect output directory"),
  make_option(c("-c", "--classification_bed_dir"),           type='character', help="Top-level AmpliconArchitect output directory", default='classification_bed_files'),
  make_option(c("-C", "--collapse_to_amplicons"),            type='logical',   help="Use full amplicon structure, instead of decomposed cycles", default=FALSE),
  make_option(c("-o", "--out_file"),                         type='character', help="Output PDF"),
  make_option(c("-O", "--out_file_txt"),                     type='character', help="Output table"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor normal pairs 
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)

## Read AA amplicon classification summary
aa = read.csv(opt$amplicon_classification_profiles, h=T, stringsAsFactors=F, sep='\t', check.names=F) 
aa$pair_name = gsub(paste0('_',opt$aa_run_flag), '', aa$sample_name, fixed=T)
aa = aa[aa$pair_name %in% tn$pair_name, ]
# aa$classification_dir_name = 'classification_bed_files'
aa$classification_dir_name = opt$classification_bed_dir

## Update PM636-Z2 in table (PM636-Z2-1-Case-WGS--PM636-EBC1-Ctrl-WGS)
aa.pm636 = read.csv(opt$pm636_z2_profiles, h=T, stringsAsFactors=F, sep='\t', check.names=F) 
aa.pm636$pair_name = gsub(paste0('_',opt$pm636_z2_aa_run_flag), '', aa.pm636$sample_name, fixed=T)
aa.pm636 = aa.pm636[aa.pm636$pair_name == 'PM636-Z2-1-Case-WGS--PM636-EBC1-Ctrl-WGS', ]
aa.pm636$classification_dir_name = 'classification_bed_files_5X_downsample10'

message('WARNING: Swapping in previous successfuly PM636-Z2 run!')
aa = aa[aa$pair_name !='PM636-Z2-1-Case-WGS--PM636-EBC1-Ctrl-WGS', ]
aa = rbind(aa, aa.pm636)

aa = aa[aa$pair_name %in% tn$pair_name, ]



###################
## Basic summary ##
###################

aa.event.counts = as.data.frame(table(aa$amplicon_decomposition_class))
colnames(aa.event.counts) = c('event_type', 'count')
aa.event.counts

num.samples.with.amplicon = length(unique(aa$pair_name[aa$amplicon_decomposition_class != 'No amp/Invalid']))
pct.samples.with.amplicon = round(num.samples.with.amplicon / nrow(tn), 4) * 100
message('Samples with an amplicon detected: ',num.samples.with.amplicon,'/',nrow(tn), ' (', pct.samples.with.amplicon,'%)')

num.samples.with.cyclic = length(unique(aa$pair_name[aa$amplicon_decomposition_class == 'Cyclic']))
pct.samples.with.cyclic = round(num.samples.with.cyclic / nrow(tn), 4) * 100
message('Samples with a CYCLIC amplicon detected: ',num.samples.with.cyclic,'/',nrow(tn), ' (', pct.samples.with.cyclic,'%)')

num.samples.with.cpx.noncyclic = length(unique(aa$pair_name[aa$amplicon_decomposition_class == 'Complex non-cyclic']))
pct.samples.with.cpx.noncyclic = round(num.samples.with.cpx.noncyclic / nrow(tn), 4) * 100
message('Samples with a COMPLEX NON-CYCLIC amplicon detected: ',num.samples.with.cpx.noncyclic,'/',nrow(tn), ' (', pct.samples.with.cpx.noncyclic,'%)')

num.samples.with.linear.amp = length(unique(aa$pair_name[aa$amplicon_decomposition_class == 'Linear amplification']))
pct.samples.with.linear.amp = round(num.samples.with.linear.amp / nrow(tn), 4) * 100
message('Samples with a LINEAR AMPLIFICATION amplicon detected: ',num.samples.with.linear.amp,'/',nrow(tn), ' (', pct.samples.with.linear.amp,'%)')

message('\n\n')


## Filter for amplicon calls
aa = aa[aa$amplicon_decomposition_class != 'No amp/Invalid', ]



## Summarize per-sample counts
message('Amplicon per-sample count summary')
summary(as.numeric(table(aa$pair_name[aa$amplicon_decomposition_class != 'No amp/Invalid'])))

message('CYCLIC per-sample count summary')
summary(as.numeric(table(aa$pair_name[aa$amplicon_decomposition_class == 'Cyclic'])))

message('COMPLEX NON-CYCLIC per-sample count summary')
summary(as.numeric(table(aa$pair_name[aa$amplicon_decomposition_class == 'Complex non-cyclic'])))

message('LINEAR AMPLIFICATION per-sample count summary')
summary(as.numeric(table(aa$pair_name[aa$amplicon_decomposition_class == 'Linear amplification'])))



## Confirm that all Cyclic calls are ecDNA+ 
if (!all(aa$`ecDNA+`[aa$amplicon_decomposition_class == 'Cyclic'] == 'Positive')) {
  stop('ERROR: Not all cyclic amplicons are ecDNA+!')
} else {
  message('All cyclic amplicons are ecDNA+!')
}


## Check to see if we detected any BFBs
if (any(aa$`BFB+` == 'Positive')) {
  stop('AmpliconArchitect detected BFB cycles!')
} else {
  message('No BFB cycles detected!')
}



#########################################
## Summarize amplicon characterisitics ##
#########################################

## Read in individual amplicon files 
amplicons = list()
for (i in 1:nrow(aa)) {

  amplicon.class = ifelse(aa$amplicon_decomposition_class[i] == 'Cyclic', 'ecDNA', aa$amplicon_decomposition_class[i])

  in.file.pattern = paste0(opt$aa_dir,'/',aa$classification_dir_name[i],'/',aa$sample_name[i],'_',aa$amplicon_number[i],'_',aa$amplicon_number[i],'_*bed')
  in.files = Sys.glob(in.file.pattern)

  ## Construct unique, somewhat concise ID
  id = gsub(paste0(opt$aa_dir,'/',aa$classification_dir_name[i],'/',aa$sample_name[i],'_',aa$amplicon_number[i],'_'), '', in.files, fixed=T)
  id = gsub(' ', '_', id, fixed=T)
  id = gsub('.bed', '', id, fixed=T)
  id = paste0(id,'_', aa$pair_name[i])
  amplicons[[i]] = read.amplicons(f=in.files, id=id)

  amplicons[[i]]$class = amplicon.class
  amplicons[[i]]$sample = aa$pair_name[i]
  amplicons[[i]]$amplicon_number = aa$amplicon_number[i]

  ## Flag walks that were classified as unknown
  amplicons[[i]]$unknown_class_flag = grepl('unknown', amplicons[[i]]$id, fixed=T)

  ## Count number of cycles for this amplicon
  in.file.cycles = paste0(opt$aa_dir,'/',aa$pair_name[i],'/',aa$sample_name[i],'_',aa$amplicon_number[i],'_cycles.txt')
  amplicons[[i]]$cycles = read.cycles(in.file.cycles)
  if(!file.exists(in.file.cycles)) {
    stop('Missing cycles file: ', in.file.cycles)
  }

}



## Collapse into single GRanges and write result
amplicons = Reduce(c, amplicons)


## Optionally collapse walks into amplicons
if (opt$collapse_to_amplicons) {

  amplicons$id = paste0(amplicons$amplicon_number,'_',amplicons$sample)

} else {

  ## Otherwise, we don't need to keep the walks classified as "unknown"
  amplicons = amplicons[!amplicons$unknown_class_flag, ]

}

write.table(as.data.frame(amplicons), opt$out_file_txt, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file_txt)



## Summarize amplicon sizes by class
size.df = c()
for (i in unique(amplicons$class)) {

  sizes = tapply(width(amplicons[amplicons$class == i]), amplicons$id[amplicons$class == i], sum) / 1E6
  round(sizes, 2)

  message(i, ' (MB)')
  print(summary(sizes))

  size.df.i = data.frame(class=i, size_mb=sizes)
  size.df = rbind(size.df, size.df.i)

}


size.df$class = factor(size.df$class, names(sort(tapply(size.df$size_mb, size.df$class, median))))


################
## Plot sizes ##
################

pdf(opt$out_file)

aa.class.cols = c(`ecDNA`='#E07694', `Linear amplification`='#2C8EC1', `Complex non-cyclic`='#45C166')

ggplot(size.df, aes(x=class, y=size_mb, fill=class)) +
  geom_violin() +
  geom_quasirandom() +
  scale_fill_manual(values=aa.class.cols) +
  scale_y_log10() + 
  xlab('AmpliconArchitect classification') +
  ylab('Event size (MB)') +
  theme_bw() +
  theme(text=element_text(size=15),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position='none')

dev.off()
message(opt$out_file)
