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
## Plot barplot summarizing the overlap between SV types and kataegis events
libs = c('optparse', 'stringr', 'ggplot2', 'reshape2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Read SigProfilerClusters annotated output
read.sigprofiler = function(f) {

  if (!file.exists(f)) {
    return(NULL)
  }

  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')

  return(x)

}



## Read jabba footprints
read.jabba.fp = function(f, sample_id) {

  if (!file.exists(f)) {
    return(NULL)
  }
  
  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  x$sample = sample_id

  return(x)

}



## Get arguments
option_list = list(
  make_option(c("-i", "--in_dir_jba"),   type='character', help="Jabba output directory"),
  make_option(c("-I", "--in_dir_ktg"),   type='character', help="Directory holding SigProfilerClusters kataegis VCFs"),
  make_option(c("-f", "--ktg_flag"),     type='character', help="Flag to set SigProfilerClusters run type"),
  make_option(c("-t", "--tn_file"),      type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-o", "--out_file"),     type='character', help="SVG output"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor-normal pairs and get input files
tn = read.csv(opt$tn, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)

tn$in_files_ktg = paste0(opt$in_dir_ktg,'/',tn$pair_name,'-mutationtimer-snv.annotated.jabba.txt')
tn$in_files_jba = paste0(opt$in_dir_jba,'/',tn$pair_name,'/jabba.events.footprints.kataegis.',opt$ktg_flag,'.txt')



## Read input 
ktg = do.call(rbind, lapply(tn$in_files_ktg, read.sigprofiler))
jba = do.call(rbind, mapply(read.jabba.fp, f=tn$in_files_jba, sample_id=tn$tumor, SIMPLIFY=F))

ec = jba[jba$final_type == 'ecDNA', ]
ec = as.data.frame(table(as.numeric(table(ec$sample))))
colnames(ec) = c('Number of ecDNA events', 'Number of samples')


## Exclude any events that only include breakpoints in their event footprints,
## as it's not possible for these to contain a kataegis event
jba = jba[!jba$type %in% c('chromoplexy', 'qrp', 'tra'), ]


## Update event names
jba$final_type = gsub('invdup','Inverted Duplication',jba$final_type)
jba$final_type = gsub('inv','Inversion',jba$final_type)
jba$final_type = gsub('del','Deletion',jba$final_type)
jba$final_type = gsub('dup','Duplication',jba$final_type)
jba$final_type = gsub('rigma','Rigma',jba$final_type)
jba$final_type = gsub('pyrgo','Pyrgo',jba$final_type)
jba$final_type = gsub('chromothripsis','Chromothripsis',jba$final_type)
jba$final_type = gsub('tic','TIC',jba$final_type)
jba$final_type = gsub('cpxdm','Complex DM',jba$final_type)
jba$final_type = gsub('bfb','BFB (non-ecDNA)',jba$final_type)
jba$final_type = gsub('dm','DM (non-ecDNA)',jba$final_type)
jba$final_type = gsub('tyfonas','Tyfonas (non-ecDNA)',jba$final_type)


## Get stats, order events by proportion of ecDNA
jba.prop = as.data.frame.matrix(table(jba$final_type, jba$has_kataegis))
jba.prop$prop = jba.prop$`TRUE` / rowSums(jba.prop)
event.order = rownames(jba.prop)[order(jba.prop$prop, decreasing=T)]

jba.stats = melt(table(jba$final_type, jba$has_kataegis))
colnames(jba.stats) = c('Event', 'Has kataegis', 'count')
jba.stats$Event = factor(jba.stats$Event, levels=event.order)


jba.stats$prop_label = paste0(round(jba.prop$prop[match(as.character(jba.stats$Event), rownames(jba.prop))]*100, 1),'%')
jba.stats$prop_label[jba.stats$`Has kataegis`] = ''


## Get number of kataegis events associated with ecDNA
ktg = ktg[!duplicated(ktg$kataegis_id), ] 
ktg.ecnda.prop = round(sum(ktg$ecdna) / nrow(ktg), 3)
message('Kataegis events in ecDNA: ', sum(ktg$ecdna), '/',nrow(ktg), '(', round(ktg.ecnda.prop, 3)*100,'%)')


## Get number of APOBEC kataegis events associated with ecDNA
ktg.ecdna = ktg[!duplicated(ktg$kataegis_id) & ktg$ecdna, ] 
ktg.ecnda.apobec.prop = round(sum(ktg.ecdna$apobec_associated) / nrow(ktg.ecdna), 3)
message('APOBEC Kataegis events in ecDNA: ', sum(ktg.ecdna$apobec_associated), '/',nrow(ktg.ecdna), '(', round(ktg.ecnda.apobec.prop, 3)*100,'%)')



##########
## Plot ##
##########

svg(opt$out_file, width=10.5, height=6)

col = c(`TRUE`='#E5809B', `FALSE`='#FFFFFF')

## Filter for events that happen at least 5 number of times for visual clarity
## Still include these in the source/supplemental data
event.counts = tapply(jba.stats$count, jba.stats$Event, sum)
events.keep = names(event.counts)[event.counts > 5]
jba.stats = jba.stats[jba.stats$Event %in% events.keep, ]


levels(jba.stats$Event) = paste0(levels(jba.stats$Event), ' [n=',event.counts[levels(jba.stats$Event)],']')


## Proportion bar plot
ggplot(jba.stats, aes(x=Event, y=count, fill=`Has kataegis`)) +
  geom_bar(stat='identity', position='fill') +
  geom_text(aes(label=prop_label), position=position_fill(vjust=0.1), size=5) +
  scale_fill_manual(values=col) +
  scale_y_continuous(expand=c(0,0)) +
  ylab('Proportion of events with Kataegis') +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        text=element_text(size=17),
        legend.position='none')


dev.off()
message(opt$out_file)



###########################
## Write out source data ##
###########################

out.file.txt = gsub('\\.svg$', '.txt', opt$out_file)

write.table(jba.prop, out.file.txt, row.names=T, col.names=T, quote=F, sep='\t')
message(out.file.txt)
