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
## Summarize kyklonas events (and SBS31/35)
libs = c('optparse', 'data.table')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)

LOW_VAF_CUTOFF = 0.333
HIGH_VAF_CUTOFF = 0.667

## Read SigProfilerClusters annotated output
read.sigprofiler = function(f) {

  if (!file.exists(f)) {
    return(NULL)
  }

  x = fread(f, header=T, stringsAsFactors=F, sep='\t', data.table=F)

  return(x)

}



## Get arguments
option_list = list(
  make_option(c("-I", "--in_dir_ktg"),  type='character', help="Directory holding SigProfilerClusters VCFs"),
  make_option(c("-n", "--in_dir_nc"),   type='character', help="Directory holding SigProfilerClusters VCFs (unclustered)"),
  make_option(c("-t", "--tn_file"),     type='character', help="Tab delimited tumor-normal pairing file"))
opt = parse_args(OptionParser(option_list=option_list))



## Read TN pairs and get input files
tn = read.csv(opt$tn, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)

tn$in_files_ktg = paste0(opt$in_dir_ktg,'/',tn$pair_name,'-mutationtimer-snv.annotated.jabba.txt')
tn$in_files_nc = paste0(opt$in_dir_nc,'/',tn$pair_name,'-mutationtimer-snv.annotated.jabba.txt')


## Read input 
ktg.all = Reduce(rbind, lapply(tn$in_files_ktg, read.sigprofiler))
nc.all = Reduce(rbind, lapply(tn$in_files_nc, read.sigprofiler))
nc.all = nc.all[nc.all$tricontext != 'DBS', ]



#######################################
## Non-kyklonas breakpoint distances ##
#######################################

ktg.jabba = ktg.all[!ktg.all$jabba_event == 'none' & !ktg.all$ecdna, ]
dist.per.event = tapply(ktg.jabba$distance_to_nearest_junction, ktg.jabba$kataegis_id, mean)
message('\nMean jabba-associated kataegis distance to nearest breakpoint: ', mean(dist.per.event)/1000, ' KB')
message('Median jabba-associated kataegis distance to nearest breakpoint: ', median(dist.per.event)/1000, ' KB')


ktg.none = ktg.all[ktg.all$jabba_event == 'none' & !ktg.all$ecdna & !is.na(ktg.all$distance_to_nearest_junction), ]
dist.per.event = tapply(ktg.none$distance_to_nearest_junction, ktg.none$kataegis_id, mean)
message('\nMean NON-SV-associated kataegis distance to nearest breakpoint: ', mean(dist.per.event)/1E6, ' MB (events on non-rearranged chromosomes excluded)')
message('Median NON-SV-associated kataegis distance to nearest breakpoint: ', median(dist.per.event)/1E6, ' MB (events on non-rearranged chromosomes excluded)')



##############
## Kyklonas ##
##############

ktg = ktg.all[ktg.all$ecdna, ]

message('\nKyklonas mutations: ', nrow(ktg))
message('Kyklonas events: ', length(unique(ktg$kataegis_id)))


## Compute mean VAF per event 
mean.ktg.vaf = tapply(ktg$vaf, ktg$kataegis_id, mean, na.rm=T)
vaf.gt.80 = sum(mean.ktg.vaf > HIGH_VAF_CUTOFF)
vaf.gt.80.pct = round(vaf.gt.80 / length(mean.ktg.vaf), 3) * 100

vaf.lt.50 = sum(mean.ktg.vaf < LOW_VAF_CUTOFF)
vaf.lt.50.pct = round(vaf.lt.50 / length(mean.ktg.vaf), 3) * 100


message('\nKyklonas with VAF > ',HIGH_VAF_CUTOFF,': ', vaf.gt.80, ' (', vaf.gt.80.pct,'%)')
message('Kyklonas with VAF < ',LOW_VAF_CUTOFF,': ', vaf.lt.50, ' (', vaf.lt.50.pct,'%)')



#########################################
## Kyklonas                            ##
## Only patients highlighted in figure ##
#########################################

ktg = ktg.all[ktg.all$ecdna & ktg.all$sample %in% c('PM264-Z1-1-Case-WGS', 'PM911-Z1-1-Case-WGS'), ]

message('\nPM264-Z1/PM911-Z1 Kyklonas mutations: ', nrow(ktg))
message('PM264-Z1/PM911-Z1 Kyklonas events: ', length(unique(ktg$kataegis_id)))


## Compute mean VAF per event 
mean.ktg.vaf = tapply(ktg$vaf, ktg$kataegis_id, mean, na.rm=T)
vaf.gt.80 = sum(mean.ktg.vaf > HIGH_VAF_CUTOFF)
vaf.gt.80.pct = round(vaf.gt.80 / length(mean.ktg.vaf), 3) * 100

vaf.lt.50 = sum(mean.ktg.vaf < LOW_VAF_CUTOFF)
vaf.lt.50.pct = round(vaf.lt.50 / length(mean.ktg.vaf), 3) * 100


message('\nPM264-Z1/PM911-Z1 Kyklonas with VAF > ',HIGH_VAF_CUTOFF,': ', vaf.gt.80, ' (', vaf.gt.80.pct,'%)')
message('PM264-Z1/PM911-Z1 Kyklonas with VAF < ',LOW_VAF_CUTOFF,': ', vaf.lt.50, ' (', vaf.lt.50.pct,'%)')




## Distance to breakpoint 
ktg=ktg.all
dist.per.event = tapply(ktg$distance_to_nearest_junction, ktg$kataegis_id, mean)

message('\nMean kyklonas distance to nearest breakpoint: ', mean(dist.per.event)/1000, 'KB')
message('Median kyklonas distance to nearest breakpoint: ', median(dist.per.event)/1000, 'KB')
message('Total number of kyklonas events within 10KB of a breakpoint (<=10KB): ', sum(dist.per.event <= 1E4))
message('Total number of kyklonas events within 1MB of a breakpoint (<=1MB): ', sum(dist.per.event <= 1E6))



#########################
## Kataegis near ecDNA ##
#########################

ktg = ktg.all[!is.na(ktg.all$distance_to_nearest_ecdna) & !ktg.all$ecdna & ktg.all$distance_to_nearest_ecdna <= 1E4, ]

message('\nNear-ecDNA kataegis mutations: ', nrow(ktg))
message('Near-ecDNA kataegis events: ', length(unique(ktg$kataegis_id)))


## Compute mean VAF per event 
mean.ktg.vaf = tapply(ktg$vaf, ktg$kataegis_id, mean, na.rm=T)
vaf.gt.80 = sum(mean.ktg.vaf > HIGH_VAF_CUTOFF)
vaf.gt.80.pct = round(vaf.gt.80 / length(mean.ktg.vaf), 3) * 100

vaf.lt.50 = sum(mean.ktg.vaf < LOW_VAF_CUTOFF)
vaf.lt.50.pct = round(vaf.lt.50 / length(mean.ktg.vaf), 3) * 100


message('\nNear-ecDNA kataegis with VAF > ',HIGH_VAF_CUTOFF,': ', vaf.gt.80, ' (', vaf.gt.80.pct,'%)')
message('Near-ecDNA kataegis with VAF < ',LOW_VAF_CUTOFF,': ', vaf.lt.50, ' (', vaf.lt.50.pct,'%)')



###################
## Non-clustered ##
###################

nc = nc.all[nc.all$jabba_event == 'ecDNA', ]
message('Median non-clustered ecDNA mutation distance to nearest breakpoint: ', median(nc$distance_to_nearest_junction)/1000, 'KB')

nc = nc.all[!nc.all$jabba_event %in% c('none', 'ecDNA'), ]
message('Median non-clustered jabba mutation distance to nearest breakpoint: ', median(nc$distance_to_nearest_junction)/1E6, 'MB')

nc = nc.all[nc.all$jabba_event == 'none' & !is.na(nc.all$distance_to_nearest_junction), ]
message('Median non-clustered non-SV mutation distance to nearest breakpoint: ', median(nc$distance_to_nearest_junction)/1E6, 'MB')

nc = nc.all[nc.all$jabba_event == 'ecDNA', ]

message('\nNonclustered mutations on ecDNA: ', nrow(nc))

vaf.gt.80 = sum(nc$vaf > HIGH_VAF_CUTOFF)
vaf.gt.80.pct = round(vaf.gt.80 / nrow(nc), 3) * 100

vaf.lt.50 = sum(nc$vaf < LOW_VAF_CUTOFF)
vaf.lt.50.pct = round(vaf.lt.50 / nrow(nc), 3) * 100


message('\nNonclustered on ecDNA with VAF > ',HIGH_VAF_CUTOFF,': ', vaf.gt.80, ' (', vaf.gt.80.pct,'%)')
message('Nonclustered on ecDNA with VAF < ',LOW_VAF_CUTOFF,': ', vaf.lt.50, ' (', vaf.lt.50.pct,'%)')



#########################################
## Nonclustered SBS31/SBS35            ##
## Only patients highlighted in figure ##
#########################################

nc = nc.all[nc.all$sample == 'PM264-Z1-1-Case-WGS' & nc.all$jabba_event == 'ecDNA' & nc.all$assigned_signature %in% c('SBS31', 'SBS35'), ]

message('\nPM264-Z1-1-Case-WGS SBS31/35 nonclustered mutations on ecDNA: ', nrow(nc))

vaf.gt.80 = sum(nc$vaf > HIGH_VAF_CUTOFF)
vaf.gt.80.pct = round(vaf.gt.80 / nrow(nc), 3) * 100

vaf.lt.50 = sum(nc$vaf < LOW_VAF_CUTOFF)
vaf.lt.50.pct = round(vaf.lt.50 / nrow(nc), 3) * 100


message('\nPM264-Z1-1-Case-WGS SBS31/35 nonclustered on ecDNA with VAF > ',HIGH_VAF_CUTOFF,': ', vaf.gt.80, ' (', vaf.gt.80.pct,'%)')
message('PM264-Z1-1-Case-WGS SBS31/35 nonclustered on ecDNA with VAF < ',LOW_VAF_CUTOFF,': ', vaf.lt.50, ' (', vaf.lt.50.pct,'%)')



##############################
## Nonclustered SBS31/SBS35 ##
##############################
nc = nc.all[!is.na(nc.all$jabba_event) & nc.all$jabba_event == 'ecDNA' & nc.all$assigned_signature %in% c('SBS31', 'SBS35'), ]

message('\nSBS31/35 nonclustered mutations on ecDNA: ', nrow(nc))

vaf.gt.80 = sum(nc$vaf > HIGH_VAF_CUTOFF)
vaf.gt.80.pct = round(vaf.gt.80 / nrow(nc), 3) * 100

vaf.lt.50 = sum(nc$vaf < LOW_VAF_CUTOFF)
vaf.lt.50.pct = round(vaf.lt.50 / nrow(nc), 3) * 100


message('\nSBS31/35 nonclustered on ecDNA with VAF > ',HIGH_VAF_CUTOFF,': ', vaf.gt.80, ' (', vaf.gt.80.pct,'%)')
message('SBS31/35 nonclustered on ecDNA with VAF < ',LOW_VAF_CUTOFF,': ', vaf.lt.50, ' (', vaf.lt.50.pct,'%)')



#############################
## Nonclustered SBS2/SBS13 ##
#############################
nc = nc.all[!is.na(nc.all$jabba_event) & nc.all$jabba_event == 'ecDNA'& nc.all$assigned_signature %in% c('SBS2', 'SBS13'), ]

message('\nSBS2/13 nonclustered mutations on ecDNA: ', nrow(nc))

vaf.gt.80 = sum(nc$vaf > HIGH_VAF_CUTOFF)
vaf.gt.80.pct = round(vaf.gt.80 / nrow(nc), 3) * 100

vaf.lt.50 = sum(nc$vaf < LOW_VAF_CUTOFF)
vaf.lt.50.pct = round(vaf.lt.50 / nrow(nc), 3) * 100


message('\nSBS2/13 nonclustered on ecDNA with VAF > ',HIGH_VAF_CUTOFF,': ', vaf.gt.80, ' (', vaf.gt.80.pct,'%)')
message('SBS2/13 nonclustered on ecDNA with VAF < ',LOW_VAF_CUTOFF,': ', vaf.lt.50, ' (', vaf.lt.50.pct,'%)')
