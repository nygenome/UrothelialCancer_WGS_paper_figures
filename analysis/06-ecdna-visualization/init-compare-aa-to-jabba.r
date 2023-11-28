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
## Compare AmpliconArchitect results to Jabba results
libs = c('optparse', 'GenomicRanges', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


JABBA_AMPS = c('bfb', 'dm', 'tyfonas')


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



## Convert a footrpint table into a single GRanges
fp.table.to.gr = function(x) {

  x$id = make.unique(x$type)

  res = mapply(FUN=fp2gr, 
               x=x$footprint, 
               id=x$id, 
               type=x$type, 
               SIMPLIFY=F)

  res = Reduce(c, res)

  return(res)

}



## Get arguments
option_list = list(
  make_option(c("-t", "--tn_file"),         type='character', help="Tumor-normal pairs file"),
  make_option(c("-a", "--in_file_aa"),      type='character', help="Table output of init-summarize-aa-results.r"),
  make_option(c("-j", "--in_dir_jabba"),    type='character', help="AmpliconArchitect run flag"),
  make_option(c("-o", "--out_file"),        type='character', help="Output TSV"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor normal pairs 
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','gender','patient'))
tn$pair_name = paste0(tn$tumor,'--',tn$normal)


## Read AA summary as GRanges, split on event call
aa = makeGRangesFromDataFrame(read.csv(opt$in_file_aa, h=T, stringsAsFactors=F, sep='\t'), keep.extra=T)

aa$sample = gsub('.*_(?=PM)', '', aa$id, perl=T)
aa = split(aa, aa$id)


## Construct table to hold results (1 row per AA event)
res = data.frame(aa_id=names(aa),
                 sample=unname(sapply(aa, function(x) unique(x$sample))),
                 class=unname(sapply(aa, function(x) unique(x$class))),
                 cycles=unname(sapply(aa, function(x) unique(x$cycles))),
                 width=unname(sapply(aa, function(x) sum(width(x)))),
                 intervals=unname(sapply(aa, function(x) length(x))),
                 fp=NA,
                 overlapping_jabba_event=NA,
                 overlapping_jabba_event_fp=NA, 
                 jabba_event_width = NA,
                 jabba_event_overlap_bp = NA,
                 jabba_event_overlap_fraction = NA,
                 aa_event_overlap_fraction = NA, 
                 reciprocal_overlap=NA)



## Read jabba events 
tn$in_files = paste0(opt$in_dir_jabba, '/', tn$pair_name, '/jabba.events.genes.tab')
tn = tn[file.info(tn$in_files)$size > 0, ]
jba = lapply(tn$in_files, read.csv, h=T, stringsAsFactors=F, sep='\t')
jba = lapply(jba, fp.table.to.gr)
names(jba) = tn$pair_name



## For each AA event, compute overlap statistics with jabba 
for (i in 1:nrow(res)) {

  aa.event = aa[[res$aa_id[i]]]
  jba.sample = jba[[res$sample[i]]]

  ## Filter for jabba events overlapping this amplicon any %
  id.keep = jba.sample$id[jba.sample %^% aa.event]
  if (length(id.keep) == 0) {
    next
  } 
  jba.sample = jba.sample[jba.sample$id %in% id.keep]

  ## Compute bp overlap 
  jba.sample$bp_overlap = jba.sample %o% aa.event

  ## Sum within each event call
  event.bp.overlap = tapply(jba.sample$bp_overlap, jba.sample$id, sum)
  event.size = tapply(width(jba.sample), jba.sample$id, sum)
  event.frac.overlap = event.bp.overlap / event.size
  event.reciprocal.overlap = (2 * event.bp.overlap) / (sum(width(aa.event)) + event.size)
  aa.event.overlap.fraction = event.bp.overlap / sum(width(aa.event))

  ## Concatenate results
  res$overlapping_jabba_event[i] = paste(names(event.bp.overlap), collapse=',')
  res$jabba_event_width[i] = paste(event.size, collapse=',')
  res$jabba_event_overlap_bp[i] = paste(event.bp.overlap, collapse=',')
  res$jabba_event_overlap_fraction[i] = paste(event.frac.overlap, collapse=',')
  res$aa_event_overlap_fraction[i] = paste(aa.event.overlap.fraction, collapse=',')
  res$reciprocal_overlap[i] = paste(event.reciprocal.overlap, collapse=',')
  res$overlapping_jabba_event_fp[i] = do.call(function(...) paste(..., collapse='|'), lapply(split(jba.sample, jba.sample$id), function(x) paste(grl.string(x), collapse=',')))
  res$fp[i] = paste(grl.string(aa.event), collapse=',')

}



## Select jabba amp that best overlaps 
res$overlapping_jabba_amp = NA
res$overlapping_jabba_amp_size = NA
res$overlapping_jabba_amp_id = NA
res$overlapping_jabba_amp_ro = NA
res$overlapping_jabba_event_amp_fp = NA
res$jabba_amp_overlap_fraction = NA
res$aa_event_amp_overlap_fraction = NA

for (i in 1:nrow(res)) {

  events = gsub('\\.[0-9]+', '', res$overlapping_jabba_event[i])
  events = unlist(strsplit(events, ','))

  event.ids = unlist(strsplit(res$overlapping_jabba_event[i], ','))
  event.sizes = as.numeric(unlist(strsplit(res$jabba_event_width[i], ',')))
  event.ro = as.numeric(unlist(strsplit(res$reciprocal_overlap[i], ',')))
  event.fp = unlist(strsplit(res$overlapping_jabba_event_fp[i], '|', fixed=T))

  event.overlap = as.numeric(unlist(strsplit(res$jabba_event_overlap_fraction[i], ',')))
  aa.overlap = as.numeric(unlist(strsplit(res$aa_event_overlap_fraction[i], ',')))


  if (any(JABBA_AMPS %in% events)) {

    if (sum(JABBA_AMPS %in% events) > 1) {
      message('More than one overlapping jabba amplification. Choosing event with greatest reciprocal overlap')
      print(res[i, ])
    }

    event.idx = which(events %in% JABBA_AMPS)
    events = events[event.idx]
    event.ids = event.ids[event.idx]
    event.sizes = event.sizes[event.idx]
    event.ro = event.ro[event.idx]
    event.overlap = event.overlap[event.idx]
    aa.overlap = aa.overlap[event.idx]

    res$overlapping_jabba_amp[i] = events[which.max(event.ro)]
    res$overlapping_jabba_amp_id[i] = event.ids[which.max(event.ro)]
    res$overlapping_jabba_amp_size[i] = event.sizes[which.max(event.ro)]
    res$overlapping_jabba_amp_ro[i] = event.ro[which.max(event.ro)]
    res$overlapping_jabba_event_amp_fp[i] = event.fp[which.max(event.ro)]

    res$jabba_amp_overlap_fraction[i] = event.overlap[which.max(event.ro)]
    res$aa_event_amp_overlap_fraction[i] = aa.overlap[which.max(event.ro)]

  }


}


## Write result
write.table(res, opt$out_file, row.names=F, col.names=T, quote=F, sep='\t')
message(opt$out_file)
