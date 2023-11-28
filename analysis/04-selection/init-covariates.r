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
## Initialize all covariates used in the FishHook analysis 
libs = c('optparse', 'GenomicRanges', 'gUtils', 'data.table', 'BSgenome.Hsapiens.UCSC.hg38', 'rtracklayer')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Get arguments
option_list = list(
  make_option(c("-c", "--chr_len"),       type='character', help="Chromosome lengths"),
  make_option(c("-e", "--encode_dir"),    type='character', help="Dir with ENCODE data"),
  make_option(c("-r", "--rep_timing"),    type='character', help="Replication timing bigWig file"),
  make_option(c("-f", "--fragile_sites"), type='character', help="Fragile sites as annotated in HGNC biomart"),
  make_option(c("-y", "--cytoband"),      type='character', help="UCSC Cytoband file"),
  make_option(c("-m", "--repeatmasker"),  type='character', help="UCSC RepeatMasker file"),
  make_option(c("-b", "--bin_size"),      type='numeric',   help="Bin size", default=1E5),
  make_option(c("-s", "--step_size"),     type='numeric',   help="Bin step size", default=5E4),
  make_option(c("-o", "--out_dir"),       type='character', help="Output dir"))
opt = parse_args(OptionParser(option_list=option_list))



## Read chromosomes of interest
chr.len = read.csv(opt$chr_len, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','end'))
chr.len$start = 1
chr.len = chr.len[!chr.len$chr %in% c('chrM', 'chrY'), ]
chr.len = makeGRangesFromDataFrame(chr.len)
bins = gr.tile(chr.len, opt$bin_size)

if (opt$step_size < opt$bin_size) {

  message('Using step size: ', opt$step_size)

  ## Original bins, plus those shifted right by the step size
  ## Make sure to exclude bins that fall off the edge of the chromosome
  bins = c(bins, bins %+% opt$step_size)
  bins = bins[bins %O% chr.len == 1]

}



############################
## Nucleotide frequencies ##
############################

nuc = bins

message('Computing nucleotide frequency...')
mono = oligonucleotideFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, bins), width = 1, as.prob = TRUE)

message('Computing dinucleotide frequency...')
di = dinucleotideFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, bins), as.prob = TRUE)

message('Computing trinucleotide frequency...')
tri = trinucleotideFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, bins), as.prob = TRUE)

mcols(nuc) = cbind(mono, di, tri)

out.file = paste0(opt$out_dir,'/covariate_nuc_',opt$bin_size,'.rds')
saveRDS(nuc, out.file)



############################################
## Histone marks / DNase hypersensitivity ## 
############################################

## Borrowed from Will Liao's code: 
## /gpfs/commons/groups/compbio/projects/Project_Cornell_Collaboration_Bladder_Cancer/Project_KIM_14248_B01_SOM_WGS/compbio/noncoding/fishHook/fishHook.ipynb
hm.in.files = Sys.glob(paste0(opt$encode_dir,c("/*ChIP*narrowPeak", "/GRCh38.urinary_bladder.ENCFF823HYK.DNase-seq.p0.05_peaks.2020-10-24.narrowPeak")))
names(hm.in.files) = gsub("ChIP\\-seq\\.", "", gsub(".*ENC\\w+\\.(.+?)\\.p.*", "\\1", basename(hm.in.files)))

narrow.col.names =   col_names = c("chrom", "start", "end", "name", "score", "strand",
                                   "signalValue", "pValue", "qValue", "peak")

message('Processing histone marks...')
for (i in 1:length(hm.in.files)) {

  hm = names(hm.in.files)[i]
  narrow = read.csv(hm.in.files[i], h=F, stringsAsFactors=T, sep='\t', col.names=narrow.col.names)
  narrow = makeGRangesFromDataFrame(narrow)
  mcols(narrow)[[hm]] <- hm

  out.file = paste0(opt$out_dir,'/','covariate_',hm,'_',opt$bin_size,'.rds')
  saveRDS(narrow, out.file)
  message(out.file)

}



########################
## Replication timing ##
########################

message('Processing replication timing...')

rep.timing = rtracklayer::import(opt$rep_timing, format='bigWig')
colnames(mcols(rep.timing))[1] = 'replication_timing'

out.file = paste0(opt$out_dir,'/covariate_reptiming_',opt$bin_size,'.rds')
saveRDS(rep.timing, out.file)
message(out.file)



###################
## Fragile sites ##
###################

message('Processing fragile sites...')

# ## Read fragile sites and cytoband
fragile = read.csv(opt$fragile_sites, h=T, stringsAsFactors=F, sep='\t')
cyto = read.csv(opt$cytoband, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','start','end','cytoband','stain'))
cyto = makeGRangesFromDataFrame(cyto, keep.extra.columns=T)

cyto$cytoband = paste0(gsub('^chr','', as.character(seqnames(cyto))), cyto$cytoband)

## Subset to fragile sites only
cyto = granges(cyto[cyto$cytoband %in% fragile$Chromosome.location])
cyto$fragile_site = 'fragile'

out.file = paste0(opt$out_dir,'/covariate_fragile_',opt$bin_size,'.rds')
saveRDS(cyto, out.file)
message(out.file)



##################
## RepeatMasker ##
##################

message('Processing RepeatMasker annotations')

rmsk.keep = c('LINE', 'SINE', 'LTR', 'Simple_repeat', 'DNA')

## Read repeatmasker annotations
rmsk = fread(cmd=paste('zcat', opt$repeatmasker), header=F, stringsAsFactors=F, sep='\t')[,c(6:8, 12)]
colnames(rmsk) = c('chr', 'start', 'end', 'name')
rmsk = makeGRangesFromDataFrame(rmsk, keep.extra.columns=T)

## Output each RM type individually
for (i in rmsk.keep) {

  res = granges(rmsk[rmsk$name == i])
  mcols(res)[[i]] = i

  out.file = paste0(opt$out_dir,'/','covariate_',i,'_',opt$bin_size,'.rds')
  saveRDS(res, out.file)
  message(out.file)

}
