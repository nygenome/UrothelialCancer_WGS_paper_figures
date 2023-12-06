#!/nfs/sw/R/R-4.0.0/bin/Rscript
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: Timothy R. Chu 

################################################################# /COPYRIGHT ###
################################################################################
libs = c('optparse', 'deconstructSigs')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)

option_list = list(
  make_option(c("-d","--dir"), type="character", default=NULL, help="[REQUIRED either -d or -f] vcf file directory", metavar="character"),
  make_option(c("-f","--file"), type="character", default=NULL, help="[REQUIRED either -d or -f] single vcf file", metavar="character"),
  make_option(c("-r","--ref"), type="character", default="hg38", help="reference: [hg19,hg38]", metavar="character"),
  make_option("--highconf", type="logical", default=FALSE, action="store_true", help="vcf file directory", metavar="logical"),
  make_option(c("-o","--output"), type="character", default=NULL, help="[REQUIRED] output prefix", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


if (!is.null(opt$dir)){
  files = dir(opt$dir, pattern="*.snv.txt")
} else if (!is.null(opt$file)){
  files = c(opt$file)
} else {
  print("error. -d or -f option is REQUIRED")
}


if(opt$ref == "hg38"){
  library(BSgenome.Hsapiens.UCSC.hg38)
  load("signatures.genome.cosmic.v3.2.march2021.grch38.rda")
  sig_ref = signatures.genome.cosmic.v3.2.march2021.grch38
  bsg_ref = BSgenome.Hsapiens.UCSC.hg38
} else if (opt$ref =="hg19"){
  library(BSgenome.Hsapiens.UCSC.hg19,lib='/gpfs/commons/home/tchu/R')
  load("signatures.genome.cosmic.v3.2.march2021.grch37.rda")
  sig_ref = signatures.genome.cosmic.v3.2.march2021.grch37
  bsg_ref = BSgenome.Hsapiens.UCSC.hg19
} else{
  print("error. ref must be hg38 or hg19")
}

sbs96.input = data.frame()
tmb = c()
for(i in files){
  print(i)
  input = read.table(paste(opt$dir,i,sep=""), fill=T, col.names=c('chr','start','ref','alt','vaf'), sep="\t",header=T)
  samplename = strsplit(i,"/")[[1]][length(strsplit(i,"/")[[1]])]
  samplename = gsub(".snv.txt","",samplename)
  input[,'SAMPLE'] = samplename
  if(opt$highconf){
    input = input[grepl("HighConfidence",input[,"INFO"]),]
  }
  sigs.input <- mut.to.sigs.input(input, 
                                  sample.id = "SAMPLE", 
                                  chr = "chr", 
                                  pos = "start", 
                                  ref = "ref", 
                                  alt = "alt",
                                  bsg = bsg_ref,
                                  sig.type = 'SBS')
  sbs96.input = rbind(sbs96.input,sigs.input)
  tmb = append(tmb, nrow(input))
}
  

sbs96.output = data.frame()
sbs96.tumor = data.frame()
sbs96.product = data.frame()
sbs96.diff = data.frame()
for(samplename in rownames(sbs96.input)){
  sbssigs = whichSignatures(tumor.ref=sbs96.input, signatures.ref=sig_ref, sample.id=samplename, contexts.needed=T, signature.cutoff = 0.00, tri.counts.method = 'default')
  sbs96.output = rbind(sbs96.output, sbssigs$weights)
  sbs96.tumor = rbind(sbs96.tumor, sbssigs$tumor)
  sbs96.product = rbind(sbs96.product, sbssigs$product)
  sbs96.diff = rbind(sbs96.diff, sbssigs$diff)
}
sbs96.tmb.output = sweep(sbs96.output, 1, tmb, FUN="*")

sbs96.output = cbind(Sample = rownames(sbs96.output), sbs96.output)
sbs96.tmb.output = cbind(Sample = rownames(sbs96.tmb.output), sbs96.tmb.output)
sbs96.tumor = cbind(Sample = rownames(sbs96.tumor), sbs96.tumor)
sbs96.product = cbind(Sample = rownames(sbs96.product), sbs96.product)
sbs96.diff = cbind(Sample = rownames(sbs96.diff), sbs96.diff)


write.table(sbs96.output, file=paste(opt$output,".txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
write.table(sbs96.tmb.output, file=paste(opt$output,".counts.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
write.table(sbs96.tumor, file=paste(opt$output,".input.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
write.table(sbs96.product, file=paste(opt$output,".reconstructed.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
write.table(sbs96.diff, file=paste(opt$output,".diff.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)


