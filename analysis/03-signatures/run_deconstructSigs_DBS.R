library(deconstructSigs)
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg38)

option_list = list(
  make_option(c("-d","--dir"), type="character", default=NULL, help="[REQUIRED] vcf file directory", metavar="character"),
  make_option("--highconf", type="logical", default=FALSE, action="store_true", help="vcf file directory", metavar="logical"),
  make_option(c("-o","--output"), type="character", default=NULL, help="[REQUIRED] output", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(opt$version == "v3.4"){
  load("signatures.dbs.cosmic.v3.4.oct2023.grch38.rda")
  sig_ref = signatures.dbs.cosmic.v3.4.oct2023.grch38
}
else if(opt$version=="v3.2"){
  load("signatures.dbs.cosmic.v3.2.march2021.grch38.rda")
  sig_ref = signatures.dbs.cosmic.v3.2.march2021.grch38
}else{
  print("error. version must be 'v3.2' or 'v3.4")
}

files = dir(opt$dir, pattern="*.vcf")

dbs78.input = data.frame()
for(i in files){
  print(i)
  input = read.table(paste(opt$dir,i,sep=""), fill=T, col.names=c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR'))
  input[,'SAMPLE'] = gsub(".vcf","",i)
  if(opt$highconf){
    input = input[grepl("HighConfidence",input[,"INFO"]),]
  }
  sigs.input <- mut.to.sigs.input(input, 
                                  sample.id = "SAMPLE", 
                                  chr = "CHROM", 
                                  pos = "POS", 
                                  ref = "REF", 
                                  alt = "ALT",
                                  bsg = BSgenome.Hsapiens.UCSC.hg38,
                                  sig.type = 'DBS')
  dbs78.input = rbind(dbs78.input,sigs.input)
}


dbs78.output = data.frame()
for(samplename in rownames(dbs78.input)){
  dbssigs = whichSignatures(tumor.ref=dbs78.input, signatures.ref=sig_ref, sample.id=samplename, contexts.needed=T, signature.cutoff = 0)
  dbs78.output = rbind(dbs78.output, dbssigs$weights)
}
dbs78.output = cbind(Sample = rownames(dbs78.output), dbs78.output)


write.table(dbs78.output, file=opt$output, sep="\t", quote=FALSE, row.names=FALSE)

warnings()
