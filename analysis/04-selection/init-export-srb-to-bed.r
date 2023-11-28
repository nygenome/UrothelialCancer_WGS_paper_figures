## Run fishhook to check for significantly recurrent breakpoints
libs = c('optparse','gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)


## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),   type='character', help="Output of annotate-fishhook-results.r"),
  make_option(c("-f", "--fdr"),       type='numeric', help="FDR cutoff"),
  make_option(c("-d", "--distance"),  type='numeric', help="Distance to nearest gene cutoff", default=5E5),
  make_option(c("-o", "--out_file"),  type='character', help="Output BED"))
opt = parse_args(OptionParser(option_list=option_list))


## Read data, filter
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')
dta = dta[dta$fdr < opt$fdr, ]
dta$nearest[dta$nearest_distance >= opt$distance] = ''


message('Total hits: ', nrow(dta))
message('Total hits within ', opt$distance,' of a CGC gene: ', sum(dta$nearest != ''))
table(dta$nearest)

if (nrow(dta)) {

  dta$cytoband = paste0(dta$chr, dta$cytoband)
  dta$nearest = paste0(dta$cytoband, ' (',dta$nearest,')')
  dta$nearest = gsub(' ()', '', dta$nearest, fixed=T)


}



## Collapse adjacent hits
dta = makeGRangesFromDataFrame(dta, keep.extra.columns=T)
dta = gr.reduce(dta, by='nearest')


## Select columns and write out 
if (!is.null(dta)) {

  dta = as.data.frame(dta)[ ,c('seqnames', 'start', 'end', 'nearest')]
  print(dta)
  write.table(dta, opt$out_file, row.names=F, col.names=F, quote=F, sep='\t')

} else {

  cmd = paste0('touch ', opt$out_file)
  system(cmd, intern=T)

}
message(opt$out_file)
