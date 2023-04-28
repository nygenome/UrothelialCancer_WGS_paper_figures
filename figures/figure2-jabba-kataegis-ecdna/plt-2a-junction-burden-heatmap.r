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
## Plot junction burden heatmap

libs = c('optparse', 'ComplexHeatmap', 'circlize', 'dendextend', 'ggplot2', 'RColorBrewer')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



## Generate named vector for use in metadata color tracks
build.colormap = function(var, col) {
  
  var.col = structure(col, names=as.character(var))
  var.col = na.omit(var.col[!duplicated(as.character(var))])
  return(var.col)
  
}



## Oncoprint-style sorting -- borrowed from somewhere within 
## the ComplexHeatmap source
oncoprint_column_order = function(count_matrix) {
  scoreCol = function(x) {
    score = 0
    x = as.numeric(x > 0)

    for(i in 1:length(x)) {
      if(x[i]) {
        score = score + 2^(length(x)-i*1/x[i])
      }
    }
    return(score)
  }
  scores = apply(count_matrix[, ,drop = FALSE], 2, scoreCol)
  order(scores, decreasing=TRUE)
}



AMP = c('tyfonas', 'dm', 'bfb')
MUT_COL = c(`NONE`='white', `MODERATE`='firebrick1', `HIGH`='firebrick')



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),       type='character', help="Input junction summary by type"),
  make_option(c("-p", "--pp"),            type='character', help="Tab-delimited headerless file with columns tumor,normal,patient,site,purity,ploidy"),
  make_option(c("-m", "--metadata"),      type='character', help="Metadata file"),
  make_option(c("-f", "--fga"),           type='character', help="Tab-delimited file with columns sample,FGA"),
  make_option(c("-x", "--tp53"),          type='character', help="Tab-delimited file with columns sample,TP53_status"),
  make_option(c("-o", "--out_file"),      type='character', help="Output SVG"))
opt = parse_args(OptionParser(option_list=option_list))


## Read tumor-normal pairs
pp = read.csv(opt$pp, h=F, stringsAsFactors=F, col.names=c('tumor','purity','ploidy'))
pp$name = pp$tumor

## Read in metadata, align with tumor-normal pairs
mta = read.csv(opt$metadata, h=T, stringsAsFactors=F, sep='\t')
mta = mta[match(pp$tumor, mta$tumor), ]
mta$sample = paste0(mta$tumor,'--',mta$normal)

## Read FGA, add to metadata
fga = read.csv(opt$fga, h=T, stringsAsFactors=F, sep='\t')
mta$fga = fga$FGA[match(mta$sample, fga$sample)]

## Read TP53/TERTp status, add to metadata
tp53 = read.csv(opt$tp53, h=T, stringsAsFactors=F, sep='\t')
mta$tp53 = tp53$TP53_status[match(mta$sample, tp53$sample)]

## Read data, select tumor-normal pairs, convert to matrix
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F)
rownames(dta) = gsub('--.*', '', dta$sample)
dta = dta[pp$tumor, setdiff(colnames(dta), 'sample')]
njunc = dta$njunc
dta$dm = dta$dm + dta$cpxdm 
dta = as.matrix(dta[, setdiff(colnames(dta), c('unclassified', 'njunc', 'cpxdm'))])


## Take ln(x+1) and transpose
dta.orig = as.data.frame(dta)
dta = t(log1p(dta))


## Remove simple events
events.rm = c('del', 'dup', 'inv', 'tra', 'invdup')
dta = dta[!rownames(dta) %in% events.rm, ]

## Update event names
rownames(dta) = gsub('invdup','Inverted Duplication',rownames(dta))
rownames(dta) = gsub('tra','Translocation',rownames(dta))
rownames(dta) = gsub('inv','Inversion',rownames(dta))
rownames(dta) = gsub('del','Deletion',rownames(dta))
rownames(dta) = gsub('dup','Duplication',rownames(dta))
rownames(dta) = gsub('chromoplexy','Chromoplexy',rownames(dta))
rownames(dta) = gsub('rigma','Rigma',rownames(dta))
rownames(dta) = gsub('pyrgo','Pyrgo',rownames(dta))
rownames(dta) = gsub('chromothripsis','Chromothripsis',rownames(dta))
rownames(dta) = gsub('tic','TIC',rownames(dta))
rownames(dta) = gsub('qrp','QRP',rownames(dta))
rownames(dta) = gsub('cpxdm','Complex DM',rownames(dta))
rownames(dta) = gsub('bfb','BFB',rownames(dta))
rownames(dta) = gsub('dm','Double Minute',rownames(dta))
rownames(dta) = gsub('tyfonas','Tyfonas',rownames(dta))



#################
## Annotations ##
#################

mta$patient_original = mta$patient
mta$patient[mta$patient_color == '#FFFFFF'] = 'Singleton'
mta$platinum_chemotherapy = factor(mta$platinum_chemotherapy, levels=c('Pre-chemo', 'Post-chemo'))
mta$localization = factor(mta$localization, levels=c('Primary', 'Metastatic'))

top.ano = ComplexHeatmap::HeatmapAnnotation(#`FGA`=ComplexHeatmap::anno_barplot(mta$fga, baseline=0, ylim=0:1, gp=gpar(fill='black')),
                                            #Ploidy=ComplexHeatmap::anno_barplot(pp$ploidy, baseline=2, gp=gpar(fill='black')),
                                            `Chemotherapy`=mta$platinum_chemotherapy,
                                            `Localization`=mta$localization,
                                            `Somatic TP53 SNV/INDEL Impact`=mta$tp53,
                                            `Patient`=mta$patient,
                                            col=list(`Localization`=build.colormap(mta$localization, mta$localization_color),
                                                     `Chemotherapy`=build.colormap(mta$platinum_chemotherapy, mta$platinum_chemotherapy_color),
                                                     `Family History`=build.colormap(mta$family_history, mta$family_history_color),
                                                     `Tobacco`=build.colormap(mta$smoking_status, mta$smoking_status_color),
                                                     `Patient`=build.colormap(mta$patient, mta$patient_color),
                                                     `Somatic TP53 SNV/INDEL Impact`=MUT_COL),
                                            gap=unit(c(0, 0, 0, 0, 0), "mm"),
                                            gp = gpar(col="black"))


## Generate colors
quantiles = seq(0,1,by=0.2)
tmp = brewer.pal(n = length(quantiles), name = "YlOrRd")
tmp = c('grey95', tmp)
burden.col = circlize::colorRamp2(c(0, quantile(dta[dta > 1], quantiles)), tmp)


## Compute dendrogram
k = 5
dend = hclust(dist(t(dta), 'euclidean'), method='complete')
dend = dendextend::color_branches(dend, k=k)


## Total junction summary
njunc = matrix(njunc, ncol=length(njunc), nrow=1)
rownames(njunc) = 'Total SV\nBurden'
colnames(njunc) = mta$sample_id[match(colnames(dta), mta$tumor)]


## Manually order rows
row.order = c('Tyfonas', 'BFB', 'Double Minute', 'Chromothripsis', 'Chromoplexy', 'TIC', 'QRP', 'Rigma', 'Pyrgo')
row.split = factor(c('Putative ecDNA', 'Putative ecDNA', 'Putative ecDNA', '', '', '', '', '', ''), levels=c('Putative ecDNA', ''))
dta = dta[row.order, ]

col.order = oncoprint_column_order(rbind(dta, njunc))
dta = dta[, col.order]
njunc = njunc[, col.order, drop=F]


## Generate heatmap
njunc.row = Heatmap(log1p(njunc),
                    col=burden.col,
                    rect_gp = gpar(col = "black", lwd=1),
                    border=T,
                    row_names_side='left',
                    heatmap_legend_param = list(title = "Total SV Burden"),
                    show_column_names=T,
                    cluster_rows=F,
                    cluster_columns=F,
                    show_heatmap_legend=F,
                    cell_fun = function(j, i, x, y, w, h, col) grid.text(njunc[j], x, y, gp=gpar(fontsize=6)))

main.heatmap = ComplexHeatmap::Heatmap(dta,
                                      col=burden.col,
                                      border=T,
                                      row_split=row.split,
                                      top_annotation=top.ano,
                                      row_names_side='left',
                                      column_title=NULL,
                                      heatmap_legend_param = list(title = "ln(SV burden)"),
                                      show_column_names=T,
                                      show_row_names=T,
                                      show_column_dend=F,
                                      show_row_dend=F, 
                                      cluster_rows=F,
                                      cluster_columns=F,
                                      show_heatmap_legend=T)


## Draw heatmap
svg(opt$out_file, height=10, width=18)

draw(main.heatmap %v% njunc.row)

dev.off()
message(opt$out_file)



## Export table with all relevant info  
colnames(dta.orig) = paste0(colnames(dta.orig), '_junctions')
dta.orig$total_junction_count = unlist(njunc[1, ])
mta$patient = mta$patient_original

jba.cols = colnames(dta.orig)
mta.cols = c('sample_id', 'patient', 'platinum_chemotherapy', 'localization', 'tp53', 'fga')
dta.orig[, mta.cols] = mta[match(rownames(dta.orig), mta$tumor), mta.cols]
dta.orig[, c('purity','ploidy')] = pp[match(rownames(dta.orig), pp$name), c('purity','ploidy')]

dta.orig = dta.orig[, c(mta.cols, 'purity','ploidy', jba.cols)]

colnames(dta.orig) = tolower(colnames(dta.orig))
out.file.txt = gsub('\\.svg', '.txt', opt$out_file)
write.table(dta.orig, out.file.txt, row.names=F, col.names=T, quote=F, sep='\t')
message(out.file.txt)
