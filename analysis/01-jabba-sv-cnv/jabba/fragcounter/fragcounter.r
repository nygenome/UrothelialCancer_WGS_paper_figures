## FragCounter is developed by the Marcin Imielinski group. For full attribution and usage
## details, please see https://github.com/mskilab-org/fragCounter
library(optparse)
options(bitmapType='cairo')

if (!exists('opt'))
    {
      option_list = list(
          make_option(c("-b", "--bam"),
                      type = "character", help = "Path to .bam file"),
          make_option(c("-c", "--cov"),
                      type = "character", help = "Path to existing coverage rds or bedgraph"),
          make_option(c("-m", "--midpoint"),
                      type = "character", default = "TRUE",
                      help = "If TRUE only count midpoint if FALSE then count bin footprint of every fragment interval"),
          make_option(c("-w", "--window"),
                      type = "integer", default = 1000, help = "Window / bin size"),
        make_option(c("-d", "--gcmapdir"), type = "character", help = "Mappability / GC content dir"),
        make_option(c("-q", "--minmapq"), type = "integer", default = 1, help = "Minimal map quality"),
        make_option(c("-p", "--paired"), type = "logical", default = TRUE, help = "Is paired"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "Directory to dump output into"),
        make_option(c("-l", "--libdir"), type = "character", default = paste(Sys.getenv('GIT_HOME'), 'isva', sep = '/'), help = "Directory containing this R file")
      )

      parseobj = OptionParser(option_list=option_list)
      opt = parse_args(parseobj)

      if (is.null(opt$libdir) | (is.null(opt$bam) & is.null(opt$cov)))
        stop(print_help(parseobj))
      
      ## keep record of run
      writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
      saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }


print(opt)

print(Sys.getenv('LD_LIBRARY_PATH'))

print(.libPaths())
# .libPaths(c(.libPaths(), opt$rlibdir))

require(GenomicRanges)
library(skitools)
library(gUtils)
library(data.table)
library(rtracklayer)
library(bamUtils)

set.seed(42) #' twalradt Friday, Apr 20, 2018 01:05:40 PM

if (!is.null(opt$gcmapdir))
    if (!file.exists((opt$gcmapdir)))
        opt$gcmapdir = NULL
                   
if (is.null(opt$gcmapdir))
    {
        GC.WIG.DIR = paste(opt$libdir, 'gccontent/', sep = '/') ## hg19 gc content
        MAP.WIG.DIR = paste(opt$libdir, 'mappability/', sep = '/') ## hg19 mappability
    } else
        {
            GC.WIG.DIR = paste(opt$gcmapdir, 'gccontent/', sep = '/') ## hg19 gc content
            MAP.WIG.DIR = paste(opt$gcmapdir, 'mappability/', sep = '/') ## hg19 mappability
        }



correctcov_stub = function(out.rds, ## root of outfile which will be RDS .rds file and .bedgraph bedgraph files for uncorrected and corrected read counts
    cov.wig, ## wig file of coverage tiles of width W or pointer to rds file of sorted GRanges object or GRanges object
    mappability = 0.9,
    samplesize = 5e4,
    verbose = T,
    gc.wig.dir = GC.WIG.DIR, ## for tiles of width W will look here for a file named gc{W}.wig in this directory
    map.wig.dir = MAP.WIG.DIR ## for tiles of width W will look here for a file named map{W}.wig in this directory
    )
    {
        if (is.character(cov.wig))
            {
                if (grepl('(\\.bedgraph$)|(\\.wig$)|(\\.bw$)', cov.wig))
                    {
                        cov = import.ucsc(cov.wig)
                    }
                else if (grepl('(\\.rds)', cov.wig))
                    {
                        cov = gr.stripstrand(readRDS(cov.wig))
                        names(values(cov))[1] = 'score' ## will take the first value column as the (uncorrected) read count
                    }
                else
                    stop("Unsupported coverage format (should be UCSC bedgraph, wig, bw, or R Bioconductor GRanges .rds file")
            }
        else ## assume it is a sorted GRanges
            cov = gr.stripstrand(cov.wig)
                
        n = length(cov)
        
        wid = as.numeric(names(sort(-table(width(cov))))[1])

        gc.wig = paste(gc.wig.dir,'/gc', wid, '.wig', sep = '')
        map.wig = paste(map.wig.dir,'/map', wid, '.wig', sep = '')
        
        gc.rds = paste(gc.wig.dir,'/gc', wid, '.rds', sep = '')
        map.rds = paste(map.wig.dir,'/map', wid, '.rds', sep = '')

        cat('Loaded GC and mappability\n')
        
        
        if (!file.exists(gc.rds) | !file.exists(map.rds))
            {        
                if (!file.exists(gc.wig))
                    stop(sprintf('GC wig file %s not found, either directory is wrong or file needs to be generated for this particular window width', gc.wig))
                
                if (!file.exists(map.wig))
                    stop(sprintf('mappability wig file %s not found, either directory is wrong or file needs to be generated for this particular window width', map.wig))
                
                gc = import.ucsc(gc.wig)
                map = import.ucsc(map.wig)
            }
        else ## if we have rds files then let's use these to avoid using rtracklayer
            {
                gc = readRDS(gc.rds)
                map = readRDS(map.rds)
            }

        if (is.null(cov$score)) ## if $score field is empty then just assume that the first column of coverage is the "score" i.e. read count
            names(values(cov))[1] = 'score'

        cov = gr.sub(cov, 'chr', '')
        map = gr.sub(map, 'chr', '')
        gc = gr.sub(gc, 'chr', '')
        
        gc.str = gr.string(gc)
        map.str = gr.string(map)
        cov.str = gr.string(cov)
      
        all.str = intersect(intersect(map.str, cov.str), gc.str) ## in case we have mismatches in the ordering / genome definition

        
        cat(sprintf('length cov is %s, length gc is %s, length map is %s\n',
                length(cov),
                length(gc),
                length(gc)))
        
        print('about to fail?')
        
        print(class(match(all.str, map.str)))
        
        map = map[match(all.str, map.str)]
        
        print('failed?')
        
        cov = cov[match(all.str, cov.str)]
        gc = gc[match(all.str, gc.str)]
                
        
        
        
        if (length(cov) != length(gc) | length(gc) != length(map) | (length(cov))/n<0.1)
            stop('Mismatch / problem in cov, gc, or map definition.  Check if they come from the same width tiling')
    
        cov$reads = cov$score
        cov$gc = gc$score
        cov$gc[cov$gc<0] = NA
        cov$map = map$score
        cov$score = NULL

        rm(gc)
        rm(map)
        gc()
        
        cat('Synced coverage, GC, and mappability\n')

        imageroot = gsub('.rds$', '', out.rds)
        
        .correctReadCount = function(x, mappability = 0.9, samplesize = 5e4, verbose = T)
            {
                if (!is.data.frame(x))
                    x = as.data.frame(x)
                if (length(x$reads) == 0 | length(x$gc) == 0 | length(x$map) == 0)          
                    {
                        stop("Missing one of required columns: reads, gc, map")
                    }
                if (verbose) {
                    message("Applying filter on data...")
                }
                x$valid <- TRUE
                x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
                x$ideal <- TRUE
                routlier <- 0.01
                range <- quantile(x$reads[x$valid], prob = c(0, 1 - routlier), 
                                  na.rm = TRUE)
                doutlier <- 0.001
                domain <- quantile(x$gc[x$valid], prob = c(doutlier, 1 - 
                                                               doutlier), na.rm = TRUE)
                x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] | 
                            x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
                if (verbose) {
                    message("Correcting for GC bias...")
                }
                set <- which(x$ideal)
                select <- sample(set, min(length(set), samplesize))
                rough = loess(x$reads[select] ~ x$gc[select], span = 0.03)
                i <- seq(0, 1, by = 0.001)
                final = loess(predict(rough, i) ~ i, span = 0.3)

                if (!is.null(imageroot))
                    {
                        out.png = paste(imageroot, ifelse(grepl("/$", imageroot), '', '.'), 'og_gc_correction.png', sep = '')
                        if (verbose) {
                            cat("Dumping figure to", out.png, "\n")
                        }
                        png(out.png, height = 1000, width = 1000)
                        x2s = x[select, ]
                        df = data.frame(gc = seq(domain[1], domain[2], 0.001))
                        plot(x2s$gc, x2s$reads, col = alpha('black', 0.1), pch = 19, cex = 0.4, xlim = domain, ylim = range, ylab = 'signal before GC correction', xlab = 'gc');
                        lines(df$gc, predict(final, df$gc), col = 'red', lwd = 2)
                        dev.off()
                    }
                
                x$cor.gc <- x$reads/predict(final, x$gc)
                if (verbose) {
                    message("Correcting for mappability bias...")
                }
                coutlier <- 0.01
                range <- quantile(x$cor.gc[which(x$valid)], prob = c(0, 1 - coutlier), na.rm = TRUE)
                set <- which(x$cor.gc < range[2])
                select <- sample(set, min(length(set), samplesize))
                final = approxfun(lowess(x$map[select], x$cor.gc[select]))
                return(x$cor.gc/final(x$map))
            }

        cov = sort(gr.fix(cov))
        
        cat('Modified gc / mappability correction\n')
        cov$reads.corrected = multicoco(cov, numlevs = 1, base = max(10, 1e5/wid), mc.cores = 1, fields = c('gc', 'map'), iterative = T, mono = T, imageroot = imageroot)$reads.corrected
        gc()
        
        cat('Original gc / mappability correction\n')
        cov$reads.corrected.og = .correctReadCount(cov, mappability, samplesize, verbose)
 
        gc()
        
        out.og = paste(gsub('.rds$', '', out.rds), '.original.bw', sep = '')
        out.corr = paste(gsub('.rds$', '', out.rds), '.corrected.bw', sep = '')

        if (!is.null(tryCatch({library(rtracklayer); 'success'}, error = function(e) NULL)))
            {
                cov.og.out = cov
                cov.og.out$score =cov$reads
                cov.og.out$score[is.na(cov.og.out$score)] = -1
                cov.og.out = cov.og.out[width(cov.og.out)==wid] # remove any funky widths at end of chromosome
                export(cov.og.out[, 'score'], out.og, 'bigWig', dataFormat = 'fixedStep')
                
                cov.corr.out = cov
                cov.corr.out$score = cov$reads.corrected
                cov.corr.out$score[is.na(cov.corr.out$score)] = -1
                cov.corr.out = cov.corr.out[width(cov.corr.out)==wid] ## remove any funky widths at end of chromosome
                export(cov.corr.out[, 'score'], out.corr, 'bigWig', dataFormat = 'fixedStep')
            }
        
        saveRDS(cov, paste(gsub('.rds$', '', out.rds), '.rds', sep = ''))
 #       bedgraph = data.frame(chrom = as.character(seqnames(cov)), chromStart = start(cov), chromEnd = end(cov))
        
 #       write.table(cbind(bedgraph, score = cov$reads.corrected),  out.corr, sep = '\t', col.names = F, row.names = F, quote = F)
 #       write.table(cbind(bedgraph, score = cov$reads ), out.og, sep = '\t', col.names = F, row.names = F, quote = F)
    }

######################################
# multicoco
#
# multi-scale coverage correction
#
# Given gc and mappability coverage correction at k "nested" scales finds the coverage
# assignment at the finest scale that yields the best correction at every scale
#
# FUN = is a function that takes in a data frame / granges with
# $reads and other covariate functions $gc, $map and outputs a vector of corrected read counts
#
# cov is a constant with GRanges of coverage samples with (by default) fields $reads, $map, $gc 
#
# base = is the multiple with which to perform "numlevs" additional scales of correction
#
######################################
multicoco = function(cov,
    numlevs = 1, ## numbers of scales at which to correct
    base = 100, ## scale multipler
    fields = c("gc", "map"), # fields of gc which to use as covariates
    iterative = TRUE,
    presegment = TRUE, ## whether to presegment
    min.segwidth = 5e6, ## when presegmenting, min seg width
    mono = FALSE, ## just do single iteration at finest scale
    verbose = T,
    imageroot = NULL, ## optional file root to which to dump images of correction
    FUN = NULL, ## function with which to correct coverage (by default loess correction modified from HMMcopy that takes in granges with fields $reads and other fields specified in "fields"
    ..., ## additional args to FUN
    mc.cores = 1)
    {
        if (verbose)
           cat('Converting to data.table\n')

        WID = max(width(cov))

        library(data.table)
        cov.dt = gr2dt(cov)
        
        sl = structure(as.numeric(1:length(seqlevels(cov))), names = seqlevels(cov))       

        if (verbose)
            cat('Grouping intervals\n')

        ## compute level means
        ## lev 0 = raw data
        ## lev 1 = base-fold collapsed
        ## lev 2 = base^2-fold collapsed
        ## etc
        parentmap= list() ## data.tables that map lev labels at level k  to parent lev labels
        cov.dt[, lev0:=as.character(1:length(seqnames))]
        for (k in 1:numlevs)
            {
                if (verbose)
                    cat('Defining', base^k, 'fold collapsed ranges\n')
                cov.dt[, eval(paste("lev", k, sep = '')) := as.character(sl[seqnames] + as.numeric(Rle(as.character(1:length(start)), rep(base^k, length(start)))[1:length(start)])/length(start)), by = seqnames]
                parentmap[[k]] = data.table(parent = cov.dt[, get(paste("lev", k, sep = ''))], child = cov.dt[, get(paste("lev", k-1, sep = ''))], key = 'child')[!duplicated(child), ]
            }

        if (presegment) ## perform rough segmentation at highest level
            {
                seg = NULL
                sl = seqlengths(cov)
                if (verbose)
                    cat('Presegmenting at ', as.integer(WID*base^(numlevs)), ' bp scale \n')
                require(DNAcopy)
                tmp.cov = seg2gr(cov.dt[,list(chr = seqnames[1], start = min(start), end = max(end), strand = strand[1], reads = mean(reads, na.rm = T)), by = get(paste("lev", numlevs, sep = ''))][end>start, ], seqlengths = sl)
                ix = which(!is.na(values(tmp.cov)[, 'reads']))
                tmp = data.frame()

                tryCatch({
                    cna = CNA(log(values(tmp.cov)[, 'reads'])[ix], as.character(seqnames(tmp.cov))[ix], start(tmp.cov)[ix], data.type = 'logratio')
                    tmp = print(segment(smooth.CNA(cna), alpha = 1e-5, verbose = T))
                    tmp = tmp[!is.na(tmp$loc.start) & !is.na(tmp$chrom) & !is.na(tmp$loc.end), , drop = F]
                }, error = function(e) warning('DNACopy error moving on without segmenting'))
                
                if (nrow(tmp)>0)
                    {
                        seg = sort(seg2gr(tmp, seqlengths = sl))
                        seg = seg[width(seg)>min.segwidth] ## remove small segments
                        seg.dt = gr2dt(seg)
                        if (nrow(seg.dt)>0)
                            {
                                seg = seg2gr(seg.dt[, list(seqnames = seqnames,
                                    start = ifelse(c(FALSE, seqnames[-length(seqnames)]==seqnames[-1]), c(1, start[-1]), 1),
                                    end = ifelse(c(seqnames[-length(seqnames)]==seqnames[-1], FALSE), c(start[-1]-1, Inf), seqlengths(seg)[as.character(seqnames)]))], seqlengths = sl)
                                seg = gr.val(seg, tmp.cov, 'reads') ## populate with mean coverage
                                seg$reads = seg$reads/sum(as.numeric(seg$reads*width(seg))/sum(as.numeric(width(seg)))) ## normalize segs by weigthed mean (so these become a correction factor)
                            }
                        else
                            seg = NULL
                    }
                else
                    seg = NULL
            }
        else
            seg = NULL
        
        if (verbose)
            cat('Aggregating coverage within levels \n')
        
        ## list of data frames showing scales of increasing collapse

        cov.dt[, ix := 1:nrow(cov.dt)]
        
        cmd1 = paste('list(ix.start = ix[1], ix.end = ix[length(ix)], reads = mean(reads, na.rm = T),', paste(sapply(fields, function(f) sprintf("%s = mean(%s, na.rm = T)", f, f)), collapse = ','), ')', sep = '')
        
        cmd2 = paste('list(lab = lev0, reads,', paste(fields, collapse = ','), ', seqnames, start, end)', sep = '')

        if (mono)
            {
                if (verbose)
                    cat('Mono scale correction \n')
                 grs = list(cov.dt[, eval(parse(text=cmd2))])
                 numlevs = 1
             }
        else
            {
                grs = c( list(cov.dt[, eval(parse(text=cmd2))]),
                    lapply(1:numlevs, function(x)
                        {
                            if (verbose)
                                cat('Aggregating coverage in level', x,  '\n')
                            out = cov.dt[, eval(parse(text=cmd1)), keyby = list(lab = get(paste('lev', x, sep = '')))]
                            out[, ":="(seqnames = cov.dt$seqnames[ix.start], end = cov.dt$end[ix.start], start = cov.dt$start[ix.start])]
                            out[, ":="(ix.start= NULL, ix.end = NULL)]
                            return(out)
                        }))
            }
        
        setkey(grs[[1]], 'lab')
                               
        ## modified from HMMCopy to 
        ## (1) take arbitrary set of covariates, specified by fields vector
        ## (2) employ as input an optional preliminary (coarse) segmentation with which to adjust signal immediately prior to loess
        ## NOTE: (this only impacts the loess fitting, does not impose any segmentation on thed ata)
        ##
        if (is.null(FUN))            
            FUN = function(x, fields = fields, samplesize = 5e4, seg = NULL, ## seg is a Granges with meta data field $reads
                verbose = T, doutlier = 0.001, routlier = 0.01)
                {                    
                    if (!all(fields %in% names(x)))
                        stop(paste('Missing columns:', paste(fields[!(fields %in% names(x))], collapse = ',')))
                    
                    x$valid <- TRUE
                    for (f in fields)
                        {
                            x$valid[is.na(x[, f])] = FALSE
                            x$valid[which(is.infinite(x[, f]))] = FALSE
                        }
                   
                    if (verbose)
                        cat('Quantile filtering response and covariates\n')

                    range <- quantile(x$reads[x$valid], prob = c(routlier, 1 - routlier), na.rm = TRUE)

                    if (verbose)
                        cat(sprintf("Response min quantile: %s max quantile: %s\n", round(range[1],2), round(range[2],2)))
                    
                    domains = lapply(fields, function(f) quantile(x[x$valid, f], prob = c(doutlier, 1 - doutlier), na.rm = TRUE))
                    names(domains) = fields
                    
                    x$ideal <- x$valid
                    x$ideal[x$reads<=range[1] | x$reads>range[2]] = FALSE
                                        
                    for (f in fields)
                        x$ideal[x[, f] < domains[[f]][1] | x[, f] > domains[[f]][2]] = FALSE

                    if (verbose)
                        cat(sprintf('Nominated %s of %s data points for loess fitting\n', sum(x$ideal), nrow(x)))

                    set <- which(x$ideal)

                    if (length(set)<=10)
                        {
                            warning("Not enough samples for loess fitting - check to see if missing or truncated data?")
                            return(x$reads)
                        }
                    
                    for (f in fields)
                        {                            
                            if (verbose)
                                message(sprintf("Correcting for %s bias...", f))
                            
                            select <- sample(set, min(length(set), samplesize))

                            x2 = x[, c('reads', f)]
                            x2$covariate = x2[, f]

                            x2s = x2[select, ]
                            
                            if (!is.null(seg)) ## here we apply a prelmiinary segmentation to correct for large scale copy variation allow more power to reveal the covariate signal
                                {
                                    if (verbose)
                                        message('Applying preliminary segmentation prior to loess fitting')
                                    
                                    x.grs = gr.val(seg2gr(x[select, ], seqlengths = NULL), seg, 'reads')
                                    x2s$reads = x2s$reads/x.grs$reads
                                }
                            
                            fit = tryCatch(loess(reads ~ covariate, data = x2s, span = 0.3), error = function(e) NULL)

                            x$reads = NA
                            
                            if (!is.null(fit))
                                if (is.na(fit$s))
                                    {
                                        warning("Using all points since initial loess failed")
                                        fit = loess(reads ~ covariate, data = x2[select, ], span = 1)
                                    }


                            tryCatch(
                                {
                                    if (!is.na(fit$s))
                                        {
                                            domain = domains[[f]]

                                            yrange <- quantile(x2s$reads, prob = c(routlier, 1 - routlier), na.rm = TRUE)                                            
                                            df = data.frame(covariate = seq(domain[1], domain[2], 0.001))
       
                                            if (!is.null(imageroot))
                                                {
                                                    out.png = paste(imageroot, ifelse(grepl("/$", imageroot), '', '.'), f,'_correction.png', sep = '')
                                                    if (verbose) {
                                                        cat("Dumping figure to", out.png, "\n")
                                                    }

                                                    png(out.png, height = 1000, width = 1000) 
                                                    plot(x2s$covariate, x2s$reads, col = alpha('black', 0.1), pch = 19, cex = 0.4, xlim = domain, ylim = yrange, ylab = sprintf('signal before %s correction', f), xlab = f);
                                                    lines(df$covariate, predict(fit, df), col = 'red', lwd = 2)
                                                    dev.off()
                                                }                               
                                            x$reads = x2$reads/predict(fit, x2) ## apply correction
                                        }
                                    else
                                        print("Loess failed, yielding NA loess object, continuing without transforming data")
                                }, error = function(e) print("Unspecified loess or figure output error"))
                        }
                    return(x$reads)
                }

        if (verbose)
            cat('Correcting coverage at individual scales\n')
        
        ## level 1,2, ..., k corrections
        ## these are the computed corrected values that will be input into the objective function
        
        correction = NULL
        for (i in rev(1:length(grs)))
                    {
                        cat('Correcting coverage at ', WID*base^(i-1), 'bp scale, with', nrow(grs[[i]]), 'intervals\n')
                        if (i != length(grs))                            
                            grs[[i]]$reads = grs[[i]]$reads/correction[parentmap[[i]][grs[[i]]$lab, parent], cor]
                        
                        if (WID*base^(i-1) > 1e5) ## for very large intervals do not quantile trim, only remove NA
                            grs[[i]]$reads.corrected = FUN(as.data.frame(grs[[i]]), fields, doutlier = 0, seg = seg)
                        else
                            grs[[i]]$reads.corrected = FUN(as.data.frame(grs[[i]]), fields, seg = seg);  
                        
                        
                        if (is.null(correction))
                            correction = data.table(lab = grs[[i]]$lab, cor = grs[[i]]$reads / grs[[i]]$reads.corrected, key = 'lab')
                        else
                            {
                                ## multiply new correction and old correction
                                old.cor = correction[parentmap[[i]][grs[[i]]$lab, parent], cor]
                                new.cor = grs[[i]]$reads / grs[[i]]$reads.corrected                                     
                                correction = data.table(lab = grs[[i]]$lab,  cor = old.cor * new.cor, key = 'lab') ## relabel with new nodes
                            }                        
                    }
        
        cov.dt$reads.corrected = grs[[1]][cov.dt$lev0, ]$reads.corrected
        cov.dt[reads.corrected < 0, reads.corrected := NA]
        rm(grs)
        gc()
        
        if (verbose)
            cat('Converting to GRanges\n')
        
        gc()
        
        out = seg2gr(as.data.frame(cov.dt), seqlengths = seqlengths(cov))
        
        if (verbose)
            cat('Made GRanges\n')
        
        gc()
        return(out)
                                        #                cov.dt[, reads.corrected := grs[[1]][lev0, reads.corrected]]
    }

#'
#' Takes provided colors and gives them the specified alpha (ie transparency) value
#' @param col RGB color
#' @param alpha value of alpha
#' @name alpha
#' @rdname trackData-class
#' @export
alpha = function(col, alpha)
  {    
    col.rgb = col2rgb(col)
    out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
    names(out) = names(col)
    return(out)
  }

out.raw.rds = paste(opt$outdir, '/cov.raw.rds', sep = '')
out.raw.bedgraph = paste(opt$outdir, 'cov.bedgraph', sep = '/')
out.rds = paste(opt$outdir, '/cov.rds', sep = '')


if (is.null(opt$bam))
    opt$bam = ''
    
if (is.null(opt$cov))
    {
        opt$cov = out.raw.rds
} else if (!file.exists(opt$cov))
      {
          opt$cov = out.raw.rds
      }

opt$midpoint = grepl("(T)|(TRUE)", opt$midpoint, ignore.case = T)

if (file.exists(opt$bam) & !file.exists(opt$cov))
    {
        if (!opt$midpoint)
            cat("Running without midpoint!!!\n")

        print('Doing it!')

        if (is.null(opt$paired))
            opt$paired = TRUE

        if (opt$paired)
          cov = bam.cov.tile(opt$bam, window = opt$window, chunksize = 1e6, midpoint = opt$midpoint, min.mapq = opt$minmapq)  ## counts midpoints of fragments
        else
            {
                sl = seqlengths(BamFile(opt$bam))
                tiles = gr.tile(sl, opt$window)
                cov = bam.cov.gr(opt$bam, tiles, isPaired = NA, isProperPair = NA, hasUnmappedMate = NA)  ## counts midpoints of fragments
                cov$count = cov$records/width(cov)
            }

        saveRDS(cov, out.raw.rds)
        cat('Finished saving coverage RDS\n')
#        bedgraph = data.frame(chrom = as.character(seqnames(cov)), chromStart = start(cov), chromEnd = end(cov), score = cov$count)
#        write.table(bedgraph, out.raw.bedgraph, sep = '\t', col.names = F, row.names = F, quote = F)
#        cat('Finished saving coverage bedgraph\n')
        
    } else if (file.exists(opt$cov))
          {
              cov = readRDS(opt$cov)
     } else
     stop("Can't locate either bam or coverage input file")

gc()

cat('Finished acquiring coverage\n')

correctcov_stub(out.rds, cov)

cat('done\n')
