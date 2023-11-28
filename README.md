# UrothelialCancer_WGS_paper_figures
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

This repository contains the code necessary for reproducing the figures and custom analysis in the following manuscript:  
[The Interplay between Mutagenesis and Extrachromosomal DNA Shapes Urothelial Cancer Evolution](https://www.biorxiv.org/content/10.1101/2023.05.07.538753v1)

## Dependencies
* R 4.0.0 or greater
* The following R packages (*NOTE*: not all of these are required for every component. Please see scripts for individual dependencies)
    * optparse
    * ComplexHeatmap
    * circlize
    * viridisLite
    * data.table
    * ggplot2
    * reshape2
    * ggpubr
    * ggrepel
    * gGnome
    * gUtils
    * rtracklayer
    * VariantAnnotation
    * fishhook
    * MutationTimeR
    * plotly
    * deconstructSigs
    * JaBbA
    * parallel
    * Matrix
    * stringr
    * dryclean
    * magrittr
    * GenomicRanges
    * skitools
    * bamUtils
    * BSgenome.Hsapiens.UCSC.hg38

* The following packages
    * bedtools
    * samtools
 