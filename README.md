# UrothelialCancer_WGS_paper_figures
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

This repository contains the code necessary for reproducing the figures and custom analysis in the following manuscript:  
Nguyen, D.D., Hooper, W.F., Liu, W. et al. The interplay of mutagenesis and ecDNA shapes urothelial cancer evolution. Nature 635, 219–228 (2024).   
https://doi.org/10.1038/s41586-024-07955-3

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

* Other utilities:
    * bedtools
    * samtools
    * LICHeE
    * AmpliconArchitect/AmpliconClassifier
    * minimap2
    * flye
 