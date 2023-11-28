#!/bin/bash
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
set -euo pipefail
source /etc/profile.d/modules.sh
module purge
module load R/4.0.0 bedtools/2.29.0



##################
# Memory presets #
##################

GC_MEM='8G'
MAP_MEM='16G'



## Collect arguments
if [ $# -ne 6 ]; then
    echo ""
    echo "Usage: $0 GENOME MAPPABILITY WORKDIR BINSIZE SRCDIR LOGDIR" >&2
    echo ""
    exit 1
fi

GENOME=$1
MAPPABILITY=$2
WORKDIR=$3
BINSIZE=$4
SRCDIR=$5
LOGDIR=$6

mkdir -p $WORKDIR/gccontent $WORKDIR/mappability

echo 'Generating fragCounter resource files...'

## Generate chromosome intervals to tile over
UNPARSED_CHR_FILE=$GENOME.fai
PARSED_CHR_FILE=$WORKDIR/chrom_lengths.txt
awk -v OFS='\t' '{print $1,0,$2}' $UNPARSED_CHR_FILE | egrep '^chr([0-9]+|X|Y)\s'  > $PARSED_CHR_FILE || true

if [[ $(cat $PARSED_CHR_FILE | wc -l) == 0 ]]; then
  awk -v OFS='\t' '{print $1,0,$2}' $UNPARSED_CHR_FILE | egrep '^([0-9]+|X|Y)\s'  > $PARSED_CHR_FILE || true
fi



######################
##    GC Content    ##
######################

UNPARSED_GC_FILE=$WORKDIR/gccontent/grch38.gc.${BINSIZE}.bedtools_nuc.tab.gz
GC_OUTFILE=$WORKDIR/gccontent/gc$BINSIZE.rds

## Generate windows (bedtools makewindows) and compute gc content (bedtools nuc)
## This is quick: ~90s on vm
bedtools nuc -fi $GENOME -bed <(bedtools makewindows -w $BINSIZE -b $PARSED_CHR_FILE) | \
tail -n +2 | awk -v OFS='\t' '{print $1,$2,$3,$5}' | gzip > $UNPARSED_GC_FILE

## Convert to GRanges object
jobname=generate_gc${BINSIZE}_file
cmd="Rscript $SRCDIR/fragcounter/init-bin-scores.r \
          --score_file=$UNPARSED_GC_FILE \
          --chr_intervals=$PARSED_CHR_FILE \
          --window_size=$BINSIZE \
          --out_file=$GC_OUTFILE"

echo 'Generating GC content file...'
sbatch \
  -o $LOGDIR/$jobname.o%A \
  --wckey compbio \
  -n 1 \
  --mem=$GC_MEM \
  --job-name=$jobname \
  --wait \
  --wrap="$cmd"
    



######################
##    Mappability   ##
######################

MAP_OUTFILE=$WORKDIR/mappability/map$BINSIZE.rds

## Bin mappability and convert to GRanges
jobname=generate_map${BINSIZE}_file
cmd="Rscript $SRCDIR/fragcounter/init-bin-scores.r \
          --score_file=$MAPPABILITY \
          --chr_intervals=$PARSED_CHR_FILE \
          --window_size=$BINSIZE \
          --out_file=$MAP_OUTFILE"

echo 'Generating mappability content file...'
sbatch \
  -o $LOGDIR/$jobname.o%A \
  --wckey compbio \
  -n 1 \
  --mem=$MAP_MEM \
  --job-name=$jobname \
  --wait \
  --wrap="$cmd"
