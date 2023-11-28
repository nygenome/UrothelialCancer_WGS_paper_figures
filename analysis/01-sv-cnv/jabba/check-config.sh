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
## Check workflow configuration

USAGE="Usage: $0 CONFIG"

if [ $# -ne 1 ]; then
    echo ""
    echo $USAGE >&2
    echo ""
    exit 1
fi

CONFIG=$1

echo "Checking config..."



###############################################
## Check to make sure all fields are present ##
###############################################

FIELDS=('project_directory' 'output_directory' 'tumor_normal_pairs' 'genome' 'mappability' 'blacklist' 'gencode' 'gene_bed' 'purity_ploidy' 'cores' 'binsize' 'use_tiers')
for f in ${FIELDS[@]}; do
  if [[ $(grep $f $CONFIG | wc -l) == 0 ]]; then
    echo "ERROR: Config missing field $f" >&2
    exit 1
  fi
done



#########################
## Sanity check fields ##
#########################

project_dir=$(awk '($1 =="project_directory") {print $2}' $CONFIG)
output_dir=$(awk '($1 =="output_directory") {print $2}' $CONFIG)
tumor_normal_pairs=$(awk '($1 =="tumor_normal_pairs") {print $2}' $CONFIG)
genome=$(awk '($1 =="genome") {print $2}' $CONFIG)
mappability=$(awk '($1 =="mappability") {print $2}' $CONFIG)
blacklist=$(awk '($1 =="blacklist") {print $2}' $CONFIG)
gencode=$(awk '($1 =="gencode") {print $2}' $CONFIG)
gene_bed=$(awk '($1 =="gene_bed") {print $2}' $CONFIG)
purity_ploidy=$(awk '($1 =="purity_ploidy") {print $2}' $CONFIG)
coverage=$(awk '($1 =="coverage") {print $2}' $CONFIG)
cores=$(awk '($1 =="cores") {print $2}' $CONFIG)
binsize=$(awk '($1 =="binsize") {print $2}' $CONFIG)
use_tiers=$(awk '($1 =="use_tiers") {print $2}' $CONFIG)



## Project and output directories
if [[ ! -d $project_dir ]]; then
  echo "ERROR: project_directory doesn't exist: $project_dir" >&2
  exit 1
fi

if [[ ! -d $output_dir ]]; then
  echo "ERROR: output_directory doesn't exist: $output_dir" >&2
  exit 1
fi


## Resource files
if [[ ! -f $tumor_normal_pairs ]]; then
  echo "ERROR: Can't find tumor_normal_pairs file: $tumor_normal_pairs" >&2
  exit 1
fi

if [[ ! -f $genome ]]; then
  echo "ERROR: Can't find reference genome: $genome" >&2
  exit 1
fi

if [[ ! -f $mappability ]]; then
  echo "ERROR: Can't find mappability file: $mappability" >&2
  exit 1
fi

if [[ -z $blacklist ]]; then
  echo "WARNING: No blacklist provided" >&2
elif [[ ! -f $blacklist ]]; then
  echo "ERROR: Can't find blacklist file: $blacklist" >&2
  exit 1
fi

if [[ -z $gencode ]]; then
  echo "WARNING: No gencode GTF provided. This is required for fusion calling." >&2
elif [[ ! -f $gencode ]]; then
  echo "ERROR: Can't find gencode GTF file: $gencode" >&2
  exit 1
fi

if [[ -z $gene_bed ]]; then
  echo "ERROR: No gene BED file provided." >&2
  exit 1
elif [[ ! -f $gene_bed ]]; then
  echo "ERROR: Can't find gene BED file: $gene_bed" >&2
  exit 1
fi


## Purity ploidy
if [[ ! -f $purity_ploidy && $purity_ploidy != 'ppgrid' && $purity_ploidy != 'sequenza' ]]; then
  echo "ERROR: purity_ploidy should be 'ppgrid', 'sequenza', or a path to a purity/ploidy table: $purity_ploidy"
  exit 1
fi

## Coverage input
if [[ ! -f $coverage && $coverage != 'ratio' && $coverage != 'dryclean' ]]; then
  echo "ERROR: coverage should be 'ratio' or 'dryclean'. Given: $coverage"
  exit 1
fi

## Cores and bin size
num='^[0-9]+$'
if ! [[ $cores =~ $num ]] ; then
   echo "ERROR: cores must be a number: $cores" >&2
   exit 1
fi

if ! [[ $binsize =~ $num ]] ; then
   echo "ERROR: binsize must be a number: $binsize" >&2
   exit 1
fi

## Tiers flag
if [[ $use_tiers != 'true' && $use_tiers != 'false' ]] ; then
   echo "ERROR: use_tiers must be set to true or false: $use_tiers" >&2
   exit 1
fi

echo "Config OK!"
exit 0
