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
## Run fusion and event calling 

set -euo pipefail
source /etc/profile.d/modules.sh
module purge
module load jabba/1.1

slurm_name_to_ids() {
        SQOUT=squeue_${RANDOM}_$$
        squeue --noheader --format "%i %j" > $SQOUT
        grep $1 $SQOUT | cut -d ' ' -f1 | tr '\n' ':' | sed -e 's/.$//'
        rm -f $SQOUT
}



####################
## Memory presets ##
####################

EVENT_MEM='8G'
FUSION_MEM='8G'
MERGE_MEM='1G'



#######################
## Argument checking ##
#######################
USAGE="Usage: $0 config"

if [ $# -ne 1 ]; then
    echo ""
    echo $USAGE >&2
    echo ""
    exit 1
fi

config=$1


## Get source directory, check config
srcdir=$(realpath $(dirname $0))
bash $srcdir/check-config.sh $config


## Read config 
output_dir=$(grep 'output_directory' $config | awk '{print $2}')
tn_file=$(grep 'tumor_normal_pairs' $config | awk '{print $2}')
gencode=$(grep 'gencode' $config | awk '{print $2}')

## Set dirs
jabbadir=$output_dir/jabba
logdir=$output_dir/logs


dependency="afterok"
while read line; do
  [[ $line =~ ^#.* ]] && continue

  tumor=$(echo $line | awk '{print $1}')
  normal=$(echo $line | awk '{print $2}')

  ## Inputs
  jabba_gg=$jabbadir/$tumor--$normal/jabba.simple.gg.rds
  jabba_cnv_vcf=$jabbadir/$tumor--$normal/jabba.simple.cnv.vcf

  ## Outputs
  jabba_calls=$jabbadir/$tumor--$normal/jabba.events.rds
  jabba_fusions=$jabbadir/$tumor--$normal/jabba.fusions.rds
  jabba_cnv_bed=$jabbadir/$tumor--$normal/jabba.simple.cnv.bed


  if [[ ! -e $jabba_gg ]]; then
    echo "WARNING: Couldn't find JaBbA simplified genome graph: $jabba_gg" >&2
  fi


  ## Call events
  jobname=ggnome.eventcalling.$tumor--$normal
  sbatch_call="Rscript $srcdir/calling/call-events.r \
                --jabba_gg=$jabba_gg \
                --out_file=$jabba_calls"

  sbatch \
    -o $logdir/$jobname.o%A \
    --wckey compbio \
    -n 1 \
    --mem=$EVENT_MEM \
    --job-name=$jobname \
    --wrap="$sbatch_call"

  dependency+=":$(slurm_name_to_ids $jobname)"

  ## Call fusions
  ## Requires gencode GTF
  if [[ ! -z $gencode ]]; then
    jobname=ggnome.fusioncalling.$tumor--$normal
    sbatch_call="Rscript $srcdir/calling/call-fusions.r \
                  --jabba_gg=$jabba_gg \
                  --gencode_gtf=$gencode \
                  --out_file=$jabba_fusions"

    sbatch \
      -o $logdir/$jobname.o%A \
      --wckey compbio \
      -n 1 \
      --mem=$EVENT_MEM \
      --job-name=$jobname \
      --wrap="$sbatch_call"
  fi

done < $tn_file



## Merge calling summaries
jobname=merge_calling_summaries
sbatch_call="Rscript $srcdir/calling/merge-events.r \
              --jabba_dir=$jabbadir \
              --tn_file=$tn_file \
              --out_dir=$output_dir"

sbatch \
  -o $logdir/$jobname.o%A \
  --wckey compbio \
  -n 1 \
  --mem=$EVENT_MEM \
  --dependency=$dependency \
  --job-name=$jobname \
  --wrap="$sbatch_call"

