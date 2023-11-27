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
##
## Description: Run preprocessing to generate the coverage files necessary for a Jabba run. 
##
## STEP 1: Run fragcounter, generating coverage in bins <config: binsize>. Note that the minimum binsize allowed by 
##         dryclean in STEP 2 is 1kb 
##
## STEP 2: Generate a dryclean panel of normals using the normals as specified in the tumor-normal pair file. This 
##         step can be skipped if you already have a dryclean PON that reflects the technical/batch biases that you
##         expect to be present in the tumor samples. If you have this, symlink the detergent.rds file in the 
##         dryclean output directory: <config: output_directory>/dryclean. If using many samples for the PON, 
##         you may have to update the memory provided 
## 
## STEP 3: Run dryclean on tumor samples using the drycleaned PON
## 
set -euo pipefail
source /etc/profile.d/modules.sh
module load R/4.0.0


slurm_name_to_ids() {
        SQOUT=squeue_${RANDOM}_$$
        squeue --noheader --format "%i %j" > $SQOUT
        grep $1 $SQOUT | cut -d ' ' -f1 | tr '\n' ':' | sed -e 's/.$//'
        rm -f $SQOUT
}



##########################################################
## Memory presets                                       ##
## DRYCLEAN_PON_MEM may have to be increased if working ##
## with a very large number of samples                  ##
##########################################################

FRAGCOUNTER_MEM='45G'
RATIO_MEM='4G'
DRYCLEAN_PON_MEM='150G'
DRYCLEAN_MEM='24G'



#######################
## Argument checking ##
#######################
USAGE="Usage: $0 <fragcounter|dryclean-pon|dryclean-germline|dryclean> config"

if [ $# -ne 2 ]; then
    echo ""
    echo $USAGE >&2
    echo ""
    exit 1
fi

if [[ $1 != 'dryclean' && $1 != 'dryclean-pon' && $1 != 'dryclean-germline' && $1 != 'fragcounter' ]]; then
    echo ""
    echo $USAGE >&2
    echo ""
    exit 1
fi

cmd=$1
config=$2

## Get source directory, check config
srcdir=$(realpath $(dirname $0))
bash $srcdir/check-config.sh $config

## Read config 
project_dir=$(grep 'project_directory' $config | awk '{print $2}')
output_dir=$(grep 'output_directory' $config | awk '{print $2}')
bam_dir=$(awk '($1 =="bam_directory") {print $2}' $config)
bam_suffix=$(awk '($1 =="bam_suffix") {print $2}' $config)
tn_file=$(grep 'tumor_normal_pairs' $config | awk '{print $2}')
genome=$(grep 'genome' $config | awk '{print $2}')
mappability=$(grep 'mappability' $config | awk '{print $2}')
purity_ploidy=$(grep 'purity_ploidy' $config | awk '{print $2}')
cores=$(grep 'cores' $config | awk '{print $2}')
binsize=$(grep 'binsize' $config | awk '{print $2}')

## Make top-level directories
logdir=$output_dir/logs
fragdir=$output_dir/fragcounter
drydir=$output_dir/dryclean

mkdir -p $logdir $fragdir $drydir/tumor-decomposition $drydir/normal-decomposition



#################
## FRAGCOUNTER ##
#################

if [[ $cmd == 'fragcounter' ]]; then
  
  ## Generate resource files
  bash $srcdir/fragcounter/init-gc-map-files.sh $genome $mappability $output_dir $binsize $srcdir $logdir
  
  module load samtools/1.6
  
  while read line; do
    [[ $line =~ ^#.* ]] && continue
  
    tumor=$(echo $line | awk '{print $1}')
    normal=$(echo $line | awk '{print $2}')
    gender=$(echo $line | awk '{print $3}')

    ## Make sure a gender was supplied
    if [[ -z $gender || ! "$gender" =~ ^(male|female) ]]; then
      echo "ERROR: Third column in $tn_file file must be of value: male|female. Value given for $tumor was $gender"
      exit 1
    fi
    
    ## Specify input files 
    if [[ ! -z $bam_dir ]]; then

      tumor_bam=$bam_dir/$tumor$bam_suffix
      normal_bam=$bam_dir/$tumor$bam_suffix

    else
    
      tumor_bam=$project_dir/Sample_$tumor/analysis/$tumor.final.bam
      normal_bam=$project_dir/Sample_$normal/analysis/$normal.final.bam
      
      ## Input bams may be in another location
      if [[ ! -e $tumor_bam ]]; then
        tumor_bam=$project_dir/Sample_$tumor/compbio/analysis/$tumor.final.bam
        normal_bam=$project_dir/Sample_$normal/compbio/analysis/$normal.final.bam
      fi
    fi
    
    ## Or they might not exist
    if [[ ! -e $tumor_bam ]]; then
      echo "WARNING: Could not find tumor BAM: $tumor_bam" >&2
    fi
    
    if [[ ! -e $normal_bam ]]; then
      echo "WARNING: Could not find normal BAM: $normal_bam" >&2
    fi
      
    mkdir -p $fragdir/$tumor $fragdir/$normal
      
      
      
    ## Tumor
    jobname=fragcounter.$tumor
    sbatch_call="Rscript $srcdir/fragcounter/fragcounter.r \
                --bam=$tumor_bam \
                --midpoint=TRUE \
                --window=$binsize \
                --paired=TRUE \
                --gcmapdir=$output_dir \
                --outdir=$fragdir/$tumor"
    
    [[ ! -e $tumor_bam ]] && continue
    
    sbatch \
        -o $logdir/$jobname.o%A \
        --wckey compbio \
        -n 1 \
        --mem=$FRAGCOUNTER_MEM \
        --job-name=$jobname \
        --wrap="$sbatch_call"
      
    jobid_tumor=$(slurm_name_to_ids $jobname)
        
        
        
    ## Normal
    jobname=fragcounter.$normal
    sbatch_call="Rscript $srcdir/fragcounter/fragcounter.r \
                --bam=$normal_bam \
                --midpoint=TRUE \
                --window=$binsize \
                --paired=TRUE \
                --gcmapdir=$output_dir \
                --outdir=$fragdir/$normal"
    
    
    ## Might have already kicked this normal off
    jobid_normal=$(slurm_name_to_ids $jobname)
    
    if [[ -z $jobid_normal ]]; then
    
      sbatch \
          -o $logdir/$jobname.o%A \
          --wckey compbio \
          -n 1 \
          --mem=$FRAGCOUNTER_MEM \
          --job-name=$jobname \
          --wrap="$sbatch_call"
          
      jobid_normal=$(slurm_name_to_ids $jobname)
      
    fi
          
    
    ## TN ratio
    tumor_rds=$fragdir/$tumor/cov.rds
    normal_rds=$fragdir/$normal/cov.rds
    coverage_file=$fragdir/${tumor}--${normal}_coverage_ratio.csv
  
    jobdependency="afterok:$jobid_tumor:$jobid_normal"
    jobname=covratio.$tumor--$normal
    cmd="Rscript $srcdir/fragcounter/compute-tn-depth-ratio.r \
         --tumor=$tumor_rds \
         --normal=$normal_rds \
         --gender=$gender \
         --outfile=$coverage_file"
    
    sbatch \
      -o $logdir/$jobname.o%A \
      --wckey compbio \
      -n 1 \
      --mem=$RATIO_MEM \
      --job-name=$jobname \
      --dependency=$jobdependency \
      --wrap="$cmd"
    
  done < $tn_file
fi



#############################
## DRYCLEAN PON GENERATION ##
#############################

if [[ $cmd == 'dryclean-pon' ]]; then
  
  ## Make list of normal samples and their fragcounter output files
  normal_file=$drydir/normal-samples.tab
  grep -v '^#' $tn_file | awk '{print $2}' | sort -u | xargs -I {} echo -e "{}\t$fragdir/{}/cov.rds" > $normal_file
  
  jobname="dryclean.detergent_generation"
  sbatch_call="Rscript $srcdir/dryclean/init-detergent.r \
          --normals=$normal_file \
          --outdir=$drydir \
          --cores=$cores"

  sbatch \
    -o $logdir/$jobname.o%A \
    --wckey compbio \
    --partition=bigmem \
    -n 4 \
    --mem=$DRYCLEAN_PON_MEM \
    --job-name=$jobname \
    --wrap="$sbatch_call"

fi



if [[ $cmd == 'dryclean-germline' ]]; then
  
  ## Make list of normal samples and their fragcounter output files
  normal_file=$drydir/normal-samples.tab
  grep -v '^#' $tn_file | awk '{print $2}' | sort -u | xargs -I {} echo -e "{}\t$fragdir/{}/cov.rds\t$drydir/normal-decomposition/{}.dryclean-decomp.rds" > $normal_file
  
  
  ## Decompose each normal
  dependency='afterok'
  while read normal; do 
    
    ## Specify input and output files
    cov_file=$fragdir/$normal/cov.rds
    detergent=$drydir/detergent.rds
    out_file=$drydir/normal-decomposition/$normal.dryclean-decomp.rds
    
    jobname="dryclean.decomp.$normal"
    sbatch_call="Rscript $srcdir/dryclean/run-dryclean.r \
            --coverage=$cov_file \
            --detergent=$detergent \
            --germline_mode \
            --out_file=$out_file"

    sbatch \
      -o $logdir/$jobname.o%A \
      --wckey compbio \
      -n 1 \
      --mem=$DRYCLEAN_MEM \
      --job-name=$jobname \
      --wrap="$sbatch_call"
    
    dependency="$dependency:$(slurm_name_to_ids $jobname)"
    
  done < <(grep -v '^#' $tn_file | awk '{print $2}' | sort -u)
  
  ## Identify germline events
  jobname="dryclean.identify-germline"
  sbatch_call="Rscript $srcdir/dryclean/identify-germline.r \
    --normals=$normal_file \
    --out_file=$drydir/germline_markers.rds \
    --cores=1"
  
   sbatch \
     -o $logdir/$jobname.o%A \
     --wckey compbio \
     -n 4 \
     --mem=30G \
     --job-name=$jobname \
     --dependency=$dependency \
     --wrap="$sbatch_call"
  
  
fi



##############
## DRYCLEAN ##
##############

if [[ $cmd == 'dryclean' ]]; then
  
  while read line; do
    [[ $line =~ ^#.* ]] && continue
    
    tumor=$(echo $line | awk '{print $1}')
    
    ## Specify input and output files
    cov_file=$fragdir/$tumor/cov.rds
    detergent=$drydir/detergent.rds
    out_file=$drydir/tumor-decomposition/$tumor.dryclean-decomp.rds
    
    jobname="dryclean.decomp.$tumor"
    sbatch_call="Rscript $srcdir/dryclean/run-dryclean.r \
            --coverage=$cov_file \
            --detergent=$detergent \
            --germline_file=$drydir/germline_markers.rds \
            --out_file=$out_file"

    sbatch \
      -o $logdir/$jobname.o%A \
      --wckey compbio \
      -n 1 \
      --mem=$DRYCLEAN_MEM \
      --job-name=$jobname \
      --wrap="$sbatch_call"
    
  done < $tn_file
  
fi
