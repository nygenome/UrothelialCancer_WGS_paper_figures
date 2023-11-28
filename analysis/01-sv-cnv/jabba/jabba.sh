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
##  Run JaBbA

set -euo pipefail
source /etc/profile.d/modules.sh
module load jabba/1.1

slurm_name_to_ids() {
        SQOUT=squeue_${RANDOM}_$$
        squeue --noheader --format "%i %j" > $SQOUT
        grep $1 $SQOUT | cut -d ' ' -f1 | tr '\n' ':' | sed -e 's/.$//'
        rm -f $SQOUT
}



##########################################################
## Memory presets                                       ##
## This should be fine but may have to be increased for ##
## highly rearranged samples                            ##
##########################################################

JABBA_MEM='60G'



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
project_dir=$(awk '($1 =="project_directory") {print $2}' $config)
output_dir=$(awk '($1 =="output_directory") {print $2}' $config)
junc_dir_conf=$(awk '($1 =="junction_directory") {print $2}' $config)
junc_suffix=$(awk '($1 =="junction_suffix") {print $2}' $config)
seg_dir=$(awk '($1 =="segmentation_directory") {print $2}' $config)
seg_suffix=$(awk '($1 =="segmentation_suffix") {print $2}' $config)
tn_file=$(awk '($1 =="tumor_normal_pairs") {print $2}' $config)
blacklist=$(awk '($1 =="blacklist") {print $2}' $config)
gene_bed=$(awk '($1 =="gene_bed") {print $2}' $config)
purity_ploidy=$(awk '($1 =="purity_ploidy") {print $2}' $config)
coverage=$(awk '($1 =="coverage") {print $2}' $config)
cores=$(awk '($1 =="cores") {print $2}' $config)
use_tiers=$(awk '($1 =="use_tiers") {print $2}' $config)
additional_args=$(awk -v FS=':' '($1=="additional_args") {print $2}' $config)

## Make top-level directories
logdir=$output_dir/logs
hetdir=$output_dir/hets
jabbadir=$output_dir/jabba

mkdir -p $logdir $hetdir $jabbadir



#####################################
## Make male and female nseg files ##
#####################################

echo "Generating --nseg files..."

chrlen=$output_dir/chrom_lengths.txt
Rscript $srcdir/utils/init-nseg.r --chrlen=$chrlen

module unload R



###############
## Run jabba ##
###############

if [[ $use_tiers == 'true' ]]; then 
  echo 'Marking high-confidence junctions as tier 2, and everything else as tier 3...'
fi

while read line; do 
  [[ $line =~ ^#.* ]] && continue

  tumor=$(echo $line | awk '{print $1}')
  normal=$(echo $line | awk '{print $2}')
  gender=$(echo $line | awk '{print $3}')
  
  ## Make sure a gender was supplied
  if [[ -z $gender || ! "$gender" =~ ^(male|female) ]]; then
    echo "ERROR: Third column in $tn_file file must be of value: male|female. Value given for $tumor was $gender. Defaulting to MALE"
    gender='male'
  fi

  
  ## Set input files, output directory
  jabba_outdir=$jabbadir/$tumor--$normal
  nseg=$output_dir/${gender}_nseg.rds
  
  mkdir -p $jabba_outdir
  

  ## Coverage data
  if [[ $coverage == 'ratio' ]]; then
    cov_file=$output_dir/fragcounter/${tumor}--${normal}_coverage_ratio.rds
    field='ratio'
  else
    cov_file=$output_dir/dryclean/tumor-decomposition/$tumor.dryclean-decomp.rds
    field='foreground'
  fi
  
  
  ## Junction files
  if [[ ! -z $junc_dir_conf ]]; then

    junctions=$junc_dir_conf/$tumor--$normal$junc_suffix

  else 

    junc_dir=$project_dir/Somatic_analysis/$tumor--$normal
    junc_file=$tumor--$normal.sv.annotated.v6.somatic.high_confidence.final.bedpe

    ## Junctions might be in another location
    if [[ ! -e $junc_dir/$junc_file ]]; then
      junc_dir=$project_dir/compbio/Somatic_analysis/$tumor--$normal
    fi

    junctions=$junc_dir/$junc_file

  fi
  
  
  ## Quit if we can't find junctions or coverage
  if [[ ! -e $junctions ]]; then
    echo "WARNING: Couldn't find junction file: $junctions, skipping..." >&2
    continue
  fi
  
  if [[ ! -e $cov_file ]]; then
    echo "WARNING: Couldn't find coverage file: $cov_file, skipping..." >&2
    continue
  fi
  
  
  ## If we're using tiered junctions, encode in bedpe score column
  if [[ $use_tiers == 'true' ]]; then
    
    echo "TODO: reimplement tier system"
    exit 0

  else

    tierstring=''

  fi
  
  
  ## Exclude chrY junctions
  junctions_yfiltered=$jabba_outdir/$tumor--$normal.junctions.yfiltered.bedpe
  awk '$1 != "chrY" && $4 != "chrY" && $1 != "Y" && $4 != "Y"' $junctions > $junctions_yfiltered
  junctions=$junctions_yfiltered
  
  
  ## Try to find heterozygous germline variants
  germline_file_unparsed=$project_dir/Sample_$normal/compbio/variants/tmp/$tumor--$normal.haplotypecaller.v3.5.0.alleles.txt
  germline_file_parsed=$hetdir/$tumor--$normal.haplotypecaller.v3.5.0.alleles.txt
  
  if [[ -e $germline_file_unparsed ]]; then
    echo -e 'seqnames\tstart\tend\tref.count.t\talt.count.t\talt.frac.t\tref.count.n\talt.count.n' > $germline_file_parsed

    tail $germline_file_unparsed -n +2 | awk -F'\t' -v OFS='\t' '{if ($8+$9 == 0) d=1; else d=$8+$9}''{print $1,$2,$2,$8,$9,$9/d,$19,$20}' \
    | egrep 'chr([0-9]+|X|Y)\s' >> $germline_file_parsed
  fi

  if [[ -e $germline_file_parsed ]]; then
    hetstring="--hets=$germline_file_parsed"
  else
    hetstring=""
  fi
  
  
  
  ## If method is sequenza we need to grab heterozygous germline variants
  if [[ $purity_ploidy == 'sequenza' ]]; then
    
    if [[ -e $germline_file_parsed ]]; then
      ppstring="--ppmethod=sequenza"
    else
      echo "WARNING: Couldn't fine germline het file: $germline_file_unparsed"
      echo "WARNING: Using ppgrid estimation for $tumor--$normal"
      ppstring="--ppmethod=ppgrid"
    fi
  
  ## If ppgrid, just specify as such
  elif [[ $purity_ploidy == 'ppgrid' ]]; then
    ppstring="--ppmethod=ppgrid"
  
  ## Extract manual estimates from the csv 
  else
      
      ## Make sure we have estimates for this sample
      if [[ $(grep $tumor $purity_ploidy | wc -l) == 0 ]]; then
        echo "Couldn't find purity/ploidy estimate for $tumor, skipping..."
        continue
      fi
  
      purity=$(grep $tumor $purity_ploidy | cut -d',' -f2)
      ploidy=$(grep $tumor $purity_ploidy | cut -d',' -f3)
      ppstring="--purity=$purity --ploidy=$ploidy"
  fi
  
  
  
  ## Coverage blacklist is optional
  if [[ -z $blacklist ]]; then
    blstring=''
  else
    blstring="--blacklist.coverage=$blacklist"
  fi
  

  ## Segmentation is optional
  if [[ -z $seg_dir ]]; then
    segstring=''
    echo "no segmentation found!"
  else 
    segstring="--seg=$seg_dir/$tumor--$normal$seg_suffix"
    blstring=''
  fi

  ## Run jabba
  jobname=jabba.$tumor--$normal
  sbatch_call="jba --cores=$cores \
              --outdir=$jabba_outdir \
              --field=$field \
              --nseg=$nseg \
              --lp=TRUE \
              --ism=TRUE \
              $blstring \
              $ppstring \
              $hetstring \
              $tierstring \
              $segstring \
              $additional_args \
              $junctions \
              $cov_file"
  
  echo $sbatch_call > $logdir/$jobname.$(date +"%m%d%Y_%H%S")
            
  sbatch \
    -o $logdir/$jobname.o%A \
    --wckey compbio \
    -n $cores \
    --mem=$JABBA_MEM \
    --chdir=$jabba_outdir \
    --job-name=$jobname \
    --wrap="$sbatch_call"
  jobid=$(slurm_name_to_ids $jobname)
  
  
  ## Convert jabba CNV VCF to BED 
  jobname=vcf2bed.$tumor--$normal
  jobdependency="afterok:$jobid"

  sbatch_call="Rscript $srcdir/utils/cnv-vcf-to-bed.r \
                --in_file=$jabba_outdir/jabba.cnv.vcf \
                --gene_bed=$gene_bed \
                --out_file=$jabba_outdir/jabba.cnv.bed"

  sbatch \
    -o $logdir/$jobname.o%A \
    --wckey compbio \
    -n 1 \
    --mem=1G \
    --dependency=$jobdependency \
    --kill-on-invalid-dep=yes \
    --job-name=$jobname \
    --wrap="$sbatch_call"
  
done < $tn_file
