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

## Assemble nanopore reads


########################
## Collect QC metrics ##
########################

## Aggregate in tabular format 
set +e 
while read sample uuid; do 
    output_path=gs://nygc-comp-s-1856-output/ReadsAlignment/$uuid
    stats_json=$output_path/call-lraStats/$sample.lra_stats_aggregate.json

    ## Could be under attempt2, 3, etc.. 
    gsutil -q stat $stats_json 
    exists=$?
    if [[ $exists != 0 ]]; then 
        for i in $(seq 2 4); do

            stats_json=$output_path/call-lraStats/attempt-$i/$sample.lra_stats_aggregate.json
            gsutil -q stat $stats_json 
            exists=$?
            [[ $exists == 0 ]] && break 

        done 
    fi 

    gcloud storage cat $stats_json | \
    jq -r '. as $in | 
        reduce leaf_paths as $path ({};
        . + { ($path | map(tostring) | join(".")): $in | getpath($path) }) | 
        to_entries[] | [.key,.value] | @tsv' | \
    xargs -I {} echo -e "$sample\t$uuid\t{}"

done < $ALIGN_UUIDS > $METADATA_DIR/lra-stats-aggregated.pcr_free.txt
set -e 



#############################################################
## Extract reads overlapping amps of interest and assemble ##
#############################################################

PS1=
. /gpfs/commons/groups/nygcfaculty/kancero/envs/miniconda3/etc/profile.d/conda.sh
conda activate flye 



## For each ONT sample 
while read ont_sample ilmn_dir ilmn_tumor ilmn_normal; do 

    cram=$WORKING_DIR/cram/pcr-free/$ont_sample.cram
    jba_fp=$ilmn_dir/jabba/$ilmn_tumor--$ilmn_normal/jabba.events.genes.tab


    ## Extract relevant footprints for this sample
    amp_fp_list=$(Rscript $SRCDIR/init-export-amp-footprints.r \
                    --in_file=$jba_fp \
                    --sample=$ilmn_tumor \
                    --padding=0 \
                    --out_dir=$WORKING_DIR/footprints) 


    while read amp_bed; do 

        amp_id=$(basename $amp_bed .bed)
        
        
        ## Create output dirs
        flye_dir=$WORKING_DIR/flye/$amp_id.pcr_free
        minimap_dir=$flye_dir/minimap2

        mkdir -p $minimap_dir $canu_dir $miniasm_dir
        

        ## Extract reads
        fastq=$flye_dir/$amp_id.fastq.gz
        samtools view -h -T $REFERENCE -@ 16 -L $amp_bed $cram | samtools fastq | gzip -c > $fastq


        ## CCND1 ecDNAs seem to have bad alignment behavior at ALT contig chr11_KI270827v1_alt
        ## Pull down reads mapping to the alt contig if we need to and use them in assembly
        if [[ $(awk '$1 == "chr11"' $amp_bed | wc -l) > 0 ]]; then 

            ALT_CONTIG="chr11_KI270827v1_alt"
            samtools view -h -T $REFERENCE -@ 16 $cram $ALT_CONTIG | samtools fastq | gzip -c >> $fastq

        fi


        ## Assemble
        flye \
            --nano-hq $fastq \
            --read-error 0.03 \
            --meta \
            --threads 32 \
            --iterations 1 \
            --out-dir $flye_dir


        ## Align assembly to reference 
        minimap2
            -t 9 \
            -ax asm5 \
            $REFERENCE \
            $flye_dir/assembly.fasta | \
        samtools sort - > $minimap_dir/$amp_id-assembly.bam && \
        samtools index $minimap_dir/$amp_id-assembly.bam

    done < $amp_fp_list

done < <(tail -n +2 $ONT_ILMN_MAP)
