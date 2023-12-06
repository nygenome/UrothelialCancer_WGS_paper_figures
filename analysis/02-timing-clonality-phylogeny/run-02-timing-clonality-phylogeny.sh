#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper & Timothy R. Chu 

################################################################# /COPYRIGHT ###
################################################################################

## Use "pileup method" for joint calling, use that to get mutation timing 
## and phylogenetic trees 



###################
## Joint calling ##
###################

## Merge calls 
while read patient; do 

    union_vcf=$VARIANT_DIR/$patient.merged.highconf.vcf

    bcftools merge --force-samples $PROJECT_DIR/compbio/Somatic_analysis/$patient*/*.snv.indel.final.v6.annotated.ffpePON_filtered.vcf | \
    bcftools filter -i INFO/HighConfidence=1 -o $union_vcf

done < (cut -f4 $TNFILE | sort -u)


## For each tumor, select calls with evidence in pileup
while read tumor normal gender patient; do 

    tn=$tumor--$normal
    union_vcf=$VARIANT_DIR/$patient.merged.highconf.vcf
    vcf=$PROJECT_DIR/compbio/Somatic_analysis/$tn/$tn.snv.indel.final.v6.annotated.ffpePON_filtered.vcf

    python $SRCDIR/get_vaf.py \
        -t $tumor.final.bam \
        -n $normal.final.bam \
        --tumor $tumor \
        --normal $normal \
        -v $vcf \
        -u $union_vcf \
        -o $VARIANT_DIR/$patient/$tn.highconf.ffpefilter.union.vcf

done < $TNFILE



#####################
## Mutation timing ##
#####################

while read tumor normal gender patient; do

    tn=$tumor--$normal

    ## Inputs
    vcf=$VARIANT_DIR/$patient/$tn.highconf.ffpefilter.union.vcf
    jba=$JABBA_DIR/jabba/$tn/jabba.simple.rds
    purity=$(awk -v FS=',' -v tumor=$tumor '($1 == tumor) {print $2}' $PURITY_PLOIDY)

    ## Case where this patient only has a single sample
    if [[ ! -f $vcf ]]; then
        vcf=$PROJECT_DIR/compbio/Somatic_analysis/$tn/$tn.snv.indel.final.v6.annotated.ffpePON_filtered.vcf
    fi

    [[ ! -f $jba ]] && continue

    ## Prepped input files
    ascn=$WORKING_DIR/cn/$tn-jabba-ascn.rds

    ## Output files
    out_file_vcf=$WORKING_DIR/mutationtimer/$tn-mutationtimer-snv.vcf
    out_file_bed=$WORKING_DIR/mutationtimer/$tn-mutationtimer-cnv.bed
    out_file_rds=$WORKING_DIR/mutationtimer/$tn-mutationtimer.rds
    out_file_txt=$WORKING_DIR/mutationtimer/$tn-mutationtimer.summary.txt
    out_file_png=$FIG_DIR/$tn-mutationtimer.png

    Rscript $SRCDIR/init-mutation-timer-input.r \
    --jba=$jba \
    --purity=$purity \
    --tumor=$tumor \
    --gender=$gender \
    --out_file_cn=$ascn

    Rscript $SRCDIR/run-mutation-timer.r \
    --cn=$ascn \
    --purity=$purity \
    --tumor=$tumor \
    --gender=$gender \
    --vcf=$vcf \
    --out_file=$out_file_rds

    Rscript $SRCDIR/postprocess-mutation-timer.r \
    --cn=$ascn \
    --mt=$out_file_rds \
    --tumor=$tumor \
    --purity=$purity \
    --gender=$gender \
    --vcf=$vcf \
    --out_file_vcf=$out_file_vcf \
    --out_file_bed=$out_file_bed \
    --out_file_summary=$out_file_txt \
    --out_file_pdf=$out_file_png

done < $TNFILE



##############################
## Build phylogenetic trees ##
##############################

while read patient; do 

    ## Preprocessing 
    bcftools merge *.vcf.gz | bcftools query -H -f '%CHROM\t%POS\t%REF:%ALT:%INFO/CSQ\t[%CCF\t]\n' > $patient.merge.ccf.txt
    awk -F $'\t' 'BEGIN {OFS = FS} $3 = $3 FS "0.0"' PM63.merge.ccf.txt | sed 's/\.\t/0.0\t/g' | sed '1s/0.0/Normal/' | sed '1s/\[.\{1,2\}\]//g'> $patient.lichee.input.txt

    ## Run lichee 
    $LICHEE_REPO/lichee.sh -build \
        -i $patient.lichee.input.txt \
        -dot \
        -cp \
        --minVAFPresent 0.05 \
        --maxVAFAbsent 0.0 \
        -n 0  \
        -o $output_lichee

done < (cut -f4 $TNFILE | sort -u)