#!/bin/bash

if [ $# -ne 5 ]; then
    echo "Usage: $0 POSITIONS_FILE BAM SAMPLE_NAME OUT REF" >&2
    exit 1
fi

positions_file=$1
bam=$2
sample_name=$3
out=$4
ref=$5

#module load samtools

## b38
/nfs/sw/samtools/samtools-1.9/bin/samtools mpileup -f $ref -q 10 -Q 10 --ff UNMAP --ff SECONDARY --ff QCFAIL --ff DUP -l $positions_file -B $bam | sequenza-utils pileup2acgt -p - -o $out
