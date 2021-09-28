#!/bin/bash
###usage: bash ./run_delly_Macaca_mulatta.sh sample
#####Macaca_mulatta

echo "script: $0"
echo "first: $1"
cd /data1/nyy/offtarget_addfq/delly_SV
sample=$1
genome=/data/dsx/Genome/bwaindex/Macaca_mulatta.fa
bam=/data1/nyy/offtarget/bamfile/${sample}_chrok.bam   #08431WT.sorted_dupm.bam

delly call -t INS -o ${sample}_INSnew.bcf -g ${genome} ${bam} >logINS.${sample} 2>&1

# delly call -t INS -o 08431WT_INSnew.bcf -g /data/dsx/Genome/bwaindex/Macaca_mulatta.fa /data1/nyy/offtarget/bamfile/08431WT.sorted_dupm.bam >logINS.08431WT 2>&1
