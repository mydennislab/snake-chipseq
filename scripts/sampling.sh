#!/bin/bash

# usage: bash sampling.sh input output_prefix

samtools view $1 | sed 's/ZW:f:/ZW:f:\t/g' | awk -F"\t" -v OFS="\t" -v seed=$RANDOM 'BEGIN {srand(seed)}; {$(NF+1)=int(rand()*100000000)}1' | sort -k1,1 -k15,15rg -k16,16rg | sort -k1,1 -u | cut -f1-13 > $2.sam
samtools view -H $1 > $2.header
cat $2.header $2.sam | samtools view -Sb > $2.bam
