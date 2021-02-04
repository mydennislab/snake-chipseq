#!/bin/bash
# usage: bash sampling.sh input output_prefix
set -euxo

samtools view $1 | sed 's/ZW:f:/ZW:f:\t/g' | python3 /share/dennislab/projects/noncoding/long_bowtie2_PE/scripts/csem_bestPE.py | sed 's/ZW:f:\t/ZW:f:/g' > $2.tmp
samtools view -H $1 > $2.header
cat $2.header <(cut -f1-11 $2.tmp) | samtools view -Sb | samtools sort > $2.bam

