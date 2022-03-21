#!/bin/bash

# this script will print the number of read mapped to the virral chromosome to stdout and create a coverage file for read density along the viral chromosome
sample=$1 ## sample name for output file
bam=$2 ## path to bam file
chrom=$3 ## which chromosome in index (e.g., NC_045512.2 for SARS-CoV-2)

samtools view -b ${bam} ${chrom} > ${sample}.bam
samtools index ${sample}.bam
echo "${sample}	$(samtools flagstat ${sample}.bam | head -n1 | cut -f 1 -d " ")" 
genomeCoverageBed -ibam ${sample}.bam -dz > ${sample}_cov.txt
