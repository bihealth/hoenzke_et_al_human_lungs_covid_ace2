#!/bin/bash

sample=$1
ncells=$2
material=$3

mkdir -p data/${material}/cellbender
cd data/${material}/cellbender
cellbender remove-background --input ../cellranger/${sample}_raw.h5 --output ${sample}.h5 --expected-cells ${ncells} --total-droplets-included 20000 --fpr 0.01 --epochs 150 --cuda
