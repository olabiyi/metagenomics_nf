#!/usr/bin/env bash

# If the fastq files exist locally and are contained in 
# folders where each folder is named with the corresponding sample name

# Get the samples names i.e. the directory names
SAMPLES=($(ls -1 01.raw_data/ |sort -uV ))

# Automatically generate files.csv
(echo -ne "sample,forward,reverse\n"; \
 for sample in ${SAMPLES[*]}; do echo -ne "${sample},01.raw_data/${sample}/${sample}_R1.fastq.gz,01.raw_data/${sample}/${sample}_R2.fastq.gz\n"; done) > files.csv
