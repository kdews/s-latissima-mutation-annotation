#!/bin/bash
#SBATCH --partition=cegs
#SBATCH --time=00:30:00

reformat.sh in=bbduk/${1}_R1.bad.fastq.gz in2=bbduk/${1}_R2.bad.fastq.gz out=blast/${1}_R1.bad.fa out2=blast/${1}_R2.bad.fa

