#!/bin/bash
#SBATCH --mem=20gb
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=12

blastn -num_threads 12 -query blast/${1}_R2.bad.fa -db nt -perc_identity 95 -out blast/${1}_R2.blast.out -outfmt "6 qseqid sseqid pident length"

