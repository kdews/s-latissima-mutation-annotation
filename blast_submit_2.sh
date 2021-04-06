#!/bin/bash
#SBATCH --mem=20gb
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=32

blastn -num_threads 32 -query blast/${1}_R2.bad.fa -db nt -perc_identity 95 -out blast/${1}_R2.blast.out -outfmt "6 qseqid sseqid pident length"

