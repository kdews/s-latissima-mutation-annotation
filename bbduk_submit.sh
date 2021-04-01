#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=125gb
#SBATCH --cpus-per-task=12

contaminants=/project/noujdine_61/common_resources/all_contaminants_ref_lib.fasta.gz
data_dir=${1}

bbduk.sh -Xmx106g threads=12 rskip=6 prealloc=t ordered=t ref=${contaminants} in=${1}/${2}_R1.fastq.gz in2=${1}/${2}_R2.fastq.gz out=bbduk/${2}_R1.fastq.gz out2=bbduk/${2}_R2.fastq.gz outm=bbduk/${2}_R1.bad.fastq.gz outm2=bbduk/${2}_R2.bad.fastq.gz

