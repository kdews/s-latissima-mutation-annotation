#!/bin/bash
#SBATCH -p cegs
#SBATCH -t 10:00:00
#SBATCH --mem 50gb

# Help message
if [[ $1 = "-h" ]] || [[ $1 = "--help" ]] || [[ -z "$1" ]]; then
  printf "Premade for running Kellys_code_with_indels.py to \
annotate \e[3mSaccharina latissima\e[0m VCFs.\n\
Usage: sbatch Kellys_code_with_indels_submit.sh input.vcf\n\n\
Requires:   biopython\n"
  exit 0
fi


# Path to script
script_path=/home1/kdeweese/scripts/s-latissima-mutation-annotation
# CDS coordinates file
coords=assembly_ST_collapse_with_short_genes.complete_ORFs.cds_coords
assembly=assembly_ST_collapse_with_short_genes.fa
vcf=$1

source activate biopython
python ${script_path}/Kellys_code_with_indels.py $coords $assembly $vcf

