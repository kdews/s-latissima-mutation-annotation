#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=25gb
#SBATCH --time=2-0
#SBATCH -J vcf_extract
#SBATCH -o %x.log

# Source configuration file
[[ $1 ]] && config_file=$1 || config_file=mut_annot.config
if [[ $config_file ]] && [[ -f $config_file ]]
then
	source $config_file
else
	echo "Error - please provide config file. $config_file not found."
	exit 1
fi

# Optional: Anaconda configuration
[[ $conda_sh ]] && source_conda $conda_sh

# Define BED file
bed=extract_from_vcf.bed

# Create BED file to use with VCFtools (if needed)
if [[ -f $bed ]] && (( $(cat $bed | wc -l) > 1 ))
then
	echo "Using BED file ${bed}..."
else
	if [[ -f $gene_list ]]
	then
		echo "Creating BED file $bed from ${gene_list}..."
		awk '{print $4,$5,$6}' $gene_list > $bed
	else
		echo "Error - no BED file ($bed) detected, and no gene list \
($gene_list) to create it from. Exiting... "
		exit 1
	fi
fi

# Subset regions of interest from annotated VCF into new VCF file
vcftools --vcf snpEff/${vcf_base}.ann.vcf --bed $bed \
--recode --recode-INFO-all \
--out ${vcf_base}.ann.gene_list

# --max-maf 0.1 if VCF is massive
