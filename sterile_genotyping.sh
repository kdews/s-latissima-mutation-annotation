#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH -J sterile_genotyping
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

# Test R installation
R --version
[[ $? -ne 0 ]] && { echo "Error - check R installation. Exiting..."; exit 1; }

# Create R-friendly version of gene list before running Rscript
if [[ -f $gene_list ]] && [[ -f $R_gene_list ]]
then
	echo "Found $R_gene_list"
elif [[ -f $gene_list ]]
then
	sed "s/#//g" $gene_list > $R_gene_list
else
	echo "Error - $gene_list not found. Exiting..."
	exit 1
fi
# Run Rscript
echo "Running Rscript on $simple_summ and ${R_gene_list}..."
Rscript --vanilla ${scripts_dir}sterile_genotyping.R $simple_summ $R_gene_list \
$annot_summ $validate_bed $indiv_file

[[ $? -eq 0 ]] && echo "Job completed. Find results in:
$annot_summ
$validate_bed
$indiv_file" || \
echo "Error running Rscript. (Error code $?)."

