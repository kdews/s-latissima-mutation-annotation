#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH -J validate_variants
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

# Subset regions of interest from genome with BED file


