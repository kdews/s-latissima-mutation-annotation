#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=25gb
#SBATCH --time=2-0
#SBATCH -J ann_SnpEff
#SBATCH -o %x.log

# Source configuration file
[[ $1 ]] && config_file=$1 || config_file=mut_annot.config
if [[ $config_file ]] && [[ -f $config_file ]]
then
	source $config_file
else
	echo "Error - please provide config file. \
$config_file not found."
	exit 1
fi

# Optional: Anaconda configuration
[[ $conda_sh ]] && source_conda $conda_sh

# Test Java install (must be >= v1.8)
java --version
[[ $? -eq 0 ]] && echo "Java successfully loaded." || \
{ echo "Error arose while testing Java install. Exiting..."; exit 1; }

# Set input VCF file
input_vcf=../${vcf_base}.decomp.norm.vcf

# Default Java memory=8g, unless memory is specified by SLURM
[[ ${SLURM_MEM_PER_NODE} ]] && \
mem="$(( 85 * $SLURM_MEM_PER_NODE / 100000 ))g" || \
mem=8g
echo "Java memory set to $mem"
# Run SnpEff annotation
cd snpEff
if [[ -f $input_vcf ]]
then
	echo "Running SnpEff on $input_vcf using $genome_base database build..."
	java -Xmx${mem} -jar snpEff.jar -v ${genome_base} $input_vcf > ${vcf_base}.ann.vcf
else
	echo "Error - $input_vcf file not detected."
	exit 1
fi
