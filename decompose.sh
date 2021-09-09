#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=25gb
#SBATCH --time=05:00:00
#SBATCH -J decompose
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

# Copy genome to working directory (if it doesn't exist)
if [[ -f $genome_basename_unzip ]]
then
	echo "Using genome: $genome_basename_unzip"
elif [[ -f $genone_basename ]]
then
	echo "Unzipping genome ${genome_basename}..."
else
	echo "Copying $genome to $(pwd)"
	rsync $genome .
	gunzip $genome_basename
fi

# Copy VCF to working directory (if it doesn't exist)
if [[ -f $vcf_basename_unzip ]]
then
	echo "Detected VCF file: ${vcf_basename_unzip}."
elif [[ -f $vcf_basename ]] 
then
	echo "Detected zipped VCF file: ${vcf_basename}. Unzipping..."
	gunzip $vcf_basename
else
	echo "Copying input VCF file $vcf to $(pwd)..."
	rsync --verbose --progress $vcf .
	echo "Unzipping ${vcf}..."
	gunzip $vcf_basename
fi

# Run vt decompose on input VCF
echo "Decomposing ${vcf_basename_unzip}..."
vt decompose -s \
-o ${vcf_base}.decomp.vcf \
$vcf_basename_unzip

# Run vt normalize on decomposed VCF
echo "Normalizing ${vcf_base}.decomp.vcf..."
vt normalize \
-r $genome_basename_unzip \
-o ${vcf_base}.decomp.norm.vcf \
${vcf_base}.decomp.vcf

