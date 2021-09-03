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

# Copy VCF to snpEff directory (if it doesn't exist)
cd snpEff
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

# Run SnpEff annotation
# Default Java memory=8g, unless memory is specified by SLURM
[[ ${SLURM_MEM_PER_NODE} ]] && \
mem="$(( 85 * $SLURM_MEM_PER_NODE / 100000 ))g" || \
mem=8g
echo "Java memory set to $mem"
echo "Running SnpEff on $vcf_basename_zip using $genome_base database build..."
java -Xmx${mem} -jar snpEff.jar -v ${genome_base} $vcf_basename_unzip > ${vcf_base}.ann.vcf

