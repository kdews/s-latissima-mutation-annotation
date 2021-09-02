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

# Extract names from input files
genome_basename=$(basename -- $genome)
genome_basename_unzip=`echo $genome_basename | sed 's/\.gz//g'`
genome_base=`echo $genome_basename | sed 's/\..*//g'`
vcf_basename=$(basename -- $vcf)
vcf_basename_unzip=`echo $vcf_basename | sed 's/\.gz//g'`
vcf_base=`echo $vcf_basename | sed 's/\..*//g'`

# Optional: Anaconda configuration
# Attempt to source Anaconda from $conda_sh (if provided)
if [[ $conda_sh ]] && [[ -f $conda_sh ]]
then
	source $conda_sh
	[[ $? -eq 0 ]] && \
echo "Anaconda source successful." || \
{ echo "Error on Anaconda source from ${conda_sh}. Exiting..."; exit 1; }
	conda activate mut_annot
	[[ $? -eq 0 ]] && \
echo "Activation of conda env 'mut_annot' successful." || \
{ echo "Error activating conda env 'mut_annot'. Exiting..."; exit 1; }
else
	echo "conda.sh file not detected, expecting dependencies in PATH."
fi

# Test Java install (must be >= v1.8)
java --version
[[ $? -eq 0 ]] && echo "Java successfully loaded." || \
{ echo "Error arose while testing Java install. Exiting..."; exit 1; }

# Copy VCF to snpEff directory (if it doesn't exist)
cd snpEff
if [[ -f $vcf_basename_unzip ]]
then
	echo "Detected VCF file: $vcf_basename_unzip"
else
	echo "Copying input VCF file $vcf to $(pwd)..."
	rsync $vcf .
	gunzip $vcf_basename
fi

# Run SnpEff annotation
# Default Java memory=8g, unless memory is specified by SLURM
[[ ${SLURM_MEM_PER_NODE} ]] && \
mem="$(( 85 * $SLURM_MEM_PER_NODE / 1000 ))g" || \
mem=8g
java -Xmx${mem} -jar snpEff.jar ${genome_base} $vcf_basename_unzip > ${vcf_base}.ann.vcf

