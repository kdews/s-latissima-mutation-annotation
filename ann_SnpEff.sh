#!/bin/bash
#SBATCH -J ann_SnpEff
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p epyc-64
#SBATCH --time=2-0
#SBATCH --mem=25gb

# Script variables
# Extract job name from script
if [[ -z "${SLURM_JOB_NAME}" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
scripts_dir="s-latissima-mutation-annotation/" # scripts directory
# By default looks in scripts_dir then current directory for "mut_annot.config"
config_def="mut_annot.config"
# Set Java memory from Slurm resources, if available
mem="${SLURM_MEM_PER_NODE:-"9638"}" # unit = MiB
mem="$(( 85 * mem / 100 ))m" # only use 85% of available (default = 8GiB)

# Help message
show_help () {
  echo "\
Runs SnpEff annotation with pre-built custom database.

Usage: sbatch $job_name.sbatch [config_file]
  config_file    (optional) provide config file if default config ($config_def)
                 does not exist in either scripts_dir ($scripts_dir)
                 or current directory ($pwd)
  
Output:
  .
  ├── <genome>.fasta      # unzipped reference genome FASTA
  ├── <genome>.fasta.fai  # reference genome FAIDX index
  ├── 
  ├── 
  └── xyz      # yadayada

Requires:
  - SnpEff (https://pcingola.github.io/SnpEff/)
  - Java >=v21 (https://openjdk.java.net)"
}
if [[ $1 == "-h" || $1 == "--help" ]]
then
  show_help
  exit 0
fi

# -e  exit on error
# -o pipefail  catch errors in pipes
set -eo pipefail

# Print date and time
date
echo

# Source configuration file
# If scripts_dir exists, look there first
[[ -d "$scripts_dir" ]] && config_def="$scripts_dir/$config_def"
config_file="${1:-"$config_def"}"
if [[ -f "$config_file" ]]
then
	source "$config_file"
else
	echo "Error: Config file ($config_file) not found."
	exit 1
fi


# Optional: Anaconda configuration
[[ $conda_sh ]] && source_conda $conda_sh

# Test Java install
java --version || exit 1

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
