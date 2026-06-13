#!/bin/bash
#SBATCH -J ann_SnpEff
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p epyc-64
#SBATCH --time=2-0
#SBATCH --mem=25gb

# Script variables
# Extract job name from script
if [[ -z "$SLURM_JOB_NAME" ]]
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

Usage: sbatch $job_name.sbatch[config_file]
  config_file    (optional) provide config file if default config ($config_def)
                 does not exist in either scripts_dir ($scripts_dir)
                 or current directory (./)
  
Output:
  .
  ├── <in_vcf>.ann.vcf         # SnpEff-annotated VCF
  ├── snpEff_genes.txt         # SnpEff html
  └── snpEff_summary.html      # SnpEff HTML

Requires:
  - SnpEff (https://pcingola.github.io/SnpEff/)
  - Java >=v21 (https://openjdk.java.net)"
}
if [[ $1 == "-h" || $1 == "--help" ]]
then
  show_help
  exit 0
fi

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

# Load conda
module purge
if [[ -f "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Load SnpEff env
conda activate snpeff

# Test Java install
java --version || exit 1

# Check for local SnpEff installation directory
if [[ ! -d snpEff ]]
then
	echo "Error: no snpEff directory detected."
	exit 1
fi

# Set input VCF file
input_vcf="$vcf_base.decomp.norm.vcf"

# Run SnpEff annotation
if [[ -f $input_vcf ]]
then
	echo "Running SnpEff on $input_vcf using $genome_base database build..."
	set -x
	java -Xmx"$mem" -jar \
	  snpEff/snpEff.jar ann -v \
	  -c snpEff/snpEff.config \
	  -dataDir snpEff/data \
	  -stats snpEff_summary_$vcf_base.html \
	  "$genome_base" \
	  "$input_vcf" \
	  > "$vcf_base.ann.vcf"
	  set +x
else
	echo "Error: VCF file not detected: $input_vcf"
	exit 1
fi
