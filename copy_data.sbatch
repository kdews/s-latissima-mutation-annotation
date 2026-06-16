#!/bin/bash
#SBATCH -J copy_data
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p epyc-64
#SBATCH --mem=30g
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12

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
# Fuction to copy or uncompress input files to ref_data/ in working directory
copy_or_uncompress() {
  local src="$1" # source file
  local dest_dir="$2" # destination directory
  if [[ ! -e "$src" ]] # verify source exists
  then
    echo "Error: File not found: $src"
    exit 1
  fi
  mkdir -p "$dest_dir" # create destination, if needed
  echo "Copying and uncompressing $src to $dest_dir/"
  src_basename_unzip="$(basename "$src" .gz)" # get unzipped file name
  if [[ "$src" == *.gz ]]
  then
    set -x
    gunzip -c "$src" > "$dest_dir/$src_basename_unzip"
    set +x
  else
    set -x
    rsync -ah "$src" "$dest_dir/"
    set +x
  fi
}

# Help message
show_help () {
  echo "\
Copies and indexes input VCF, and copies and uncompresses reference genome FASTA,
annotation file, protein FASTA, CDS FASTA to working directory.

Usage: sbatch $job_name.sbatch [config_file]
  config_file    (optional) provide config file if default config ($config_def)
                 does not exist in either scripts_dir ($scripts_dir)
                 or current directory (./)

Output:
  .
  ├── <vcf>                             # local VCF copy
  ├── <vcf>.tbi                         # local VCF index
  └── ref_data/                         # reference data directory
      ├── <genome>                      # uncompressed reference genome FASTA
      ├── <annot>                       # uncompressed reference annotation file
      ├── <protein>                     # uncompressed reference protein FASTA
      └── <cds>                         # uncompressed reference CDS FASTA

Requires:
  - BCFtools (https://www.htslib.org/)
  - local install in ./ of SnpEff (https://pcingola.github.io/SnpEff/)"
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

# Multithread support
if [[ -n "$SLURM_CPUS_PER_TASK" ]]
then
  nthrd="$SLURM_CPUS_PER_TASK"
else
  nthrd="1"
fi

# Load bcftools env
conda activate bcftools
bcftools --version || exit 1
echo

# Copy and uncompress all input files

# VCF
if [[ -f "$vcf" ]]
then
  rsync -ah "$vcf" .
else
  echo "Error: VCF not found: $vcf"
  exit 1
fi
# Index VCF
set -x
bcftools index -f --tbi --threads "$nthrd" "$vcf_basename"
set +x

# Reference genome FASTA
copy_or_uncompress "$genome" "ref_data/"

# Annotation file
copy_or_uncompress "$annot" "ref_data/"

# Protein FASTA
copy_or_uncompress "$prot" "ref_data/"

# CDS FASTA
copy_or_uncompress "$cds" "ref_data/"

# Done
echo
date
