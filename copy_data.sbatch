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

# Load bcftools env
conda activate bcftools

# Copy and uncompress all input files

# VCF
rsync -ah "$vcf" .
# Index VCF
set -x
bcftools index --tbi --threads "$nthrd" "$vcf_basename"
sex +x

# Reference genome FASTA
echo "Copying and uncompressing $genome to ref_data/"
if [[ "$genome" == *.gz ]]
then
  set -x
  gunzip -c "$genome" > "ref_data/$genome_basename_unzip"
  set +x
else
  set -x
  rsync -ah "$genome" "ref_data/"
  set +x
fi

# Annotation file
echo "Copying and uncompressing $annot to ref_data/"
if [[ "$annot" == *.gz ]]
then
  set -x
  gunzip -c "$annot" > "ref_data/$annot_basename_unzip"
  set +x
else
  set -x
  rsync -ah "$annot" "ref_data/"
  set +x
fi

# Protein FASTA
echo "Copying and uncompressing $prot to ref_data/"
if [[ "$prot" == *.gz ]]
then
  set -x
  gunzip -c "$prot" > "ref_data/$prot_basename_unzip"
  set +x
else
  set -x
  rsync -ah "$prot" "ref_data/"
  set +x
fi

# CDS FASTA
echo "Copying and uncompressing $cds to ref_data/"
if [[ "$cds" == *.gz ]]
then
  set -x
  gunzip -c "$cds" > "ref_data/$cds_basename_unzip"
  set +x
else
  set -x
  rsync -ah "$cds" "ref_data/"
  set +x
fi

# Done
echo
date
