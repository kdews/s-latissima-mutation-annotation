#!/bin/bash
#SBATCH -p epyc-64
#SBATCH -J copy_vcfs
#SBATCH -o %x_logs/%x_%j.log
#SBATCH --mem=30g
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12

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

# Help message
show_help () {
  echo "\
Copies input VCF and reference files annotation file (GFF3), protein FASTA,
CDS FASTA to working directory.

Usage: sbatch $job_name.sbatch

Output:
  .
  ├── $(basename "$vcf")                  # copied VCF file
  ├── $(basename "$vcf").tbi              # copied VCF index
  ├── $snpeff_db/                          # custom SnpEff database
  ├── $(basename "$genome")               # copied reference genome FASTA
  └── $(basename "$genome").fai           # copied reference genome FAIDX index

Requires: Samtools/BCFtools/HTSlib (https://www.htslib.org/)"
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

# Load conda
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

# Copy, unzip and rename all input files

# Copy VCF and VCF index
rsync -a "$vcf" "$vcf_idx" .

# Create directory for new genome and annotation files
echo "Creating SnpEff database: $snpeff_db"
mkdir -p "$snpeff_db"

# Copy genome FASTA, uncompress, rename to "sequences.fa"
echo "Copying $genome to sequences.fa"
rsync -a "$genome" "$snpeff_db/"
gunzip -c "$snpeff_db/$genome_basename" > "$snpeff_db/sequences.fa"
# Copy annotation file, uncompress, and rename to "genes.gtf"
echo "Copying $annot to genes.${annot_filetype}"
rsync -a "$annot" "$snpeff_db/"
gunzip -c "$snpeff_db/$annot_basename" > "$snpeff_db/genes.$annot_filetype"
# Convert input GFF3 files to GTF3 (if provided)
if [[ -f "$snpeff_db/genes.gtf" ]]
then
	echo "Using $snpeff_db/genes.gtf to build SnpEff database."
elif [[ $annot_filetype = "gff3" ]]
then
  conda activate agat
  set -x
	agat_convert_sp_gff2gtf.pl \
	  --gff "$snpeff_db/genes.gff3" \
	  --output "$snpeff_db/genes.gtf" > \
	  agat_convert_sp_gff2gtf.log 2>&1 
  set +x
	rm "$snpeff_db/genes.gff3"
	conda deactivate
else
	echo \
"Error - annotation file format $annot cannot be handled. This script only
accepts GFF3 or GTF3 files as input.
REQUIRED: GFF3 files MUST end in '.gff3', not '.gff'. Sowwy :("
	exit 1
fi
# Copy protein FASTA, uncompress, then rename to "protein.fa"
echo "Copying $prot to protein.fa"
rsync "$prot" "$snpeff_db/"
gunzip -c "$snpeff_db/$prot_basename" > "$snpeff_db/protein.fa"
# Copy CDS FASTA, uncompress, then rename to "cds.fa"
echo "Copying $cds to cds.fa"
rsync "$cds" "$snpeff_db/"
gunzip "$snpeff_db/$cds_basename" > "$snpeff_db/cds.fa"

# Done
echo
date
