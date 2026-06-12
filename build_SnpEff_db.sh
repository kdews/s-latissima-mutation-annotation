#!/bin/bash
#SBATCH -J build_SnpEff_db
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
Builds SnpEff database for a reference with given genome and annotation files.

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

# Check for SnpEff installation and change into snpEff directory
if [[ -d snpEff ]]
then
	cd snpEff
else
	echo "Error: no snpEff directory detected."
	exit 1
fi
echo "Building SnpEff database in: $(pwd)"

# Load conda
if [[ -f "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Test Java install
java --version || exit 1

# Fix IDs in protein and CDS FASTAs using GTF file (if needed)
cd ../../..
# Check if SLURM sbatch can be used for submission
sbatch --help
[[ $(echo $?) -eq 0 ]] && \
submit=sbatch || submit=bash
num_prot=$(grep -c ">" snpEff/data/${genome_base}/protein.fa)
num_cds=$(grep -c ">" snpEff/data/${genome_base}/cds.fa)
[[ $num_prot -eq \
$(grep -c ">mRNA" snpEff/data/${genome_base}/protein.fa) ]] && \
[[ $num_cds -eq \
$(grep -c ">mRNA" snpEff/data/${genome_base}/cds.fa) ]] &&
echo "FASTA IDs already fixed." || \
$submit ${scripts_dir}fix_ids.sh $genome_base
until [[ $num_prot -eq \
$(grep -c ">mRNA" snpEff/data/${genome_base}/protein.fa) ]] && \
[[ $num_cds -eq \
$(grep -c ">mRNA" snpEff/data/${genome_base}/cds.fa) ]]
do
	echo "Waiting for fix_ids step to complete."
	date
	printf "%s\n" "Progress:\nprotein.fa - %s\ncds.fa - %s\n"
"$(grep -c '>mRNA' snpEff/data/${genome_base}/protein.fa) / $num_prot"
"$(grep -c '>mRNA' snpEff/data/${genome_base}/cds.fa) / $num_cds"
	sleep 3600
done

# Add genome to "Non-standard Databases" section of snpEff.config file
# Retrieve line number of section header
cd snpEff
line_num=$(grep -n "Non-standard Databases" snpEff.config | cut -d : -f 1)
line_num=$(( $line_num + 2 ))
[[ $(grep "$species" snpEff.config) ]] || \
{ echo "Adding $species to snpEff.config file..."; \
sed -i "$line_num a # $species genome, $genome_info\n\
${genome_base}.genome : ${species}\n" snpEff.config; }

# Build SnpEff database for species genome (if it doesn't exist)
if [[ -f ${genome_base}.build ]]
then
	echo "Using genome build ${genome_base}.build."
else
	[[ -d data/${genome_base} ]] && \
[[ $(grep "$species" snpEff.config) ]] && \
echo "Building SnpEff database for $species (${genome_base})..." && \
java -jar snpEff.jar build -v $genome_base 2>&1 | tee ${genome_base}.build
fi

