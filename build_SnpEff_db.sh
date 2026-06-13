#!/bin/bash
#SBATCH -J build_SnpEff_db
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
Builds SnpEff database for a reference with given genome and annotation files.

Usage: sbatch $job_name.sbatch [config_file]
  config_file    (optional) provide config file if default config ($config_def)
                 does not exist in either scripts_dir ($scripts_dir)
                 or current directory (./)
  
Output:
  .
  ├── <genome_base>.build              # SnpEff database build log
  └── snpEff/                          # local SnpEff install directory
      ├── snpEff.config                # SnpEff config (updated with custom db)
      ├── data/                        # SnpEff 'dataDir' of databases
      │   └── <genome_base>/           # custom SnpEff database (post-build)
      │   │   ├── cds.fa               # renamed reference CDS FASTA
      │   │   ├── genes.<ext>          # renamed reference annotation file
      │   │   ├── protein.fa           # renamed reference protein FASTA
      │   │   └── sequences.fa         # renamed reference genome FASTA
      │   └── GRCh37.75/...            # GRCh37.75 example database (from test)
      └── ...

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

# Create directory for new genome and annotation files
echo "Building SnpEff database: $snpeff_db"
mkdir -p "$snpeff_db"
echo

# Copy genome to "sequences.fa"
echo "Copying ref_data/$genome_basename_unzip to sequences.fa"
cp "ref_data/$genome_basename_unzip" "$snpeff_db/sequences.fa"

# Copy annotation file and rename to "genes.$annot_filetype"
echo "Copying ref_data/$annot_basename_zip to genes.$annot_filetype"
cp "ref_data/$annot_basename_zip" "$snpeff_db/genes.$annot_filetype"

# # Convert input GFF3 files to GTF3 (if provided)
# if [[ -f "$snpeff_db/genes.gtf" ]]
# then
# 	echo "Using $snpeff_db/genes.gtf to build SnpEff database."
# elif [[ "$annot_filetype" == "gff3" ]]
# then
#   conda activate agat
#   set -x
# 	agat_convert_sp_gff2gtf.pl \
# 	  --gff "$snpeff_db/genes.gff3" \
# 	  --output "$snpeff_db/genes.gtf" > \
# 	  agat_convert_sp_gff2gtf.log 2>&1
#   set +x
# 	rm "$snpeff_db/genes.gff3"
# 	conda deactivate
# else
# 	echo \
# "Error: Annotation file format $annot cannot be handled. This script only
# accepts GFF3 or GTF3 files as input.
# 
# Required: GFF3 files MUST end in '.gff3', not '.gff'. Sowwy :("
# 	exit 1
# fi

# Copy protein FASTA to "protein.fa"
echo "Copying ref_data/$prot_basename_unzip to protein.fa"
cp "ref_data/$prot_basename_unzip" "$snpeff_db/protein.fa"

# Copy CDS FASTA to "cds.fa"
echo "Copying ref_data/$cds_basename_unzip to cds.fa"
cp "ref_data/$cds_basename_unzip" "$snpeff_db/cds.fa"

# Fix IDs in protein and CDS FASTAs using GTF file (if needed)
# Check if SLURM sbatch can be used for submission
sbatch --help
[[ $(echo $?) -eq 0 ]] && \
submit=sbatch || submit=bash
num_prot=$(grep -c ">" "snpEff/data/$genome_base/protein.fa")
num_cds=$(grep -c ">" "snpEff/data/$genome_base/cds.fa")
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
echo

# Add genome to "Non-standard Databases" section of snpEff.config file
# Retrieve line number of section header
line_num=$(grep -n "Non-standard Databases" snpEff/snpEff.config | cut -d : -f 1)
line_num=$(( $line_num + 2 ))
if grep -q "$species" snpEff/snpEff.config
then
  echo "Adding $species to SnpEff config (snpEff/snpEff.config)..."
  set -x
  sed -i \
    "$line_num a # $species genome, $genome_info\n${genome_base}.genome : ${species}\n" \
    snpEff/snpEff.config
  set +x
fi
echo

# Build SnpEff database for species genome (if it doesn't exist)
if [[ -f "$snpeff_db/$genome_base.build" ]]
then
	echo "Found previous database build: $snpeff_db/$genome_base.build"
elif [[ -d "$snpeff_db" ]] && grep -q "$species" snpEff/snpEff.config
then
  echo "Building SnpEff database for $species ($genome_base)..."
  set -x
  java -Xmx"$mem" -jar \
    snpEff/snpEff.jar build -v \
    -c snpEff/snpEff.config \
    -dataDir snpEff/data \
    -gff3 \
    "$genome_base" \
    2>&1 | "tee $genome_base.build"
  set +x
fi

