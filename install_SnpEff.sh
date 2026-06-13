#!/bin/bash
#SBATCH -J install_SnpEff
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p epyc-64
#SBATCH --time=01:00:00
#SBATCH --mem=10gb

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
Downloads, unzips and tests latest SnpEff version into current directory.

Usage: sbatch $job_name.sbatch [config_file]
  config_file    (optional) provide config file if default config ($config_def)
                 does not exist in either scripts_dir ($scripts_dir)
                 or current directory (./)
  
Output:
  .
  └── snpEff/                     # uncompressed SnpEff install directory
      ├── data/                   # SnpEff 'dataDir' of databases
      │   └── GRCh37.75/...       # local GRCh37.75 example database
      ├── examples/               # SnpEff example data
      │   ├── test.chr22.vcf      # example VCF
      │   └── ...
      ├── exec/                   # dir containing only executable script links
      │   ├── snpeff              # alias for SnpEff wrapper in scripts/
      │   └── snpsift             # alias for SnpSift wrapper in scripts/
      ├── galaxy/...              # scripts for GATK & Galaxy integration
      ├── scripts/                # SnpEff core scripts
      │   ├── snpEff              # SnpEff wrapper script (handles Java)
      │   ├── snpSift             # SnpSift wrapper script (handles Java)
      │   └── ...
      ├── snpEff.config           # SnpEff config file
      ├── snpEff.jar              # SnpEff Java ARchive
      ├── SnpSift.jar             # SnpSift Java ARchive
      ├── test.chr22.ann.vcf      # test annotated VCF (generated after install)
      └── ...

Requires:
  - wget [tested with v1.19.5] (https://www.gnu.org/software/wget/)
  - UnZip [tested with v6.0] (https://infozip.sourceforge.net/UnZip.html)
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

# Test dependencies
wget --version || exit 1
echo
unzip -v || exit 1
echo
java --version || exit 1
echo

# Download latest SnpEff version
set -x
wget -N -nv https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip
unzip -o snpEff_latest_core.zip
rm snpEff_latest_core.zip
set +x
echo

# Test SnpEff installation
java -Xmx"$mem" -jar snpEff/snpEff.jar -version || exit 1
echo

# Run SnpEff on test data (creates snpEff/data directory)
echo "Testing SnpEff install and creating snpEff/data directory..."
mkdir -p test_snpEff

set -x
java -Xmx"$mem" -jar \
  snpEff/snpEff.jar ann -v \
  -c snpEff/snpEff.config \
  -stats test_snpEff/snpEff_summary_test.chr22.html \
  GRCh37.75 \
  snpEff/examples/test.chr22.vcf \
  > test_snpEff/test.chr22.ann.vcf
set +x

echo "SnpEff test results available in test_snpEff/"

# Print date and time
echo
date
