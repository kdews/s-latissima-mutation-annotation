#!/bin/bash
#SBATCH -J mut_annot_pipeline
#SBATCH -o %x_logs/%x_%j.log
#SBATCH --time=05:00:00
#SBATCH --mem=3gb

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

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || (( $# < 2 ))
then
  echo "\
Runs mutation annotation pipeline.
Usage:
  Shell
    bash $job_name.sh
    bash $job_name.sh
  Slurm
    sbatch $job_name.sh
    sbatch $job_name.sh
Requires:
  - SnpEff (https://pcingola.github.io/SnpEff/)
  - R (https://www.r-project.org)"
  exit 0
fi

# Print date and time
date
echo

# Run pipeline scripts
sbatch "$scripts_dir"/install_SnpEff.sh
sbatch "$scripts_dir"/copy_data.sbatch
sbatch "$scripts_dir"/build_SnpEff_db.sh
sbatch "$scripts_dir"/fix_missing.sbatch "$vcf"
sbatch "$scripts_dir"/ann_SnpEff.sh

# Print date and time
echo
date
