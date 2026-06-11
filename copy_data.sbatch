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
latissima_dir="/project2/noujdine_61/kdeweese/latissima/"
popgen_dir="$latissima_dir/popgen_with_feeling/"
in_vcf="$popgen_dir/vcfs/master_2026_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.vcf.gz"
in_idx="$in_vcf.tbi"
genome="$latissima_dir/Assembly/Assembled_scaffolds__masked_/SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.fasta.gz"

# Help message
show_help () {
  echo "\
Copies input VCF and reference files annotation file (GFF3), protein FASTA,
CDS FASTA to working directory.

Usage: sbatch $job_name.sbatch

Output:
  .
  ├── $(basename "$in_vcf")          # input VCF file
  ├── $(basename "$in_vcf").tbi      # input VCF index
  ├── $(basename "$genome")          # reference genome FASTA
  └── $(basename "$genome").fai      # reference genome FAIDX index

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
cond=~/.conda_for_sbatch.sh
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

# Input
# List of VCF files
# vcf_list_file="vcf_list.txt"
vcf_list_file="$1"
# Reference FASTA
# ref="/project2/noujdine_61/kdeweese/latissima/Assembly/Assembled_scaffolds__masked_/SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.fasta.gz"
ref="$2"

# Output




# Copy, unzip and rename all input files, as necessary
# Create snpEff "data" directory by initiating random download (if needed)
[[ -d data ]] || java -jar snpEff.jar download -v GRCh37.75
# Create directory for new genome and annotation files
cd data
[[ -d $genome_base ]] || mkdir $genome_base
cd $genome_base
# Copy genome FASTA, uncompress, rename to "sequences.fa" (if needed)
[[ -f sequences.fa ]] || \
{ echo "Copying $genome to sequences.fa"; \
rsync $genome .; \
gunzip $genome_basename; \
mv $genome_basename_unzip sequences.fa; }
# Copy annotation file, uncompress, and rename to "genes.gtf" (if needed)
[[ -f genes.${annot_filetype} ]] || [[ -f genes.gtf ]] || \
{ echo "Copying $annot to genes.${annot_filetype}"; \
rsync $annot .; \
gunzip $annot_basename; \
mv $annot_basename_unzip genes.${annot_filetype}; }
# Convert input GFF3 files to GTF3 (if provided)
if [[ -f genes.gtf ]]
then
	echo "Using genes.gtf to build SnpEff database."
elif [[ $annot_filetype = "gff3" ]]
then
	agat_convert_sp_gff2gtf.pl --gff genes.gff3 --output \
genes.gtf > agat_convert_sp_gff2gtf.log 2>&1 
	[[ $? -ne 0 ]] && { echo "Error on AGAT conversion step. \
Exiting..."; exit 1; }
	rm genes.gff3
else
	printf "%s\n" "\
Error - annotation file format $annot cannot be handled. This script\n\
only accepts GFF3 or GTF3 files as input.\n\n\
REQUIRED: GFF3 files MUST end in '.gff3', not '.gff'. Sowwy :("
	exit 1
fi
# Copy protein FASTA, uncompress, then rename to "protein.fa" (if needed)
[[ -f protein.fa ]] || \
{ echo "Copying $prot to protein.fa"; \
rsync $prot .; \
gunzip $prot_basename; \
mv $prot_basename_unzip protein.fa; }
# Copy CDS FASTA, uncompress, then rename to "cds.fa" (if needed)
[[ -f cds.fa ]] || \
{ echo "Copying $cds to cds.fa"; \
rsync $cds .; \
gunzip $cds_basename; \
mv $cds_basename_unzip cds.fa; }





# Verify inputs exist
# Read VCF list file
if [[ -f "$vcf_list_file" ]]
then
  echo "Using VCF list: $vcf_list_file"
  mapfile -t in_vcfs < "$vcf_list_file"
  # Check for local indexed reference, or create if needed
  if [[ -f "$ref" ]]
  then
    echo "Reference: $ref"
    ref_basename="$(basename "$ref")"
    ref_base="${ref_basename%%.*}"
    ref_last_ext="${vcf_basename##*.}"
    ref_local="$ref_base.fasta"
    ref_idx="$ref_base.fasta.fai"
    if [[ -f "$ref_local" ]] && [[ -f "$ref_idx" ]]
    then
      echo "Local indexed reference found: $ref_local"
    elif [[ "$ref_last_ext" == "gz" ]]
    then
      echo "Unzipping and indexing reference..."
      set -x
      gunzip -k -c "$ref" > "$ref_local"
      set +x
      conda activate samtools
      samtools --version || exit 1
      echo
      set -x
      samtools faidx "$ref_local"
      set +x
      conda deactivate
    fi
  else
    echo "Error: reference FASTA ($ref) not found."
  fi
  echo
else
  echo "Error: could not find VCF list file ($vcf_list_file)."
  exit 1
fi

# Load BCFtools
conda activate bcftools
bcftools --version || exit 1
echo

# Loop through VCFs
for i in "${!in_vcfs[@]}"
do
  in_vcf="${in_vcfs[i]}"
  in_vcf="$(realpath -e "$in_vcf")"
  if [[ -f "$in_vcf" ]]
  then
    echo "Input VCF: $in_vcf"
  else
    echo "Error: Input VCF file ($in_vcf) not found."
    exit 1
  fi
  echo
  
  # Output
  # Name from input
  vcf_basename="$(basename "$in_vcf")"
  vcf_last_ext="${vcf_basename##*.}"
  
  # Don't bgzip compress if already bgzipped
  if [[ "$vcf_last_ext" == "gz" ]]
  then
    out_vcf="$vcf_basename"
  else
    out_vcf="$vcf_basename.gz"
  fi
  # Need tabix format for snpEff/SnpSift
  out_idx="$out_vcf.tbi"
  
  # Skip all steps if local indexed VCF exists
  if [[ -f "$out_vcf" ]] && [[ -f "$out_idx" ]]
  then
    echo "Local bgzipped VCF found: $out_vcf"
    echo "VCF index found: $out_idx"
  elif [[ "$vcf_last_ext" == "gz" ]]
  then
    echo "Skipping VCF compression."
    echo "Copying input VCF to $TMPDIR"
    cmd="rsync -htu --progress $in_vcf $TMPDIR/"
    echo "$cmd"
    $cmd
    echo
  else
    # Compress input VCF with bgzip
    echo "bgzip compressing $in_vcf..."
    cmd=(
        bgzip
        -@ "$nthrd"
        -o "$TMPDIR/$out_vcf"
        "$in_vcf"
    )
    echo "${cmd[*]}"
    "${cmd[@]}"
    echo
  fi
  
  # Only run if VCF was copied to TMPDIR
  if [[ -f "$TMPDIR/$out_vcf" ]]
  then
    # Index bgzipped VCF with bcftools
    cmd=(
        bcftools index
        --tbi
        -f
        --threads "$nthrd"
        -o "$TMPDIR/$out_idx"
        "$TMPDIR/$out_vcf"
    )
    echo "${cmd[*]}"
    "${cmd[@]}"
    echo

    # Copy output from /tmp to current directory, preserving timestamps (-t)
    echo "Copying output in $TMPDIR/ to $CURDIR/"
    rsync -htu --progress "$TMPDIR/"* "$CURDIR/"
    echo
  fi
  
  # Get metadata from each VCF
  echo "Extracting metadata from VCF..."
  year="$(bcftools head "$out_vcf" \
    | grep "##GATKCommandLine" \
    | sed -nE 's/.*Date="[^"]*([0-9]{4})[^"]*".*/\1/p' \
    | sort -u)"
  echo "Year: $year"
  if [[ -z "$year" ]]
  then
    echo "Error: Could not parse year from input VCF ($in_vcf)."
    exit 1
  elif [[ "$year" == 2026 ]]
  then
    trimmer="fastp"
  else
    trimmer="cutadapt"
  fi
  echo "Trimmer: $trimmer"
  version="$(bcftools head "$out_vcf" \
    | grep "##GATKCommandLine" \
    | sed -nE '/##GATKCommandLine/s/.*Version="([^"]*)".*/\1/p' \
    | sort -u)"
  echo "GATK version: $version"
  if [[ -z "$version" ]]
  then
    echo "Error: Could not parse GATK4 version from input VCF ($in_vcf)."
    exit 1
  fi
  
  # Save metadata to table
  echo "Saving VCF metadata to: $vcf_tab_file"
  if [[ "$i" -eq 0 ]]
  then
    printf "%s\t%s\t%s\t%s\t%s\n" "#VCF" "Reference" "Year" "Trimmer" "GATK_Version" > "$vcf_tab_file"
    printf "%s\t%s\t%s\t%s\t%s\n" "$out_vcf" "$ref_local" "$year" "$trimmer" "$version" >> "$vcf_tab_file"
  else
    printf "%s\t%s\t%s\t%s\t%s\n" "$out_vcf" "$ref_local" "$year" "$trimmer" "$version" >> "$vcf_tab_file"
  fi
  echo
done
echo

# Done
echo "Finished at: $(date)"
