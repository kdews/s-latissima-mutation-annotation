#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=25gb
#SBATCH --time=2-0
#SBATCH -J build_SnpEff_db
#SBATCH -o %x.log

# Source configuration file
[[ $1 ]] && config_file=$1 || config_file=mut_annot.config
if [[ $config_file ]] && [[ -f $config_file ]]
then
	source $config_file
else
	echo "Error - please provide config file. \
$config_file not found."
	exit 1
fi

# Extract names from input files
genome_basename=$(basename -- $genome)
genome_basename_unzip=`echo $genome_basename | sed 's/\.gz//g'`
genome_base=`echo $genome_basename | sed 's/\..*//g'`
annot_basename=$(basename -- $annot)
annot_basename_unzip=`echo $annot_basename | sed 's/\.gz//g'`
annot_filetype=`echo $annot_basename_unzip |  sed 's/.*\.//g'`
prot_basename=$(basename -- $prot)
prot_basename_unzip=`echo $prot_basename | sed 's/\.gz//g'`
cds_basename=$(basename -- $cds)
cds_basename_unzip=`echo $cds_basename | sed 's/\.gz//g'`

# Check for SnpEff installation and change into snpEff directory
if [[ -d snpEff ]]
then
	cd snpEff
else
	echo "Error - no snpEff directory detected."
	exit 1
fi
echo "Building SnpEff database in $(pwd)"

# Optional: Anaconda configuration
# Attempt to source Anaconda from $conda_sh (if provided)
if [[ $conda_sh ]] && [[ -f $conda_sh ]]
then
	source $conda_sh
	[[ $? -eq 0 ]] && \
echo "Anaconda source successful." || \
{ echo "Error on Anaconda source from ${conda_sh}. Exiting..."; exit 1; }
	conda activate mut_annot
	[[ $? -eq 0 ]] && \
echo "Activation of conda env 'mut_annot' successful." || \
{ echo "Error activating conda env 'mut_annot'. Exiting..."; exit 1; }
else
	echo "conda.sh file not detected, expecting dependencies in PATH."
fi

# Test Java install (must be >= v1.8)
java --version
[[ $? -eq 0 ]] && echo "Java successfully loaded." || \
{ echo "Error arose while testing Java install. Exiting..."; exit 1; }

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

# Fix IDs in protein and CDS FASTAs using GTF file (if needed)
cd ../../..
num_prot=$(grep -c ">" snpEff/data/${genome_base}/protein.fa)
num_cds=$(grep -c ">" snpEff/data/${genome_base}/cds.fa)
[[ $num_prot -eq \
$(grep -c ">mRNA" snpEff/data/${genome_base}/protein.fa) ]] && \
[[ $num_cds -eq \
$(grep -c ">mRNA" snpEff/data/${genome_base}/cds.fa) ]] &&
echo "FASTA IDs already fixed." || \
sbatch ${scripts_dir}fix_ids.sh $genome_base
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

