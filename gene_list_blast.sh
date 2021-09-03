#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=20gb
#SBATCH --time=08:00:00
#SBATCH -c 12
#SBATCH -J gene_list_blast
#SBATCH -o %x.log

# Help message
if [[ $1 = "-h" ]] || [[ $1 = "--help" ]]
then
	echo \
"For a given set of FASTA files of genes of interest, runs BLASTp with a 
protein database generated from a protein FASTA of reference. Reports the 
best hit for each FASTA. 

Usage:
sbatch gene_list_blast.sh [config_file]
  [config_file]      must either exist in current directory with name 
                     'mut_annot.config' or be specified by the first 
                     positional parameter (i.e., \$1)

Enviornment variable \$fastas in \$config_file must be set to multi-FASTA 
containing sequences of genes of interest from related species.

Requires NCBI command-line BLAST (https://www.ncbi.nlm.nih.gov/books/NBK52640)"
	exit 0
fi

# Source configuration file
[[ $1 ]] && config_file=$1 || config_file=mut_annot.config
if [[ $config_file ]] && [[ -f $config_file ]]
then
	source $config_file
else
	echo "Error - please provide config file. $config_file not found."
	exit 1
fi

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

# Define and create output directory (if needed)
outdir=blast_results
[[ ! -d $outdir ]] && mkdir $outdir

# Input filenames
[[ $fastas ]] || \
{ echo "Error - no input FASTA provided. Exiting..."; \
exit 1; }
molecule_type='prot'
# Reformats databases given as filepaths and with file extensions
input_db=$prot
input_db_basename=$(basename -- $input_db)
input_db_basename_unzip=$(echo $input_db_basename | sed 's/\.gz//g')
db=$(echo $input_db_basename | sed 's/\..*//g')
db_ext=$(echo $input_db_basename | sed 's/.*\.//g')

# Create BLAST database (if needed)
echo "Checking BLAST database ${db}..."
blastdbcmd -info -db $db
if [[ $? -ne 0 ]] && [[ -f $input_db ]]
then
	echo "Creating BLAST database $db from $input_db in current directory."
	echo "Copying $input_db to current directory."
	rsync --verbose --progress $input_db .
	[[ $db_ext = "gz" ]] && [[ ! -f $input_db_basename_unzip ]] \
&& gunzip $input_db_basename
	# Replace "|" with "-" in FASTA IDs (causes BLAST error)
	sed -i "s/|/-/g" $input_db_basename_unzip
	makeblastdb -in $input_db_basename_unzip -out $db -dbtype \
$molecule_type -title "$db"
	[[ $? -ne 0 ]] && { echo "Error creating database"; exit 1; }
elif [[ $? -ne 0 ]] && [[ ! -f $input_db ]]
then
	echo "Error - file $input_db not found. Database not created."
	exit 1
else
	echo "Using BLAST database ${db}."
fi

# Handles gzipped input queries
if [[ $(echo $query | sed 's/.*\.//g') = "gz" ]]
then
	echo "Treating query as zipped FASTA."
	echo "Running $blast on $query with ${db}..."
	zcat $query | $blast -num_threads $threads -db $db \
-out ${outdir}${query_no_path_or_ext}_vs_${db}.${molecule_type}.blast.tab \
-outfmt "6 qseqid sseqid pident length"
else
	echo "Running $blast on $query with ${db}..."
	$blast -num_threads $threads -query $query -db $db \
-out ${outdir}${query_no_path_or_ext}_vs_${db}.${molecule_type}.blast.tab \
-outfmt "6 qseqid sseqid pident length"
fi

echo "$blast run finished. Results are in:
${query_no_path_or_ext}_vs_${db}.${molecule_type}.blast.tab"
