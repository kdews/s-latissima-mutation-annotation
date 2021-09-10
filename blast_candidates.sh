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

Enviornment variable \$candidate_fasta in \$config_file must be set to 
multi-FASTA containing sequences of genes of interest from related species.

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
[[ $conda_sh ]] && source_conda $conda_sh

# Optional: Multithreading settings
[[ ${SLURM_CPUS_PER_TASK} ]] && threads=${SLURM_CPUS_PER_TASK} || threads=1

# Inputs
# Query
echo "Query set to: $query"
[[ -z $query ]] && { echo "Error - no input FASTA provided. Exiting..."; \
exit 1; }
# Database
echo "Input database file set to: $input_db"
[[ -z $input_db ]] && { echo "Error - no input BLAST database provided. \
Exiting..."; exit 1; }
echo "Database set to: $db"
[[ -z $db ]] && { echo "Error - no input BLAST database provided. \
Exiting..."; exit 1; }

# Determine which blast command to use
if [[ $molecule_type = "nucl" ]]
then
	blast="blastn"
elif [[ $molecule_type = "prot" ]]
then
	blast="blastp"
else
	echo "Error, check that <molecule_type> is spelled correctly."
	exit 1
fi

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
-out ${query_no_path_or_ext}_vs_${db}.${molecule_type}.blast.tab \
-outfmt "6 qseqid sseqid evalue"
else
	echo "Running $blast on $query with ${db}..."
	$blast -num_threads $threads -query $query -db $db \
-out ${query_no_path_or_ext}_vs_${db}.${molecule_type}.blast.tab \
-outfmt "6 qseqid sseqid evalue"
fi

[[ $? -eq 0 ]] && echo "$blast run finished. Results are in:
${query_no_path_or_ext}_vs_${db}.${molecule_type}.blast.tab"
