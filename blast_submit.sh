#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=20gb
#SBATCH --time=08:00:00
#SBATCH -c 12

# Help message
if [[ $1 = "-h" ]] || [[ $1 = "--help" ]]; then
	echo \
"Usage:
sbatch blast_submit.sh <molecule_type> <query> <database>
  molecule_type    accepts either protein ('prot') or nucleotide ('nucl')
  query            must be FASTA file (can be gzipped)
  database         name of a cannonical BLAST database, e.g. 'nt' / 'nr'
                   (must be downloaded locally in a directory sourced by your
                   command-line BLAST installation; you can add the path to
                   your databases by editing ~/.ncbirc)

Runs only either blastn or blastp, depending upon <molecule_type>.

Requires NCBI command-line BLAST (https://www.ncbi.nlm.nih.gov/books/NBK52640)"
	exit 0
fi
# Source configuration file
[[ $1 ]] && config_file=$1 || config_file=mut_annot.config
if [[ $config_file ]] && [[ -f $config_file ]]

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

# Define positional arguments
molecule_type=$1
query=$2
database=$3
# Reformats queries given as filepaths and with file extensions
query_no_path_or_ext=$(basename -- $query | sed 's/\..*//g')

# Determine which blast command to use
if [[ $molecule_type == "nucl" ]]
then
	blast="blastn"
elif [[ $molecule_type == "prot" ]]
then
	blast="blastp"
else
	echo "Error, check that <molecule_type> is spelled correctly."
	exit 1
fi

# Handles gzipped input queries
if [[ $(echo $query | sed 's/.*\.//g') = "gz" ]]
then
	echo "Treating query as zipped FASTA."
	echo "Running $blast on $query with $database..."
	zcat $query | $blast -num_threads 12 -db $database \
-out ${query_no_path_or_ext}_vs_${database}.${molecule_type}.blast.tab \
-outfmt "6 qseqid sseqid pident length"
else
	echo "Running $blast on $query with $database..."
	$blast -num_threads 12 -query $query -db $database \
-out ${query_no_path_or_ext}_vs_${database}.${molecule_type}.blast.tab \
-outfmt "6 qseqid sseqid pident length"
fi

echo "$blast run finished. Results are in:
${query_no_path_or_ext}_vs_${database}.${molecule_type}.blast.tab"

