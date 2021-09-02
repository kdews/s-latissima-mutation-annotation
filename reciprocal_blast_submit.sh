#!/bin/bash

# Help message
if [[ $1 = "-h" ]] || [[ $1 = "--help" ]]
then
	echo \
"A ready-made reciprocal BLAST recipe, just for you.

Usage:
bash reciprocal_blast_submit.sh <molecule_type> <query1> <query2> [db1] [db2]
  molecule_type    accepts either protein ('prot') or nucleotide ('nucl')
  query1           must be FASTA file (can be gzipped)
  query2           must be FASTA file (can be gzipped)
  db1              BLAST database created from query 1 with makeblastdb
  db2              BLAST database created from query 2 with makeblastdb
First, checks for existence of BLAST databases built from each query or given
as input, and creates databases if necessary. Then, runs reciprocal BLASTs 
(using either blastn or blastp depending on <molecule_type>).

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

# Define positional arguments
# Data molecule type
molecule_type=$1
# Input filenames
query1=$2
query2=$3
# If no database names given, parse from query filenames
[[ -z "$4" ]] && db1=$(echo $query1 | sed 's/\..*//g') || db1=$4
[[ -z "$5" ]] && db2=$(echo $query2 | sed 's/\..*//g') || db2=$5

# Create BLAST databases, if needed
source activate $blast_env
echo "Checking for BLAST databases $db1 and ${db2}..."
blastdbcmd -info -db $db1
if [[ $? -ne 0 ]]
then
	echo "Creating BLAST database $db1 from ${query1} in current directory."
	db1=$(echo $query1 | sed 's/\..*//g')
	makeblastdb -in $query1 -out $db1 -dbtype $molecule_type -parse_seqids
else
	echo "Using BLAST database $db1."
fi
blastdbcmd -info -db $db2
if [[ $? -ne 0 ]]
then
	echo "Creating BLAST database $db2 from ${query2} in current directory."
	db2=$(echo $query2 | sed 's/\..*//g')
	makeblastdb -in $query2 -out $db2 -dbtype $molecule_type -parse_seqids
else
	echo "Using BLAST database $db2."
fi

# Submit BLAST alignment jobs to SLURM
# only if blast_submit.sh output doesn't already exist
if [[ -f ${db1}_vs_${db2}.${molecule_type}.blast.tab ]]
then
	echo "Skipping BLAST of ${db1}_vs_${db2} (output exists). 
Remove ${db1}_vs_${db2}.${molecule_type}.blast.tab to rerun BLAST."
else
	sbatch -J ${db1}_vs_${db2} \
-o blast_submit_${db1}_vs_${db2}.${molecule_type}.out \
${scripts}/blast_scripts/blast_submit.sh $molecule_type $query1 $db2
fi

if [[ -f ${db2}_vs_${db1}.${molecule_type}.blast.tab ]]
then
	echo "Skipping BLAST of ${db2}_vs_${db1} (output exists). 
Remove ${db2}_vs_${db1}.${molecule_type}.blast.tab to rerun BLAST."
else
	sbatch -J ${db2}_vs_${db1} \
-o blast_submit_${db2}_vs_${db1}.${molecule_type}.out \
${scripts}/blast_scripts/blast_submit.sh $molecule_type $query2 $db1
fi


