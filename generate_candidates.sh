#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=500mb
#SBATCH --time=01:00:00
#SBATCH -J generate_candidates
#SBATCH -o %x.log

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

# Build multi-FASTA for BLASTing candidate genes in protein FASTA
cd $candidate_genes_build
# Add gene IDs to giant kelp genes from José (if needed)
array=($(ls $gk_dir))
for i in ${array[@]}
do
	echo "Adding '$i' IDs..."
	in_pep="gk_${i}_orfs.pep"
	[[ -f $in_pep ]] || cp ${gk_dir}/${i}/longest_orfs.pep $in_pep
	[[ $(grep ">${i}_" $in_pep) ]] || sed -i "s/>/>${i}_/g" $in_pep
done
# Copy Ectocarpus protein FASTA (if needed)
if [[ -f $ecto_prot_unzip ]]
then
	echo "Detected $ecto_prot_unzip"
elif [[ -f $ecto_prot ]]
then
	echo "Unzipping ${ecto_prot}..."
	gunzip $ecto_prot
else
	echo "Downloading ${ecto_prot_link}..."
	wget $ecto_prot_link
	echo "Unzipping ${ecto_prot}..."
	gunzip $ecto_prot
fi
# Extract meiotic-annotated protein IDs from Ectocarpus (if needed)
[[ -f $ecto_list ]] || \
grep ">" $ecto_prot_unzip | \
grep -i "mei" | \
sed "s/>//g" | \
sed "s/ .*//g" > $ecto_list

# Make $query FASTA (if needed)
cd ..
if [[ -f $query ]]
then
	echo "Detected $query - exiting. Delete $query to rerun this script."
	exit 0
else
	echo "Building ${query}..."
	# Add all giant kelp FASTA sequences
	cat ${candidate_genes_build}/gk*.pep > $query
	# Add all Ectocarpus FASTA sequences
	# Keep annotation info in IDs by replacing spaces with underscores
	seqtk subseq ${candidate_genes_build}/${ecto_prot_unzip} \
${candidate_genes_build}/${ecto_list} | sed "s/ /_/g" | \
sed "s/_(.*//g" >> $query
fi
