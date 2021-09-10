#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=500mb
#SBATCH --time=01:00:00
#SBATCH -J parse_gene_list
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

# Copy and unzip annotation file (if needed)
if [[ -f $annot_basename_unzip ]]
then
	echo "Using annotation file $annot_basename_unzip (WARNING: Assuming GFF3)"
elif [[ -f $annot_basename ]]
then
	echo "Unzipping $annot_basename (WARNING: Assuming GFF3)"
	gunzip $annot_basename
else
	echo "Copying $annot to $(pwd)"
	rsync $annot .
	echo "Unzipping $annot_basename (WARNING: Assuming GFF3)"
	gunzip $annot_basename
fi

# Initialize gene list file with header
echo "Creating $gene_list in $(pwd)..."
echo "#Chr	Start	End	Source	E_value	gene_ID	\
transcript_ID	protein_ID	protein_product" > $gene_list

# Search for known patterns in annotation file
echo "Adding known reproductive annotations..."
for i in $(seq 1 1 $(grep -c -i "Spo11\|mei" $annot_basename_unzip))
do
	line=$(grep -i "Spo11\|mei" $annot_basename_unzip | sed -n ${i}p)
	coords=$(echo "$line" | awk '{print $1,$4,$5}' | sed "s/ /\t/g")
	ids=$(echo "$line" | awk 'BEGIN {FS="\t"} {print $9}')
	gene_id=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $1}' | sed "s/.*=//g")
	trans_id=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $6}' | sed "s/.*=//g")
	prot_id=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $5}' | sed "s/.*=//g")
	prot_product=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $4}' | sed "s/.*=//g")
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
"$coords" "annotation" "NA" "$gene_id" "$trans_id" "$prot_id" "$prot_product"
done >> $gene_list
# DE list
echo "Adding differentially expressed genes..."
type="transcript"
for i in $(cat $DE_list)
do
	num=$(echo "$i" | awk 'BEGIN {FS="|"} {print $3}')
	line=$(grep $'\tgene\t' $annot_basename_unzip | grep "${type}Id=${num}$")
	coords=$(echo "$line" | awk '{print $1,$4,$5}' | sed "s/ /\t/g")
	ids=$(echo "$line" | awk 'BEGIN {FS="\t"} {print $9}')
	gene_id=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $1}' | sed "s/.*=//g")
	trans_id=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $6}' | sed "s/.*=//g")
	prot_id=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $5}' | sed "s/.*=//g")
	prot_product=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $4}' | sed "s/.*=//g")
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
"$coords" "DE" "NA" "$gene_id" "$trans_id" "$prot_id" "$prot_product"
done >> $gene_list
# Input tabulated BLAST results
echo "Adding BLAST results..."
in_blast=${query_no_path_or_ext}_vs_${db}.${molecule_type}.blast.tab
type="protein"
for i in $(seq 1 1 $(cat $in_blast | wc -l))
do
	blast_line=$(sed -n ${i}p $in_blast)
	blast_source=$(echo "$blast_line" | awk '{print $1}')
	blast_eval=$(echo "$blast_line" | awk '{print $3}')
	# Test for low E-value for BLAST result
	[[ $(echo "$blast_eval" | grep "e\-") ]] || continue
	# Undo "|" replacements in FASTA IDs from BLAST run
	num=$(echo "$blast_line" | sed "s/-/|/g" | awk '{print $2}' | awk 'BEGIN {FS="|"} {print $3}')
	line=$(grep $'\tgene\t' $annot_basename_unzip | grep "${type}Id=${num};")
	coords=$(echo "$line" | awk '{print $1,$4,$5}' | sed "s/ /\t/g")
	ids=$(echo "$line" | awk 'BEGIN {FS="\t"} {print $9}')
	gene_id=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $1}' | sed "s/.*=//g")
	trans_id=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $6}' | sed "s/.*=//g")
	prot_id=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $5}' | sed "s/.*=//g")
	prot_product=$(echo "$ids" | awk 'BEGIN {FS=";"} {print $4}' | sed "s/.*=//g")
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
"$coords" "$blast_source" "$blast_eval" "$gene_id" "$trans_id" "$prot_id" "$prot_product"
done >> $gene_list

