#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=500mb
#SBATCH --time=01:00:00
#SBATCH -J high_eff_parse
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

# Define input / output filenames
invcf=${vcf_base}.ann.gene_list.recode.vcf
outvcf=${vcf_base}.ann.gene_list.HIGH_EFF.recode.vcf
summ=high_eff.tab
simple_summ=high_eff.simplest.tab

# Save header to subsetted VCF file
grep "#" $invcf > $outvcf
# Extract all high effect variants
grep "|HIGH|" $invcf >> $outvcf

# Save header to summary .tab file
grep -v "##" $outvcf | grep "#" | awk '{$3=$6=$7=$8=$9=""; print $0}' > $summ
# Keep only position, genotype and effect fields from VCF
for i in $(seq 1 1 $(grep -v "#" $outvcf | wc -l))
do
	line=$(grep -v "#" $outvcf | sed -n ${i}p)
	echo "$line" | awk '{$3=$6=$7=$8=$9=""; print $0}'
	echo "$line" | awk '{print $8}' | sed 's/.*ANN=//g' | sed "s/,/\n/g" | \
awk 'BEGIN {FS="|"} {print $2,$7}'
done | sed -z "s/\n/ /g" | sed "s/scaffold/\nscaffold/g" | tr -s " " | \
sed "s/ /\t/g" >> $summ

# Translate VCF binary encoding to human-readable yes / no / NA
sed "s/\t1:/\tyes:/g" $summ | sed "s/\t0:/\tno:/g" | sed "s/\t\.:/\tNA:/g" | \
sed "s/:[0-9,:.]\+\t/\t/g" > $simple_summ
