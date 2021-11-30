#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=10mb
#SBATCH --time=01:00:00
#SBATCH -J zach_subset
#SBATCH -o %x.log

# Subset regions of interest from annotated VCF file 
# for Zach's rotation project

# Source conda configuration
conda_sh=~/bin/anaconda3/etc/profile.d/conda.sh
source $conda_sh

# Input VCF
in_vcf=snpEff/master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.ann.vcf
# Annotation file (GFF3)
ann=SlaSLCT1FG3_1_GeneCatalog_20210608.gff3
# List of gene IDs from pathways of interest
zach_genes=zach_gene_ids.txt
# Output files
zach_bed=zach_genes.bed
zach_vcf=zach_genes.ann.vcf

# Create BED file from gene list and annotation file
echo "Creating BED file $zach_bed from $zach_genes"
printf "%s\t%s\t%s\t%s\n" "#CHR" "START" "END" "GENE" > $zach_bed
for gene in $(cat $zach_genes)
do
	grep "ID=${gene};" "$ann" | awk '{print $1,$4,$5}'
	echo "$gene"
done | sed -z "s/\n/ /g" | sed "s/scaffold/\nscaffold/g" | \
tr -s " " | sed "s/ /\t/g" | tail -n +2 >> $zach_bed

# Subset regions of interest from annotated VCF
echo "Subsetting regions in $zach_bed from $in_vcf"
conda activate vcftools
vcftools --vcf $in_vcf --bed $zach_bed --recode \
--recode-INFO-all --out $zach_vcf
echo "Output in $zach_vcf"

