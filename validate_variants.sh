#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=12
#SBATCH --time=10:00:00
#SBATCH -J bam_validate_variants
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

# Set multithreading (if specified)
[[ $SLURM_JOB_CPUS_PER_NODE ]] && { threads=$SLURM_JOB_CPUS_PER_NODE; \
echo "Threads set to ${threads}."; thread_options="-@ $threads"; }

# Create validatate BAMs output directory (if needed)
[[ -d $val_dir ]] || mkdir $val_dir

# Iterate through individuals file to subset gene regions of interest 
# from individual BAMs with lines from master BED file
for i in $(seq 1 1 $(cat $indiv_file | wc -l))
do
	gene_id=$(sed -n ${i}p $indiv_file | awk '{print $1}')
	indiv_id=$(sed -n ${i}p $indiv_file | awk '{print $2}')
	if [[ -f ${val_dir}/${gene_id}_${indiv_id}.subset.sorted.bam ]]
	then
		echo "Detected BAM subset for individual $indiv_id and gene \
$gene_id (${val_dir}/${gene_id}_${indiv_id}.subset.sorted.bam)."
	else
		echo "Subsetting $gene_id region from \
${bam_dir}/${indiv_id}.sorted.bam."
		sed -n 1p $validate_bed > temp.bed
		grep "${gene_id}$" $validate_bed >> temp.bed
		samtools view -b $thread_options -L temp.bed \
-o ${val_dir}/${indiv_id}_${gene_id}.subset.sorted.bam \
${bam_dir}/${indiv_id}_${genome_base}.sorted.marked.merged.bam
	fi
done
[[ -f temp.bed ]] && rm temp.bed

