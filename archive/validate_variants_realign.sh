#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=5gb
#SBATCH --cpus-per-task=12
#SBATCH --time=2-0
#SBATCH -J validate_variants
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
echo "Threads set to ${threads}."; thread_options="-p $threads"; }

# Create alignment output directory
[[ -d $aln_dir ]] || mkdir $aln_dir

# Iterate through BED file to subset gene regions of interest from genome into 
# separate FASTA files, and build HISAT2 indices for each FASTA
for i in $(seq 2 1 $(cat $validate_bed | wc -l))
do
	gene_id=$(sed -n ${i}p $validate_bed | awk '{print $4}')
	if [[ -f ${aln_dir}/${gene_id}.fasta ]]
	then
		echo "Detected FASTA for $gene_id in ${aln_dir}."
	else
		echo "Subsetting $gene_id region from ${genome_basename_unzip}."
		sed -n 1p $validate_bed > temp.bed
		sed -n ${i}p $validate_bed >> temp.bed
		bedtools getfasta -name -fi $genome_basename_unzip \
-bed temp.bed -fo ${aln_dir}/${gene_id}.fasta
	fi
	if [[ -f ${aln_dir}/${gene_id}.1.ht2 ]]
	then
		echo "Detected HISAT2 index for $gene_id in ${aln_dir}."
	else
		echo "Building HISAT2 index for $gene_id in ${aln_dir}."
		hisat2-build $thread_options ${aln_dir}/${gene_id}.fasta \
${aln_dir}/$gene_id
	fi
done
[[ -f temp.bed ]] && rm temp.bed

# Submit alignment step
# Check if SLURM sbatch can be used for alignment submission
sbatch --help 
[[ $(echo $?) -eq 0 ]] && \
submit=sbatch || submit=bash
if [[ $submit == "sbatch" ]] 
then
	[[ -d hisat2_logs ]] || mkdir hisat2_logs
	if [[ $max_array_size ]]
	then
		array_size=$(cat $indiv_file | wc -l)
		extra=$(( $array_size % $max_array_size ))
		[[ $extra -ne 0 ]] && plus=1 || plus=0
		iter=$(( $(( $(( $array_size - $extra )) / $max_array_size )) \
+ $plus ))
		minus=$(( $max_array_size - 1 ))
		for i in $(seq 1 1 $iter)
		do
			n=$(( $max_array_size * $i ))
			m=$(( $n - $minus ))
			if [[ $i -eq $iter ]] && [[ $extra -ne 0 ]]
			then
				p=$(( $i - 1 ))
				n=$(( $(( $max_array_size * $p )) + $extra ))
				minus=$(( $extra - 1 ))
				m=$(( $n - $minus ))
			fi
			array_sets+=("$m-$n")
		done
	else
		array_size=$(cat $indiv_file | wc -l)
		array_sets=("1-$array_size")
	fi
	echo "array sets: ${array_sets[@]}"
	for i in $(seq 0 1 ${#array_sets[@]})
	do
		if [[ $i -eq 0 ]]
		then
			echo "Submitting alignment step $i with command: $submit \
-o hisat2_logs/%x_%a.log --array=${array_sets[i]} ${scripts_dir}hisat2.sh"
			$submit -o hisat2_logs/%x_%a.log \
--array=${array_sets[i]} ${scripts_dir}hisat2.sh
		else
			num=$(echo ${array_sets[i]} | sed "s/.*-//g")
			until [[ $(ls hisat2_logs/*checkpoint | wc -l) \
-eq $num ]]
			do
				echo "Waiting for array step \
$(( $i - 1 )) to complete..."
				sleep 30m
			done
			echo "Submitting alignment step $i with command: $submit \
-o hisat2_logs/%x_%a.log --array=${array_sets[i]} ${scripts_dir}hisat2.sh"
			$submit -o hisat2_logs/%x_%a.log \
--array=${array_sets[i]} ${scripts_dir}hisat2.sh
		fi
	done
elif [[ $submit == "bash" ]]
then
	$submit ${scripts_dir}hisat2.sh
else
	echo "Error - submission type unset for alignment step. Exiting..."
	exit 1
fi

