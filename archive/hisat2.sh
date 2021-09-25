#!/bin/bash
#SBATCH -p cegs
#SBATCH --time=03:00:00
#SBATCH --mem=25g
#SBATCH -J hisat2
#SBATCH --cpus-per-task=12

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

# Define memory per thread for samtools sort (if threads set)
[[ $SLURM_JOB_CPUS_PER_NODE ]] && \
{ mem_per_thread="$(( 9 * $SLURM_MEM_PER_NODE / \
$SLURM_JOB_CPUS_PER_NODE / 10 ))M"; \
echo "Memory per thread set to ${mem_per_thread}."; \
hisat2_thread_options="-p ${SLURM_CPUS_PER_TASK}"; \
samtools_thread_options="-@ ${SLURM_CPUS_PER_TASK} -m ${mem_per_thread}"; }

if [[ $SLURM_ARRAY_TASK_ID ]]
then
	iterations=1
	# Define gene and individual IDs
	gene_id=$(awk '{print $1}' $indiv_file | sed -n ${SLURM_ARRAY_TASK_ID}p)
	indiv_id=$(awk '{print $2}' $indiv_file | \
sed -n ${SLURM_ARRAY_TASK_ID}p)
else
	iterations=$(cat $indiv_file | wc -l)
fi


for i in $(seq 1 1 $iterations)
do
	if [[ $iterations -gt 1 ]]
	then
		gene_id=$(awk '{print $1}' $indiv_file | sed -n ${i}p)
		indiv_id=$(awk '{print $2}' $indiv_file | sed -n ${i}p)
	fi
	# Define list of left and right reads for an individual
	# Use realpaths for HISAT2 inputs because directory changes
	reads1_list=($(ls ${reads_dir}/*${indiv_id}_*_R1*))
	reads2_list=($(ls ${reads_dir}/*${indiv_id}_*_R2*))
	[[ $reads1_list ]] || { echo "Error - right hand reads not \
found for ${indiv_id}. Exiting..." ; exit 1; }
	[[ $reads2_list ]] || { echo "Error - left hand reads not \
found for ${indiv_id}. Exiting..." ; exit 1; }
	echo "$gene_id $indiv_$id ${reads1_list[@]} ${reads2_list[@]}"
	[[ ${#reads1_list[@]} -ne ${#reads2_list[@]} ]] && \
{ echo "Error - at least one of a set of reads not found for \
${indiv_id}. Exiting..."; exit 1; }
	for i in $(seq 0 1 $(( ${#reads1_list[@]} - 1 )))
	do
		reads1=$(realpath ${reads1_list[i]})
		reads2=$(realpath ${reads2_list[i]})
		echo "Aligning reads associated with ${indiv_id}
($reads1 and
${reads2})
to ${gene_id} using HISAT2."
		sample_id=$(basename -- $reads1 | sed "s/_R1.*//g")
		# Define read groups for SAM headers from input $sample_id
		# and input FASTQ file read IDs
		# Input sample ID 
		# format = UniqueID_SampleID_Sequencer_Barcode_Plate_Lane_Read 
		# Create array 'rgs' from $sample_id by splitting on '_'
		rgs1=(${sample_id//_/ })
		read_id=$(zcat $reads1 | head -n1 | awk '{print $1}' | \
sed 's/@//g')
		rgs2=(${read_id//:/ })
		# Index both $sample_id array and read ID array to generate RGs
		# Read group ID = {FLOWCELL_BARCODE}.{LANE}
		rg_id="${rgs2[2]}.${rgs2[3]}"
		echo "Read group ID = ID:${rg_id}"
		# Read group SaMple name (SM) = shorten SampleID in $sample_id
		rg_sm="${rgs1[1]}"
		echo "Read group sample name = SM:${rg_sm}"
		# Read group Platform Unit (PU) 
		# format = {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_NAME}
		rg_pu="${rgs2[2]}.${rgs2[3]}.${rgs1[1]}"
		echo "Read group platform unit = PU:${rg_pu}"
		# Read group LiBrary identifier = SampleID-UniqueID
		rg_lb="${rgs1[1]}-${rgs1[0]}"
		echo "Read group library identifier = LB:${rg_lb}"
		# Read group PLatform (PL) = ILLUMINA
		rg_pl="ILLUMINA"
		echo "Read group platform = PL:${rg_pl}"
		# Run HISAT2 & pipe output to samtools
		hisat2 $hisat2_thread_options --time \
-x ${aln_dir}/$gene_id -1 ${reads1} -2 ${reads2} \
--rg-id ${rg_id} --rg SM:${rg_sm} --rg LB:${rg_lb} \
--rg PU:${rg_pu} --rg PL:${rg_pl} \
--summary-file=${aln_dir}/${sample_id}_${gene_id}.hisat2.summary | \
samtools sort $samtools_thread_options \
-o ${aln_dir}/${sample_id}_${gene_id}.sorted.bam
	done
done

# Create checkpoint
touch hisat2_logs/hisat2_${SLURM_ARRAY_TASK_ID}.checkpoint
echo "Run successfully completed and checkpoint created."
exit 0
