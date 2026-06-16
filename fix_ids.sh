#!/bin/bash
#SBATCH -J fix_ids
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p epyc-64
#SBATCH --time=05:00:00
#SBATCH --mem=25gb

# Correct FASTA IDs in protein.fa and cds.fa using genes.gtf
echo "\
---------------------FIX CDS AND PROTEIN FASTA IDS WITH GTF---------------------

IMPORTANT NOTE: This script will only perform the ID correction for 
annotations in which there is a 1 : 1 : 1 ratio bewteen the following:
  - transcript ID (e.g., mRNA_#)
  - CDS ID (e.g., jgi|SlaSLCT1FG3_1|#####|) 
  - protein ID (e.g., jgi|SlaSLCT1FG3_1|#####|)
Otherwise, it will fail with exit code 1.

This script should be smart enough not to fix the same FASTA file twice, 
but I'd recommend backing up your input FASTAs just in case. :)"
# genome_base="SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked"
genome_base="$1"
snpeff_db="snpEff/data/$genome_base"
# annot="$path_to_files/genes.gtf"
annot="$snpeff_db/genes.gff3"
prot="$snpeff_db/protein.fa"
cds="$snpeff_db/cds.fa"
ann_idx=annotation_index.txt

# Create annotation index from GTF file
echo "Creating $ann_idx from $annot..."
if [[ -f "$annot" ]]
then
	echo "Using $annot to create $ann_idx..."
	grep $'\tmRNA\t' "$annot" \
    | cut -f9 \
	  | cut -d';' -f2,1,6,4 --output-delimiter '=' \
	  | cut -d'=' -f1,3,5,7 --output-delimiter ' '\
	  | head -n1 \
	  > "$ann_idx"
	grep $'\tmRNA\t' "$annot" \
	  | cut -f9 \
	  | cut -d';' -f2,1,6,4 --output-delimiter '=' \
	  | cut -d'=' -f2,4,6,8 --output-delimiter ' ' \
	  >> "$ann_idx"
else
	echo "Error: $annot file not found."
	exit 1
fi

# Check input files for correctness
echo "Input protein FASTA: $prot"
echo "Input CDS FASTA: $cds"
echo "Validating input FASTAs against $ann_idx..."
total_transcripts="$(tail -n +2 "$ann_idx" | wc -l)"
total_prot_IDs="$(grep -c ">" "$prot")"
total_prot_mrnas="$(grep -c ">mRNA" "$prot")"
total_cds_IDs="$(grep -c ">" "$cds")"
total_cds_mrnas="$(grep -c ">mRNA" "$cds")"
# Verify that # of GTF transcripts and proteins are equivalent
if [[ "$total_transcripts" -ne "$total_prot_IDs" ]]
then
  echo "Error: Total transcripts =/= total proteins."
  echo "$ann_idx ($total_transcripts) =/= $prot ($total_prot_IDs)"
  exit 1
fi
# Verify that # of GTF transcripts and CDSs are equivalent 
if [[ "$total_transcripts" -ne "$total_cds_IDs" ]]
then
  echo "Error: Total transcripts =/= total CDS" 
  "$ann_idx ($total_transcripts) =/= $cds ($total_cds_IDs)"
  exit 1
fi

# Check if script has already been run
if [[ "$total_prot_IDs" -eq "$total_prot_mrnas" ]]
then
  echo "IDs in $prot have already been fixed. Exiting."
  exit 0
fi
if [[ "$total_cds_IDs" -eq "$total_cds_mrnas" ]]
then
  echo "IDs in $cds have already been fixed. Exiting."
  exit 0
fi

# Run ID fix program on $prot and $cds
# Grab header columns from $ann_idx
ann_header_prot="$(head -n1 "$ann_idx" | cut -d' ' -f3)"
echo "Prepending mRNA IDs to protein and CDS FASTA IDs..."
echo "Building JGI IDs using $ann_idx column: $ann_header_prot"
for ((i=1;i<=total_transcripts;i++))
do
	line="$(tail -n +2 "$ann_idx" | sed -n "${i}p")"
	mrna=">$(echo $line | cut -d' ' -f1)"
	id="jgi|SlaSLCT1FG3_1|$(echo "$line" | cut -d' ' -f3)|"
	if grep -q "$id" "$prot"
	then
	  sed -i "s/>$id/$mrna  $id/g" "$prot"
	else
	  echo "Error: $id not found in $prot"
	fi
	if grep -q "$id" "$cds"
	then
	  sed -i "s/>$id/$mrna  $id/g" "$cds"
	else
	  echo "Error: $id not found in $cds"
	fi
done
