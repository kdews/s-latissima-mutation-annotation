#!/bin/bash
#SBATCH --mem=25gb
#SBATCH --time=1-0
#SBATCH -J halp_for_jose
#SBATCH -o %x.log

# Correct FASTA IDs in protein.fa and cds.fa using genes.gtf
echo "\
---------------------FIX CDS AND PROTEIN FASTA IDS WITH GTF---------------------
Usage: bash halp_for_jose.sh <annotation file>

Programs used:
AGAT v0.8.0 (https://github.com/NBISweden/AGAT)

IMPORTANT NOTE: This script will only perform the ID correction for 
annotations in which there is a 1 : 1 : 1 ratio bewteen the following:
  - transcript ID (e.g., mRNA_#)
  - CDS ID (e.g., jgi|SlaSLCT1FG3_1|#####|) 
  - protein ID (e.g., jgi|SlaSLCT1FG3_1|#####|)
Otherwise, it will fail with exit code 1.

This script should be smart enough not to fix the same FASTA file twice, 
but I'd recommend backing up your input FASTAs just in case. :)"

# LOAD FILES
# Input annotation file
annot=$1
annot_basename=$(basename -- $annot)
annot_basename_unzip=$(echo $annot_basename | sed 's/\.gz//g')
annot_filetype=$(echo $annot_basename_unzip |  sed 's/.*\.//g')
annot_noext=$(echo $annot_basename_unzip | sed 's/\..*//g')
# Output index file
ann_idx=annotation_index.txt

# CONVERT ANNOTATION FILE
# Convert input GFF3 files to GTF3 (if needed)
if [[ $annot_filetype = "gff3" ]]
then
	agat_convert_sp_gff2gtf.pl --gff $annot --output \
${annot_noext}.gtf > agat_convert_sp_gff2gtf.log 2>&1
	annot=${annot_noext}.gtf
	[[ $? -ne 0 ]] && { echo "Error on AGAT conversion step. \
Exiting..."; exit 1; }
elif [[ $annot_filetype = "gtf" ]]
	echo "File is GTF - skipping AGAT conversion step."
else
	echo "Error - annotation file format $annot cannot be handled. This script
only accepts GFF3 or GTF3 files as input.

REQUIRED: GFF3 files MUST end in '.gff3', not '.gff'. Sowwy :("
	exit 1
fi

# CREATE INDEX
# Create annotation index from GTF file
echo "Creating $ann_idx from ${annot}..."
if [[ -f $annot ]]
then
	echo "Using $annot to create ${ann_idx}..."
	grep $'\ttranscript\t' $annot | \
awk 'BEGIN { FS = "\t" } ; {print $9}' | \
awk 'BEGIN { FS = ";" } ; {print $1,$2,$7,$9}' | \
awk 'BEGIN { FS = " " } ; {print $1,$3,$5,$7}' | \
head -n1 > $ann_idx
	grep $'\ttranscript\t' $annot | \
awk 'BEGIN { FS = "\t" } ; {print $9}' | \
awk 'BEGIN { FS = ";" } ; {print $1,$2,$7,$9}' | \
sed 's/"//g' | \
awk 'BEGIN { FS = " " } ; {print $2,$4,$6,$8}' >> $ann_idx
else
	echo "Error - $annot file not found."
	exit 1
fi


## Only run this part if you want to edit FASTA IDs in $prot and $cds
#prot=protein.fa
#cds=cds.fa
## Check input files for correctness
#echo "Validating input files $prot and $cds against ${ann_idx}..."
#total_transcripts=$(tail -n +2 $ann_idx | wc -l)
#total_prot_IDs=$(grep -c ">" $prot)
#total_prot_mrnas=$(grep -c ">mRNA" $prot)
#total_cds_IDs=$(grep -c ">" $cds)
#total_cds_mrnas=$(grep -c ">mRNA" $cds)
## Verify that # of GTF transcripts and proteins are equivalent
#[[ $total_transcripts -ne  $total_prot_IDs ]] && \
#echo "Error - total transcripts in $ann_idx ($total_transcripts) =/= total \
#proteins in $prot ($total_prot_IDs)" && exit 1
## Verify that # of GTF transcripts and CDSs are equivalent 
#[[ $total_transcripts -ne  $total_cds_IDs ]] && \
#echo "Error - total transcripts in $ann_idx ($total_transcripts) =/= total \
#CDSs ($total_cds_IDs)" && exit 1
#
## Check that script has not already been run - if so, exit without error
#[[ $total_prot_IDs -eq $total_prot_mrnas ]] && echo "IDs in $prot have already \
#been fixed. Exiting." && exit 0
#[[ $total_cds_IDs -eq $total_cds_mrnas ]] && echo "IDs in $cds have already \
#been fixed. Exiting." && exit 0
#
## Run ID fix program on $prot and $cds
## Grab header columns from $ann_idx
#ann_header_prot=$(head -n1 $ann_idx | awk '{print $3}')
#ann_header_cds=$(head -n1 $ann_idx | awk '{print $4}')
#echo "Prepending to $prot\nand $cds\nFASTA IDs with mRNA IDs from ${ann_idx} \
#that correspond to ${ann_header_prot} and ${ann_header_prot}, respectively."
#for i in $(seq 1 1 $total_transcripts)
#do
#	line=$(tail -n +2 $ann_idx | sed -n ${i}p)
#	prot_id="jgi|SlaSLCT1FG3_1|$(echo $line | awk '{print $3}')|"
#	cds_id="jgi|SlaSLCT1FG3_1|$(echo $line | awk '{print $4}')|"
#	mrna=">$(echo $line | awk '{print $2}')"
#	[[ $(grep "$prot_id" $prot) ]] && \
#sed -i "s/>$prot_id/$mrna  $prot_id/g" $prot || \
#echo "$prot_id not found in $prot"
#	[[ $(grep "$cds_id" $cds) ]] && \
#sed -i "s/>$cds_id/$mrna  $cds_id/g" $cds || \
#echo "$cds_id not found in $cds"
#done

