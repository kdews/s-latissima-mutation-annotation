#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=25gb
#SBATCH --time=1-0
#SBATCH -J fix_ids
#SBATCH -o %x.log

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
genome_base=$1
path_to_files=snpEff/data/${genome_base}
annot=${path_to_files}/genes.gtf
prot=${path_to_files}/protein.fa
cds=${path_to_files}/cds.fa
ann_idx=annotation_index.txt

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

# Check input files for correctness
echo "Validating input files $prot and $cds against ${ann_idx}..."
total_transcripts=$(tail -n +2 $ann_idx | wc -l)
total_prot_IDs=$(grep -c ">" snpEff/data/${genome_base}/protein.fa)
total_prot_mrnas=$(grep -c ">mRNA" snpEff/data/${genome_base}/protein.fa)
total_cds_IDs=$(grep -c ">" snpEff/data/${genome_base}/cds.fa)
total_cds_mrnas=$(grep -c ">mRNA" snpEff/data/${genome_base}/cds.fa)
# Verify that # of GTF transcripts and proteins are equivalent
[[ $total_transcripts -ne  $total_prot_IDs ]] && \
echo "Error - total transcripts in $ann_idx ($total_transcripts) =/= total \
proteins in $prot ($total_prot_IDs)" && exit 1
# Verify that # of GTF transcripts and CDSs are equivalent 
[[ $total_transcripts -ne  $total_cds_IDs ]] && \
echo "Error - total transcripts in $ann_idx ($total_transcripts) =/= total \
CDSs ($total_cds_IDs)" && exit 1

# Check that script has not already been run - if so, exit without error
[[ $total_prot_IDs -eq $total_prot_mrnas ]] && echo "IDs in $prot have already \
been fixed. Exiting." && exit 0
[[ $total_cds_IDs -eq $total_cds_mrnas ]] && echo "IDs in $cds have already \
been fixed. Exiting." && exit 0

# Run ID fix program on $prot and $cds
# Grab header columns from $ann_idx
ann_header_prot=$(head -n1 $ann_idx | awk '{print $3}')
ann_header_cds=$(head -n1 $ann_idx | awk '{print $4}')
echo "Prepending to $prot\nand $cds\nFASTA IDs with mRNA IDs from ${ann_idx} \
that correspond to ${ann_header_prot} and ${ann_header_prot}, respectively."
for i in $(seq 1 1 $total_transcripts)
do
	line=$(tail -n +2 $ann_idx | sed -n ${i}p)
	prot_id="jgi|SlaSLCT1FG3_1|$(echo $line | awk '{print $3}')|"
	cds_id="jgi|SlaSLCT1FG3_1|$(echo $line | awk '{print $4}')|"
	mrna=">$(echo $line | awk '{print $2}')"
	[[ $(grep "$prot_id" $prot) ]] && \
sed -i "s/>$prot_id/$mrna  $prot_id/g" $prot || \
echo "$prot_id not found in $prot"
	[[ $(grep "$cds_id" $cds) ]] && \
sed -i "s/>$cds_id/$mrna  $cds_id/g" $cds || \
echo "$cds_id not found in $cds"
done

