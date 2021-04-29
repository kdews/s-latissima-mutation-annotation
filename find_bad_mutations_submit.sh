#!/bin/bash
#SBATCH -p cegs
#SBATCH -t 10:00:00
#SBATCH --mem 50gb
#SBATCH -o find_bad_mutations.out

source activate biopython
python scripts/s-latissima-mutation-annotation/find_bad_mutations.py assembly_ST_collapse_with_short_genes.complete_ORFs.cds_coords assembly_ST_collapse_with_short_genes.fa s_latissima_wgs_50_subset_FILT_qual20_SNPs_biallelic.calls.recode.vcf > SNP_biallelic_coding_stats.txt

