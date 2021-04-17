#!/bin/bash
#SBATCH -p cegs
#SBATCH -t 00:20:00
#SBATCH --mem 10gb

python scripts/s-latissima-mutation-annotation/find_bad_mutations.py s_latissima_wgs_50_subset_FILT_qual20_SNPs_biallelic.calls.recode.vcf assembly_ST_collapse_with_short_genes.complete_ORFs.cds_coords assembly_ST_collapse_with_short_genes.fa

