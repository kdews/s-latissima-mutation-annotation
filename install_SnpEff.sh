#!/bin/bash
#SBATCH -J install_SnpEff
#SBATCH -o %x_logs/%x_%j.log
#SBATCH -p epyc-64
#SBATCH --time=01:00:00
#SBATCH --mem=10gb

# Script variables
# Extract job name from script
if [[ -z "${SLURM_JOB_NAME}" ]]
then
  job_name="$(basename "$0")"
  job_name="${job_name%%.*}"
else
  job_name="$SLURM_JOB_NAME"
fi
scripts_dir="s-latissima-mutation-annotation/" # scripts directory
# By default looks in scripts_dir then current directory for "mut_annot.config"
config_def="mut_annot.config"
# Set Java memory from Slurm resources, if available
mem="${SLURM_MEM_PER_NODE:-"9638"}" # unit = MiB
mem="$(( 85 * mem / 100 ))m" # only use 85% of available (default = 8GiB)

# Help message
show_help () {
  echo "\
Downloads, unzips and tests latest SnpEff version into current directory.

Usage: sbatch $job_name.sbatch [config_file]
  config_file    (optional) provide config file if default config ($config_def)
                 does not exist in either scripts_dir ($scripts_dir)
                 or current directory ($pwd)
  
Output:
  .
  └── snpEff/                               # unzipped snpEff install directory
      ├── examples/
      │   ├── 1kg.head_chr1.filtered.vcf.gz
      │   ├── 1kg.head_chr1.vcf.gz
      │   ├── cancer.ann.vcf
      │   ├── cancer.eff.vcf
      │   ├── cancer_pedigree.ann.vcf
      │   ├── cancer_pedigree.vcf
      │   ├── cancer.vcf
      │   ├── example_motif.vcf
      │   ├── examples.sh
      │   ├── file.vcf
      │   ├── intervals.bed
      │   ├── my_annotations.bed
      │   ├── samples_cancer_one.txt
      │   ├── samples_cancer.txt
      │   ├── test.1KG.ann_encode.vcf
      │   ├── test.1KG.ann_reg.vcf
      │   ├── test.1KG.vcf
      │   ├── test.ann.vcf
      │   ├── test.chr22.ann.filter_missense_any_TRMT2A.vcf
      │   ├── test.chr22.ann.filter_missense_any.vcf
      │   ├── test.chr22.ann.filter_missense_first.vcf
      │   ├── test.chr22.ann.filter_missense.vcf
      │   ├── test.chr22.ann.one_per_line.txt
      │   ├── test.chr22.ann.txt
      │   ├── test.chr22.ann.vcf
      │   ├── test.chr22.vcf
      │   ├── test.vcf
      │   ├── variants_1.ann.vcf
      │   ├── variants_1.vcf
      │   ├── variants_2.ann.vcf
      │   └── variants_2.vcf
      ├── exec/
      │   ├── snpeff
      │   └── snpsift
      ├── galaxy/
      │   ├── snpEff_download.xml
      │   ├── snpEffWrapper.pl
      │   ├── snpEff.xml
      │   ├── snpSift_annotate.xml
      │   ├── snpSift_caseControl.xml
      │   ├── snpSift_filter.xml
      │   ├── snpSift_int.xml
      │   ├── snpSiftWrapper.pl
      │   ├── tool_conf.xml
      │   ├── tool-data/
      │   │   ├── snpEff_genomes.loc
      │   │   └── snpEff_genomes.loc.sample
      │   └── tool_dependencies.xml
      ├── LICENSE.md
      ├── scripts/
      │   ├── buildDbNcbi.sh
      │   ├── gsa/
      │   │   ├── bayesFactor_correction_scoreCount_max10.sh
      │   │   ├── bayesFactor_correction_scoreCount.r
      │   │   ├── bayesFactor_correction_scoreCount.sh
      │   │   ├── checkGeneNames.py
      │   │   ├── create_sets.bds
      │   │   ├── geneSetOverlap.py
      │   │   ├── geneSetOverlap.sort.txt
      │   │   ├── geneSetsGtex.py
      │   │   ├── pvalue_correction_scoreCount_min10.sh
      │   │   ├── pvalue_correction_scoreCount.r
      │   │   └── pvalue_correction_scoreCount.sh
      │   ├── make_dbNSFP.sh
      │   ├── _OLD/
      │   │   ├── 1kg.sh
      │   │   ├── annotate_demo_GATK.sh
      │   │   ├── annotate_demo.sh
      │   │   ├── bedEffOnePerLine.pl
      │   │   ├── cgShore.pl
      │   │   ├── cgShore.sh
      │   │   ├── countColumns.py
      │   │   ├── db.pl
      │   │   ├── extractSequences.pl
      │   │   ├── fasta2tab.pl
      │   │   ├── fastaSample.pl
      │   │   ├── fastaSplit.pl
      │   │   ├── fastqSplit.pl
      │   │   ├── filterBy.py
      │   │   ├── gffRemovePhase.pl
      │   │   ├── isutf8.py
      │   │   ├── join.pl
      │   │   ├── joinSnpEff.pl
      │   │   ├── make_and_commit_docs.sh
      │   │   ├── nextProt_filter.pl
      │   │   ├── ped2vcf.py
      │   │   ├── plotHistogram.pl
      │   │   ├── plotLabel.pl
      │   │   ├── plotMA.pl
      │   │   ├── plot.pl
      │   │   ├── plotQQ.pl
      │   │   ├── plotQQsubsample.pl
      │   │   ├── plotSmoothScatter.pl
      │   │   ├── plotXY.pl
      │   │   ├── queue.pl
      │   │   ├── sam2fastq.pl
      │   │   ├── sample.pl
      │   │   ├── snpSift_filter_sample_to_number.pl
      │   │   ├── sortLine.py
      │   │   ├── splitChr.pl
      │   │   ├── statsNum.pl
      │   │   ├── swapCols.pl
      │   │   ├── transpose.pl
      │   │   ├── txt2fa.pl
      │   │   ├── txt2vcf.py
      │   │   ├── uniqCount.pl
      │   │   ├── uniqCut.pl
      │   │   ├── vcfAnnFirst.py
      │   │   ├── vcfBareBones.pl
      │   │   ├── vcfEffHighest.ORI.py
      │   │   ├── vcfEffOnePerLine.pl
      │   │   ├── vcfFilterSamples.pl
      │   │   ├── vcfInfoOnePerLine.pl
      │   │   ├── vcfOnlyAlts.pl
      │   │   ├── vcfReduceGenotypes.pl
      │   │   ├── vcfRefCorrect.py
      │   │   └── wigSplit.pl
      │   ├── snpEff
      │   └── snpSift
      ├── snpEff.config
      ├── snpEff.jar
      ├── SnpSift.jar
      └── test.chr22.ann.vcf                                # test annotated VCF

Requires:
  - wget [tested with v1.19.5] (https://www.gnu.org/software/wget/)
  - UnZip [tested with v6.0] (https://infozip.sourceforge.net/UnZip.html)
  - Java >=v21 (https://openjdk.java.net)"
}
if [[ $1 == "-h" || $1 == "--help" ]]
then
  show_help
  exit 0
fi

# Print date and time
date
echo

# Source configuration file
# If scripts_dir exists, look there first
[[ -d "$scripts_dir" ]] && config_def="$scripts_dir/$config_def"
config_file="${1:-"$config_def"}"
if [[ -f "$config_file" ]]
then
	source "$config_file"
else
	echo "Error: Config file ($config_file) not found."
	exit 1
fi

# Load conda
if [[ -f "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Load SnpEff env
conda activate snpeff

# Test dependencies
wget --version || exit 1
echo
unzip -v || exit 1
echo
java --version || exit 1
echo

# Download latest SnpEff version
set -x
wget -nv https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip
unzip -o snpEff_latest_core.zip
rm snpEff_latest_core.zip
set +x
echo

# Test SnpEff installation
java -Xmx"$mem" -jar snpEff/snpEff.jar -version || exit 1
echo

# Run SnpEff on test data (creates snpEff/data directory)
echo "Testing SnpEff install and creating snpEff/data directory..."
set -x
java -Xmx"$mem" -jar snpEff/snpEff.jar GRCh37.75 \
  snpEff/examples/test.chr22.vcf > snpEff/test.chr22.ann.vcf
set +x

# Print date and time
echo
date
