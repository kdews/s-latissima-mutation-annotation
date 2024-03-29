# Mutation annotation of the brown macroalgae *Saccharina latissima* (North American sugar kelp)
A pipeline to detect and annotate the effects of deleterious mutations in *S. latissima* whole genome sequencing data.

## Installation
Clone and enter repository:
```
git clone https://github.com/kdews/s-latissima-mutation-annotation.git
cd s-latissima-mutation-annotation
```

## Dependencies
* [bash](https://www.gnu.org/software/bash) shell (default on many machines; see [Note on shells](#note-on-shells))
* [Java](https://openjdk.java.net) (>=v1.8) (also installed by default on many machines, but provided in YAML for the unlucky)
* [vt](https://github.com/atks/vt) (Bioconda version 2015.11.10 or equivalent)
* [SnpEff](https://pcingola.github.io/SnpEff) (v5.0c)
* [R](https://www.r-project.org) (v4.1.1)
	* [tidyverse](https://www.tidyverse.org) (v1.3.1)
* [bedtools](https://bedtools.readthedocs.io) (v2.30.0)
### Optional 
* [SLURM](https://slurm.schedmd.com/download.html)
> For submitting scripts with SLURM `sbatch`; can allow for higher memory and time allocations
* [AGAT](https://github.com/NBISweden/AGAT) (v0.8.0)
> AGAT installation not required if annotation file is GTF

##### Note on shells
I know for certain that sourcing the configuration file breaks without full paths if `sh` is used instead of `bash` - if you must use `sh` or another shell, be sure to give the output of `realpath mut_annot.config` as the first positional argument ($1) to the steps of the pipeline.

#### Easy mode: Create Anaconda environment from provided YAML, and download latest SnpEff release as directory within your working directory
1. If you haven't already, install the lastest version of [Anaconda](https://www.anaconda.com/).
##### Note on Anaconda
This pipeline will attempt to source Anaconda from the `$conda_sh` environment variable set in `mut_annot.config`. **Please set** `$conda_sh` to `path/to/<anaconda-version>/etc/profile.d/conda.sh` in `mut_annot.config`. If you are installing Anaconda for the first time, now would be a good time to find this path.

2. Create Anaconda env from `mut_annot.yml`.
```
conda env create -f mut_annot.yml
```
3. Download SnpEff to your working directory.
```
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```

> @ Kelly: remember to add final mut\_annot.yml at end lol~

#### Hard mode: Install programs manually and add to $PATH, and download latest SnpEff release as directory within your working directory
1. Install each dependency manually, and ensure that all programs (except SnpEff) have been added to your $PATH, e.g., running `which java` produces `/path/to/java`. 
```
export PATH=$PATH:/path/to/directory/containing/java
export PATH=$PATH:/path/to/directory/containing/agat
export PATH=$PATH:/path/to/directory/containing/vt
export PATH=$PATH:/path/to/directory/containing/R
```
If you are not comfortable doing this step, I suggest [Easy mode](#easy-mode-create-anaconda-environment-from-provided-yaml-and-download-latest-snpeff-release-as-directory-within-your-working-directory)^^.

2. Download SnpEff to your working directory.
```
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```

## Configuration
Edit `mut_annot.config` with paths to:
* Anaconda conda.sh file (i.e., `path/to/<anaconda-version>/etc/profile.d/conda.sh`)
* Reference genome
* Annotation file (GTF or GFF3)
* CDS FASTA
* Protein FASTA
* VCF file

Each step of the pipeline will first take the path to a config file as the first positional argument ($1); if one is not provided, it will then look for `mut_annot.config` in your current directory. ***This pipeline REQUIRES that you edit the configuration file with your local paths and files. It will break if you do not edit the configuration file.***

## Running the pipeline
The entire pipeline can be run either with `bash` or with SLURM's `sbatch`. Each script of the pipeline has built-in options to SLURM, but you can modify these with `sbatch [options]` (recommended) or by editing the `#SBATCH` headers of each script (not recommended). See the [documentation](https://slurm.schedmd.com/sbatch.html) for more information on `sbatch` options. 
> Using SLURM to execute the pipeline is recommended, if possible, as some of these scripts can run for upwards of 10 hours and require memory allocations >8g.

In the directory *containing* the `snpEff/` directory you [just downloaded](#easy-mode-create-anaconda-environment-from-provided-yaml-and-download-latest-snpeff-release-as-directory-within-your-working-directory) (i.e., the directory just *above* it), which would be `/path/to/s-latissima-mutation-annotation/`, if you are running from within the repository (recommended):
### 1. Run vt to decompose VCF file 
```
bash/sbatch [options] decompose.sh [/path/to/mut_annot.config]
```
### 2. Build SnpEff database for *S. latissima*
```
bash/sbatch [options] build_SnpEff_db.sh [/path/to/mut_annot.config]
```
This step uses the input *S. latissima* annotation file to correlate the transcript IDs in the CDS and protein FASTAs wih the annotation file. You can view the matrix extracted from the annotation file, `annotation_index.txt`, to inspect the relationship between IDs in the annotation file and the FASTAs.

### 3. Run SnpEff variant annotation on *S. latissima* VCF
```
bash/sbatch [options] ann_SnpEff.sh [/path/to/mut_annot.config]
```
### 4. Generate list of candidate genes to BLAST with
```
bash/sbatch [options] generate_candidates.sh [/path/to/mut_annot.config]
```
### 5. BLAST candidate translated genes against *S. latissima* protein FASTA
```
bash/sbatch [options] blast_candidates.sh [/path/to/mut_annot.config]
```
### 6. Parse "gene list" file from BLAST results (and any other sources)
```
bash/sbatch [options] gene_list_parse.sh [/path/to/mut_annot.config]
```
### 7. Extract variants in regions containing genes of interest
```
bash/sbatch [options] vcf_extract.sh [/path/to/mut_annot.config]
```
### 8. Extract only variants with high effects in regions of interest 
```
bash/sbatch [options] high_eff_parse.sh [/path/to/mut_annot.config]
```
### 9. Annotate and filter high effect varaint summary 
Adds gene and annotation information to `high_eff.annot.tab`, and filters for variants present in at least one each male and female gametophyte.
```
bash/sbatch [options] sterile_genotyping.sh [/path/to/mut_annot.config]
```
### 10. Manually validate variants of interest
```
bash/sbatch [options] validate_variants.sh [/path/to/mut_annot.config]
```

## Results
### Summary of High Effect Variants
```
high_eff_annot.tab
```
Annotated summary of high effect variants and individuals affected.
### SnpEff
Results of this analysis can be found in the `snpEff` directory.
#### SnpEff annotated VCF file
```
<your_vcf_name>.ann.vcf
```
SnpEff annotated VCF file. See the [documentation](https://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files) of SnpEff values in the ANN field to understand the putative effect(s) of each variant in your VCF file.
#### SnpEff HTML summary file
```
snpEff_summary.html
```
HTML file with extensive summary of SnpEff annotations in your VCF file. Here are some example figures:

1. Plot of variant frequency by region (relative to genes)

![alt text](https://github.com/kdews/s-latissima-mutation-annotation/blob/main/images/variant_type_freqs.png)

2. Summary of functional classes of all variants

![alt text](https://github.com/kdews/s-latissima-mutation-annotation/blob/main/images/functional_class.png)

See the [documentation](https://pcingola.github.io/SnpEff/se_outputsummary/#html-summary-snpeff_summaryhtml) of the SnpEff summary HTLM for more information.
#### SnpEff gene counts summary file
```
snpEff_genes.txt
```
Tab-delimited text file with counts of number of variants affecting each transcript and gene in the reference. See [documentation](https://pcingola.github.io/SnpEff/se_outputsummary/#gene-counts-summary-snpeff_genestxt).

## Debugging & Troubleshooting
The output at each step of the pipeline can be saved to a log file by specifying one to SLURM `sbatch -o out.log <script>` or by directing the output to a file (if running with `bash`), like so: `bash <script> > out.log`.

You can use these logs to debug any errors you receive running the pipeline, or email me at kdeweese@usc.edu.

