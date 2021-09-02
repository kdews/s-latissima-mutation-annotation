# Mutation annotation of the brown macroalgae *Saccharina latissima* (North American sugar kelp)
A pipeline to detect and annotate the effects of deleterious mutations in *S. latissima* whole genome sequencing data.

## Installation
Clone repository:
```
git clone https://github.com/kellywithsword/s-latissima-mutation-annotation.git
```

## Dependencies
* Java (>=v1.8)
* [SnpEff](https://pcingola.github.io/SnpEff/) (>=v5.0c)
* [AGAT](https://github.com/NBISweden/AGAT) (v0.8.0)

#### Easy mode: Create Anaconda environment from YAML, and download latest SnpEff release as directory in working directory
1. If you haven't already, install the lastest version of [Anaconda](https://www.anaconda.com/)
> Note: This pipeline will attempt to source Anaconda from the `$conda_sh` environment variable set in `mut_annot.config`. **Please set** `$conda_sh` to `path/to/<anaconda-version>/etc/profile.d/conda.sh` in `mut_annot.config`. If you are installing Anaconda for the first time, now would be a good time to find this path.
2. Create Anaconda env from `mut_annot.yml`
```
conda env create -f mut_annot.yml
```
3. Download SnpEff to your working directory
```
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```

> @ Kelly: remember to add final mut\_annot.yml at end lol~

#### Hard mode: Install programs manually and add to $PATH, and download latest SnpEff release as directory in working directory
1. Install each dependency manually, and ensure that all programs (except SnpEff) have been added to your $PATH, e.g., `which java` produces `/path/to/java`. If you are not comfortable doing this step, I suggest **Easy mode**^^.

2. Download SnpEff to your working directory
```
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```

## Configure pipeline
Edit `mut_annot.config` with paths to:
* Anaconda conda.sh file (e.g., `path/to/<anaconda-version>/etc/profile.d/conda.sh`)
* Reference genome
* GFF3
* CDS FASTA
* Protein FASTA
* VCF file

Each step of the pipeline will first take the path to a config file as the first positional argument ($1); if one is not provided, it will then look for `mut_annot.config` in your current directory. 

## Run pipeline
### 1. Build SnpEff database for *S. latissima*
Run:
```
sbatch build_SnpEff_db.sbatch </path/to/mut_annot.config>
```

### 2. Run SnpEff on *S. latissima* VCF
Run:
```
sbatch predict_SnpEff.sbatch </path/to/mut_annot.config>
```

## Debugging & Troubleshooting
Log files are created at each step of the pipeline, named for the script or command used to generate them. You can use these to debug any errors you receive running the pipeline, or email me at kdeweese@usc.edu.

