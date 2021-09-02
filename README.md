# Mutation annotation of the brown macroalgae *Saccharina latissima* \
# (North American sugar kelp)
A pipeline to detect and annotate the effects of deleterious mutations in \
*S. latissima* whole genome sequencing data.

## The pipeline
### 1. Install dependencies
* Java (>=v1.8)
* [AGAT] (https://github.com/NBISweden/AGAT) (v0.8.0)
#### Easy mode: Create Anaconda environment from mut\_annot.yml

If you haven't already, install the lastest version of [Anaconda] \
(https://www.anaconda.com/).

The scripts will attempt to source Anaconda from the `$conda_sh` environment \
variable set in `mut_annot.config`. Please set `$conda_sh` to \
`path/to/<anaconda-version>/etc/profile.d/conda.sh` in `mut_annot.config`.

Create mut\_annot conda env:

```
conda env create -f mut_annot.yml
```

@ Kelly: remember to add final mut\_annot.yml at end~

### Hard mode: Install programs manually and add to $PATH

Ensure that each dependency has been added to your $PATH. If you are not \
comfortable doing this step, I suggest **Easy mode**.

### 2. Build SnpEff database for *S. latissima*

Edit `build_SnpEff_db.config` with paths to:
* Reference genome
* GFF3
* CDS FASTA
* protein FASTA

Run:
```
sbatch build_SnpEff_db.sbatch
```

### 2. Run SnpEff on *S. latissima* VCF

***Copy*** input VCF file to new `snpEff` directory, and run:
```
sbatch predict_SnpEff.sbatch
```
