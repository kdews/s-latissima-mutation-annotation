# Mutation annotation of the brown macroalgae *Saccharina latissima* (North American sugar kelp)
A pipeline to detect and annotate the effects of deleterious mutations in *S. latissima* whole genome sequencing data.

## The pipeline
### 1. Build SnpEff database for *S. latissima*

Edit `build_SnpEff_db.config` with PATHs to:
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
