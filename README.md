# Mutation annotation of the brown macroalgae *Saccharina latissima* (North American sugar kelp)
### A pipeline to detect deleterious mutations in *S. latissima* WGS data


## Rename original FASTQ files to have less metadata
In the directory containing FASTQ files, run:
```
python3 rename_MO_KD.py
```

## Optional: Annotation of contaminants in WGS reads with BBDuk and BLAST
The transcriptome and genome that we're using have been decontaminated, but to detect and annotate contaminating reads (e.g. bacteria, fungi, viruses) in *S. latissima* DNA reads, run:
```
mkdir bbduk
for i in $sample_ids; echo $i; sbatch -J $i -o bbduk/$i.bbduk.out bbduk_submit.sh $i; done
mkdir blast
for i in $sample_ids; echo $i; sbatch -J $i -o blast/$i.convert.out fastq_convert.sh $i; done
for i in $sample_ids; echo $i; sbatch -J $i -o blast/$i.blast_1.out blast_submit_1.sh $i; done
for i in $sample_ids; echo $i; sbatch -J $i -o blast/$i.blast_2.out blast_submit_2.sh $i; done
```


