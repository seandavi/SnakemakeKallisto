# SnakemakeKallisto

## Use

On NIH Biowulf:

```bash
git clone seandavi/SnakemakeKallisto
```

Then, create a file with columns `sample`, `read1`, and `read2` that contains the final
sample names mapped to relative (or absolute) paths to fastq files (can be gzip compressed)


```bash
module load python/3.4.3
sbatch submit.sh # needs to linked into the current directory, I think
```
