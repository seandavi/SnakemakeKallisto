import os
import logging

TRANSCRIPT_FASTA="Mus_musculus.GRCm38.cdna.all.fa"

FASTQDIR = "fastq/"
fnames = os.listdir(FASTQDIR)
samples = set([x[0:11] for x in fnames])
fnames = [FASTQDIR + f for f in fnames]

rule final:
    input: expand("results/{sample}/abundance.txt",sample=samples)

rule trim:
    input: "{base}_R1.fastq.gz", "{base}_R2.fastq.gz"
    output: "{base}_R1.trimmed.fastq.gz", "{base}_R2.trimmed.fastq.gz"
    threads: 1
    params:
        partition = "ccr",
        time = "8:00:00"
    shell: """
module load python
cutadapt -q 10 --minimum-length 25 --trim-n \
  -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
  -o {output[0]} -p {output[1]} \
  {input}
"""

rule kallisto_index:
    input: TRANSCRIPT_FASTA
    output: "kallisto.index"
    params:
        partition = "quick",
        mem = "8g",
        time = "1:00:00"
    shell:"""
module load kallisto
kallisto index -i kallisto.index {input}
"""
        
rule kallisto_quant:
    input: "fastq/{base}_R1.trimmed.fastq.gz", "fastq/{base}_R2.trimmed.fastq.gz", "kallisto.index"
    output: "results/{base}/abundance.tsv"
    threads: 32
    params:
        partition = "quick",
        mem = "8g",
        time = "1:00:00",
        outdir = "results/{base}"
    shell: """
module load kallisto
kallisto quant --bias -o {params.outdir} -i {input[2]} -t {threads} {input[0]} {input[1]}
"""
