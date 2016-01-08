import os
import logging

TRANSCRIPT_FASTA="Mus_musculus.GRCm38.cdna.all.fa"

import config

samples = set(x['sample'] for x in config.studyresults)

rule final:
    input: expand("results/{sample}/abundance.tsv",sample=samples)

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
    threads: 1
    params:
        partition = "quick",
        mem = "32000",
        time = "1:30:00",
        version = '0.42.4'
    shell:"""
module load kallisto/{params.version}
kallisto index -i kallisto.index {input}
"""
        
rule kallisto_quant:
    input: fastq=config.sample2fastq, index="kallisto.index"
    output: "results/{base}/abundance.tsv"
    threads: 32
    params:
        partition = "quick",
        mem = "8g",
        time = "1:30:00",
        outdir = "results/{base}",
        version = "0.42.4"
    shell: """module load kallisto/{params.version}
kallisto quant --bias -o {params.outdir} -i {input.index} -t {threads} {input.fastq}
"""
