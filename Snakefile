import os
import logging

GENOMEFASTA="Mus_musculus.GRCm38.dna.toplevel.fa"
GTF = "Mus_musculus.GRCm38.82.gtf"

FASTQDIR = "fastq/"
fnames = os.listdir(FASTQDIR)
samples = set([x[0:11] for x in fnames])
fnames = [FASTQDIR + f for f in fnames]

BAMS = ["bam/" + s + ".Aligned.sortedByCoord.out.bam" for s in samples]
logging.info('input filenames:')
for f in fnames:
    logging.info('   '+f)
logging.info(BAMS)
    
rule final:
    input: BAMS

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
    threads: 32
    params:
        partition = "ccr",
        mem = "8g",
        time = "4:00:00",
    shell:"""
module load kallisto
kallisto index -i kallisto.index {input}
"""
        
rule kallisto_quant:
    input: "fastq/{base}_R1.trimmed.fastq.gz", "fastq/{base}_R2.trimmed.fastq.gz", "kallisto.index"
    output: "bam/{base}/abundance.txt"
    threads: 32
    params:
        partition = "ccr",
        mem = "8g",
        time = "4:00:00"
    shell: """
module load kallisto
kallisto quant --bias -i {input[2]} -t {threads} {input[0]} {input[1]}
"""
