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

rule generateGenome:
    input: GENOMEFASTA, GTF
    output: "GenomeDir/abc"
    threads: 32
    params:
        partition = "ccr",
        mem = "58g",
        time = "24:00:00"
    shell: """
module load STAR
STAR \
    --runMode genomeGenerate \
    --runThreadN {threads} \
    --genomeFastaFiles {input[0]} \
    --sjdbGTFfile {input[1]} \
    --sjdbOverhang 75 \
    --genomeDir GenomeDir
"""    
    
rule align:
    input: "fastq/{base}_R1.trimmed.fastq.gz", "fastq/{base}_R2.trimmed.fastq.gz", "GenomeDir/abc"
    output: "bam/{base}.Aligned.sortedByCoord.out.bam"
    threads: 32
    params:
        partition = "ccr",
        mem = "58g",
        time = "24:00:00",
        outdir = "bam/{base}."
    shell: """
module load STAR
STAR \
    --runThreadN {threads} \
    --genomeDir {input[2]} \
    --readFilesIn {input[0]} {input[1]} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFileNamePrefix {params.outdir}
"""
