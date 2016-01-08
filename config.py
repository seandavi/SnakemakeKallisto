import pprint
import csv
from itertools import chain

config = {}
config['samplesheet'] = 'samplesheet.txt'

def readSampleSheet(fname=config['samplesheet']):
    res = []
    with open(fname) as f:
        reader = csv.DictReader(f,delimiter='\t')
        for line in reader:
            res.append(line)
    return(res)




studyresults = readSampleSheet() # makeSampleSheet()


#for rec in studyresults:
#    print(rec)

def sample2fastq(wc):
    sample1 = wc.base
    print(wc.base)
    # for rna-seq alignment with STAR, need all fastqfiles for the source.  
    # Turns out these are all paired-end, so form read1 and read2 lists for return
    sampleres = list(filter(lambda res: res['sample']==sample1,studyresults))
    flatlist = chain.from_iterable([[res['read1'],res['read2']] for res in sampleres])
    return(list(flatlist))
