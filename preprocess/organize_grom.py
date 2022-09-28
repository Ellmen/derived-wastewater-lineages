from os import listdir, mkdir, system
from os.path import isfile, join

from shutil import copy

"""
For all gzipped fastq files in a specified directory:
Pair the reads, align them to a reference, index the bam.
Requires cutadapt, minimap2, samtools, and a reference genome (called sars-cov-2.fasta)
"""

data_path = '...' # Fill in
samples_path = data_path

fns = [f for f in listdir(data_path) if isfile(join(data_path, f)) and f.endswith('.fastq.gz')]
fns.sort()
samples = []

# Relies on pairs being alphabetically adjacent
pairs = []
for i in range(int(len(fns) / 2)):
    sample_path = join(samples_path, 'sample{}'.format(i+1))
    read1 = fns[i*2]
    read2 = fns[i*2+1]
    name = read1[:read1.find('_')]
    samples.append(['sample{}/grom.mapped.csv'.format(i+1), name])
    try:
        mkdir(sample_path)
    except:
        print('{} already exists'.format(sample_path))
    cmd = 'python3 ./grom_map/minimap2.py {0}/{1} {0}/{2} --outdir {3} --prefix grom --ref ./grom_map/data/NC_045512.fa --replace'.format(samples_path, read1, read2, sample_path)
    system(cmd)

with open('{}/grom_samples.txt'.format(samples_path), 'w') as f:
    f.write('\n'.join(['\t'.join(line) for line in samples]))
