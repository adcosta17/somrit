import pysam
import argparse
import sys
import csv
import os
import math
from os import listdir
from os.path import isfile, join
from intervaltree import Interval, IntervalTree

from collections import defaultdict

parser = argparse.ArgumentParser( description='Get read subset fastq from tsv')
parser.add_argument('--reads', required=True)
parser.add_argument('--bam', required=True)
parser.add_argument('--fastq-folder', type=str, required=True)
args = parser.parse_args()

reads_to_use = {}
with open(args.reads, 'r') as in_reads:
    for line in in_reads:
        row = line.strip().split('\t')
        reads_to_use[row[3]] = 1

sam_reader = pysam.AlignmentFile(args.bam)
for record in sam_reader.fetch():
    reads_to_use[record.query_name] = 1

onlyfiles = [f for f in listdir(args.fastq_folder) if isfile(join(args.fastq_folder, f)) and f.endswith("fastq.gz")]
for file in onlyfiles:
    with pysam.FastaFile(filename = args.fastq_folder+ "/" + file, filepath_index_compressed = args.fastq_folder+ "/" +file + ".gzi") as in_fq:
        print(file, file=sys.stderr)
        for read in reads_to_use:
            try:
                seq = in_fq.fetch(read)
                print("@"+read)
                print(seq)
                print('+')
                print('='*len(seq))
            except KeyError:
                print(read+" not found", file=sys.stderr)
                pass

