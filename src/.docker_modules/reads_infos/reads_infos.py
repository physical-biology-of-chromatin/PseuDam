#!/usr/bin/env python

"""
TODO finish doc
Gets the mean fragment length and estimated standard deviation from a single end fastq file

Input:


Output:


"""

from Bio import SeqIO
import numpy as np

import argparse

# Defines arguments for command line call
parser = argparse.ArgumentParser()

parser.add_argument("--fastq", action="store",
                    help = "<path> Path to the fastq file containing the single end reads")
parser.add_argument("--all", action="store_true",
                    help = "<path> Does the mean of the length of all the reads (might take some time)")
parser.add_argument("--sample_size", action="store",
                    help = "<path> number of reads to do the mean from (defaults to 1,000,000)")

args = parser.parse_args()

fastq_file = args.fastq
sample_size = args.sample_size

if sample_size == None:
    sample_size = 1e6

if args.all != False:
    sample_size = len(SeqIO.parse(fastq_file, "fastq"))
    print(sample_size)
    
#fastq_file = "/datas/nathan/vscode_nextflow/nextflow-nathan/data/reads/Dam_ID/c_elegans/SRR13429076.fastq"

length_list = list()

for i, seq_record in enumerate(SeqIO.parse(fastq_file, "fastq")):
    
    length = len(seq_record.seq)
    length_list.append(length)
    
    if i == sample_size:
        break

mean_frag_length = np.mean(length_list)
std_frag_length = np.std(length_list)

print(mean_frag_length)
print(std_frag_length)
