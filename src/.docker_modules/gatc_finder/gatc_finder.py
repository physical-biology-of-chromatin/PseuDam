#!/usr/bin/env python

import re
from Bio import SeqIO
import argparse


#defines qrguments in the command line
parser = argparse.ArgumentParser()
parser.add_argument("--genome", action="store");
args = parser.parse_args();

# Gets the arguments in the command line
genome_file = args.genome


# Opening the file to write the positions (bed and GFF)
f_1 = open("sites.bed", "w")
f_2 = open("sites.gff", "w")


# Motif we are looking for
motif = "GATC"

i = 0
site = list()
sites_list = list()

# Cycles through the parsed chromosomes from the fasta file
for seq_record in SeqIO.parse(genome_file, "fasta"):
    
    # Gets the id of the chormosome in the file
    chrom = seq_record.id
    
    # Cycle throught all the motif that are found in the chromosome
    for match in re.finditer(motif, str(seq_record.seq)):
        
        i += 1
        
        start_pos = match.start() +1
        end_pos = match.end() + 1
        
        site.append(chrom)
        site.append(start_pos)
        site.append(end_pos)
        
        sites_list.append(site)
        
        site = list()
        

for i in range(2, len(sites_list)):
    
    if sites_list[i-1][0] == sites_list[i][0]:

        bin_chrom = (sites_list[i-1][0])
        bin_start = (sites_list[i-1][2] - 1)
        bin_end = (sites_list[i][1] + 1)
        
        
        # Writes the position in the .bed file (chro/start/end)
        line_f_1 = f"{bin_chrom}\t{bin_start}\t{bin_end}\n"
        f_1.write(line_f_1)
        
        
        # Writes the same thing but in gff format
        line_f_2 = f"{chrom}\tNC\texon\t{start_pos}\t{end_pos}\t.\t.\t.\t\n"
        f_2.write(line_f_2)
    
f_1.close()