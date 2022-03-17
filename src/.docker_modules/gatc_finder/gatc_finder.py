#!/usr/bin/env python

import re
from Bio import SeqIO
import argparse


#defines qrguments in the command line
parser = argparse.ArgumentParser()
parser.add_argument("--genome", action="store")
parser.add_argument("--overlap", action="store_true")
args = parser.parse_args()

# Gets the arguments in the command line
genome_file = args.genome
overlap = args.overlap

# Opening the file to write the positions (bed and GFF)
f_bed = open("sites.bed", "w")
f_gff = open("sites.gff", "w")


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
        
        # If overlap is specified, uses a bigger overlap
        if overlap == True:
            overlap_value = 10
            
        # If not sepcified takes only the whole gatc into account
        elif overlap == False:
            overlap_value = 4
        
        
        
        bin_chrom = (sites_list[i-1][0])
        bin_start = (sites_list[i-1][2] - overlap_value)
        bin_end = (sites_list[i][1] + overlap_value)
        
        
        # Writes the position in the .bed file (chro/start/end)
        line_f_gff = f"{bin_chrom}\t{bin_start}\t{bin_end}\n"
        f_bed.write(line_f_gff)
        
        # Gets the chromosome name from the file
        chrom_name = re.search("|(.+?)|", sites_list[i][0]).group(1)

        # Writes the same thing but in gff format
        line_f_gff = f"{chrom_name}\tgatc_finder\tgatc_frag\t{bin_start}\t{bin_end}\t.\t.\t.\n"
        f_gff.write(line_f_gff)
    
f_bed.close()