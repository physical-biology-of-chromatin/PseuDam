import re
from Bio import SeqIO
import argparse
import sys



#defines qrguments in the command line
parser = argparse.ArgumentParser()
parser.add_argument("--genome", action="store");
args = parser.parse_args();

# Gets the arguments in the command line
genome_file = args.genome



# Opening the file to write the positions in
f = open("/datas/nathan/vscode_nextflow/nextflow-nathan/results/GATC/sites_yeast_new.bed", "w")

# Motif we are looking for
motif = "GATC"

# Cycles through the parsed chromosomes from the fasta file
for seq_record in SeqIO.parse(genome_file, "fasta"):
    
    # Gets the id of the chormosome in the file
    chrom = seq_record.id
    
    # Cycle throught all the motif that are found in the chromosome
    for match in re.finditer(motif, str(seq_record.seq)):
        
        start_pos = match.start() +1
        end_pos = match.end() + 1
        
        # Writes the position in the .bed file (chro/start/end)
        line = f"{chrom}\t{start_pos}\t{end_pos}\n"
        f.write(line)
    
f.close()