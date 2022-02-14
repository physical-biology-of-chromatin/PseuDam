import sys
import re
from Bio import SeqIO

if len(sys.argv) < 3:
    raise IndexError("Please enter 2 arguments")

# Gets the arguments in the command line
out_file_path = "/GATC_finder/sites"
genome_file = str(sys.argv[1])

# Opening the file to write the positions in
f = open(out_file_path, "w")

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