import re
import sys
from Bio import SeqIO

if len(sys.argv) != 3:
    raise IndexError("Please enter 2 additional arguments to the function call")

f = open(sys.argv[2], "w")
motif = "GATC"
pos_list = list()

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):

    chrom = seq_record.id
    
    for match in re.finditer(motif, str(seq_record.seq)):
        start_pos = match.start() +1
        end_pos = match.end() + 1
            
        line = f"{chrom}\t{start_pos}\t{end_pos}\n"
            
        f.write(line)
