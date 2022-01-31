import re

import matplotlib.pyplot as plt
import numpy as np
import pandas
from Bio import SeqIO, motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    
    f = open("/home/nathan/projects/vscode_nextflow/nextflow-nathan/results/GATC/sites.bed", "w")
    motif = "GATC"
    pos_list = list()

    for seq_record in SeqIO.parse("/home/nathan/projects/vscode_nextflow/nextflow-nathan/data/genome/data_G.fasta", "fasta"):
        chrom = seq_record.id
        
        for match in re.finditer(motif, str(seq_record.seq)):
            start_pos = match.start() +1
            end_pos = match.end() + 1
            
            line = f"{chrom}\t{start_pos}\t{end_pos}\n"
            
            f.write(line)


if __name__ == "__main__":
    main()