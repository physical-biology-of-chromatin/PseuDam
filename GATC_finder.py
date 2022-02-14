import re

from Bio import SeqIO


def main(genome_file, out_file_path):
    """[Gets all the GATC file from the given genome or sequence and puts them in a .bed file]

    Args:
        genome_file ([string]): [full path to the fasta file]
        out_file_path ([string]): [full path to the output file]
    """
    # Opening the file to write the positions in
    f = open(out_file_path, "w")
    
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

if __name__ == "__main__":
    main()
