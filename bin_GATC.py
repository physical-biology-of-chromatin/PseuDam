import pybedtools
import pysam
from Bio import SeqIO

sites = pybedtools.BedTool("/home/nathan/projects/vscode_nextflow/nextflow-nathan/results/GATC/sites_yeast.bed")



samfile = pysam.AlignmentFile("/home/nathan/projects/vscode_nextflow/nextflow-nathan/results/mapping/data.bam","rb")

print(samfile)

print(samfile)
for read in samfile.fetch("chr1", 100, 120):
    print(read)


