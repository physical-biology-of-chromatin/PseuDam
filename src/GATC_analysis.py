import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import math

df = pd.read_csv("/home/nathan/projects/vscode_nextflow/nextflow-nathan/results/GATC/sites_yeast.bed", 
                 header = None, 
                 sep = '\t')

regions_chro = list()
chromosomes = list()
i = 0
sites_region = list()
j = 0
chrom_regions = list()
chrom_length = list()
chrom = str()
chromosome_sites = list()
sites = list()
id_list = list()

for rec in SeqIO.parse("/home/nathan/projects/vscode_nextflow/nextflow-nathan/data/genome/GCF_000146045.2_R64_genomic.fna",
                       "fasta"):
    seq = rec.seq
    chrom_length.append(len(seq))
    
    name = rec.description
    id_list.append(name[43 : -19])

for row in df.itertuples():
    if len(chrom) <= 0:
        chrom = row[1]
        
    if row[1] != chrom:
        chromosome_sites.append(sites)
        sites = list()
        chrom = row[1]
    
    sites.append(row[3])        

chromosome_sites.append(sites)

for sites in chromosome_sites:
    
    region = np.linspace(1, max(sites) + 100, round((max(sites) + 100) / 10000))
    chrom_regions.append(region)
    
    for site in sites:
        
        if site > region[i]:
            i += 1
            sites_region.append(j)
            j = 0
        j += 1
    i = 0
    sites_region.append(j)
    chromosomes.append(sites_region)
    sites_region = list()
    

fig, axes = plt.subplots(math.ceil(len(chromosomes) / 5), 5 )

i = 0
j = 0

for chrom, regions, name in zip(chromosomes, chrom_regions, id_list):
    
    if j >= 5:
        j = 0
        i += 1
        
    axes[i, j].set_title(name)
    axes[i, j].set_ylabel("site number / bin")
    axes[i, j].plot(regions, chrom)
    j += 1
    
plt.tight_layout(pad = 4.5)
plt.show()