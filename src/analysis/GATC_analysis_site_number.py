import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import math



df = pd.read_csv("/datas/nathan/vscode_nextflow/nextflow-nathan/results/GATC/sites_yeast.bed", 
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



region_size = 20000
at = 0.3085
gc = 0.1915
number_sites = region_size * gc * gc * at * at

for rec in SeqIO.parse("/datas/nathan/vscode_nextflow/nextflow-nathan/data/genome/dm6.fasta",
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

means_list = list()

for sites in chromosome_sites:
    
    region = np.linspace(1, max(sites) + 100, 
                         math.ceil((max(sites) + 100) / region_size))
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
    



i = 0
j = 0

for chrom, regions, name in zip(chromosomes, chrom_regions, range(len(chrom_regions))):
    
    
    if j >= 5:
        j = 0
        i += 1
    
    pos = np.arange(1, int(max(regions)), 1 )
    y = np.full(len(pos), number_sites)
    
    plt.title(f"Chrom {name+1}")
    plt.ylabel("site number / bin")
    plt.plot(pos, y, color = "black")
    plt.bar(regions, chrom, width = region_size, edgecolor = "black")

    j += 1

    mean_sites = round(np.mean(chrom))
    
    means_list.append(mean_sites)

    
    print(f"There is {mean_sites} sites per {region_size}b in the choromosome {name+1}")

    plt.suptitle(f"number of GATC sites per {region_size} base")
    
    plt.savefig("/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites analysis/chromosome_"
                + str(name+1)
                + "_region_"
                + str(region_size / 1000)
                + "kb"
                + ".png")
    plt.clf()

mean_total = round(np.mean(means_list))
print(f"there are {mean_total} every {region_size}b")
print(f"the theoretical average number of sites per {region_size}b"
      f" is {round(number_sites)}")