import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import math



df = pd.read_csv("/datas/nathan/vscode_nextflow/nextflow-nathan/results/GATC/sites_yeast_new.bed", 
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
number_sites = 1/(gc * gc * at * at)

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
    
    distance = 0
    previous_site = 0
    
    regions = list()
    
    for site in sites:
        
        if site > region[i]:
            
            if len(sites_region) == 0:
                sites_region.append(0)
            
            regions.append(np.mean(sites_region))
            i += 1
            
        
        distance = site - previous_site
        
        previous_site = site

        if previous_site != 0:
            
           sites_region.append(distance)
    
    i = 0
    regions.append(np.mean(sites_region))
    chromosomes.append(regions)
    
    
    sites_region = list()

i = 0
j = 0

for chrom, regions, name in zip(chromosomes, chrom_regions, range(1,len(chrom_regions)+1)):
    
    
    if j >= 5:
        j = 0
        i += 1
    
    pos = np.arange(1, int(max(regions)), 1 )
    y = np.full(len(pos), number_sites)
    
    plt.title(f"Chrom {name}")
    plt.ylabel("site number / bin")
    plt.plot(pos, y, color = "black")
    plt.bar(regions, chrom, width = region_size, edgecolor = "black")

    j += 1

    mean_sites = round(np.mean(chrom))
    
    means_list.append(mean_sites)

    
    print(f"The mean distance between sites is {mean_sites} "
          f"in bins of {region_size}b in the choromosome {name+1}")

    plt.suptitle(f"number of GATC sites per {region_size} base")
    
    plt.savefig("/datas/nathan/vscode_nextflow/nextflow-nathan"
                "/results/sites analysis/distances_chromosome_"
                + str(name)
                + "_region_"
                + str(region_size / 1000)
                + "kb"
                + ".png")
    plt.clf()

mean_total = round(np.mean(means_list))
print(f"The total mean distance between sites is {mean_total}")
print(f"La distance attendue entre les sites est de {round(number_sites)}")
