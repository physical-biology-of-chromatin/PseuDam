import numpy as np
import pandas as pd
import re
import random
import matplotlib.pyplot as plt


file_signal = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3/abundance.tsv"

file_noise = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4/abundance.tsv"

file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/test/yes/sites.bed"

chrom_list = [
    "ref|NC_001133|",
    "ref|NC_001134|",
    "ref|NC_001135|",
    "ref|NC_001136|",
    "ref|NC_001137|",
    "ref|NC_001138|",
    "ref|NC_001139|",
    "ref|NC_001140|",
    "ref|NC_001141|",
    "ref|NC_001142|",
    "ref|NC_001143|",
    "ref|NC_001144|",
    "ref|NC_001145|",
    "ref|NC_001146|",
    "ref|NC_001147|",
]

chrom_names = [
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "13",
    "14",
    "15",
    "16",
]

bootstraps_number = 10000


# Reformating the dataframe 
df_pseudo = pd.read_csv(file_signal,
                            sep = "\t",
                            header = 0)


df_sites = pd.read_csv(file_sites,
                    sep = "\t",
                    header = None)
df_pseudo.insert(0,
                 column = "chrom",
                 value = df_sites[0])
df_pseudo.insert(1,
                 column = "start",
                 value = df_sites[1])
df_pseudo.insert(2,
                 column = "stop",
                 value = df_sites[2])

# Removing the chromosome 12
df_pseudo = df_pseudo[df_pseudo["chrom"] != "ref|NC_001144|"]





for chrom, name in zip(chrom_list, chrom_names):
    
    df_chrom = df_pseudo[df_pseudo["chrom"] == chrom]
    
    # Initial sampling of the data
    sample_size = round(len(df_pseudo)*0.90)

    index_list = range(0, len(df_pseudo))

    indexes_sample = random.sample(index_list, sample_size)

    df_sampled = df_pseudo.iloc[indexes_sample]

    sample_reads = list(df_sampled["est_counts"])

    # Function applied

    sample_mean = np.mean(sample_reads)

    # Bootstrapping

    bootstrap_means_list = list()

    for i in range(bootstraps_number):
        
        bootstrap_sample = random.choices(sample_reads, k = sample_size)

        bootstrap_mean = np.median(bootstrap_sample)
        
        bootstrap_means_list.append(bootstrap_mean)
        
    plt.figure(name)
    plt.hist(bootstrap_means_list,
            bins = 50)
    plt.title(name)

plt.show()