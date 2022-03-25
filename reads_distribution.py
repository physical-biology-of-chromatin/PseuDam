import numpy as np
import pandas as pd
import re
import random
import matplotlib.pyplot as plt


file_signal = ["/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3/abundance.tsv",
               "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4/abundance.tsv"]




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


for file in file_signal:

    # Reformating the dataframe 
    df_pseudo = pd.read_csv(file,
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
    df_pseudo = df_pseudo[df_pseudo["est_counts"] > 1]

    plt.figure(file)

    plt.yscale("log")
    plt.xlim(-100,140000)

    plt.hist(df_pseudo["est_counts"] ,
            bins = 50)

plt.show()