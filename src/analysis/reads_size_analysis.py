from more_itertools import sample
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
from scipy.stats import gaussian_kde

files = ["/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1D2_/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_2RG1_/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_2RG2_EKDL220001502-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4/abundance.tsv"]

mapping_results = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/bowtie2_pe_plot.tsv"

file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/test/yes/sites.bed"

df_mapping = pd.read_csv(mapping_results,
                         sep = "\t",
                         header = 0)

ratio_mapping_list = list()
names_mapping_list = list()

ratio_pseudo_list = list()
names_pseudo_list = list()





fig, ax = plt.subplots(5, 3, 
                       sharex= True,
                       sharey= True,
                       figsize = (30, 30))

col = 0
row = 0

for file_pseudo in files:
    
    df_pseudo = pd.read_csv(file_pseudo,
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

    df_pseudo = df_pseudo[df_pseudo["chrom"] != "ref|NC_001144|"]
    
    df_pseudo = df_pseudo[df_pseudo["est_counts"] > 100]
    
    x = np.array(df_pseudo["length"])
    y = np.array(df_pseudo["est_counts"])
    
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    
    ax[row,col].scatter(x,
                     y,
                     c = z,
                     s=50)
    
    ax[row,col].set_ylim(0,20000)
    
    sample_name = re.search("N_(.+?)_", file_pseudo).group(1)
    
    ax[row,col].set_xlabel("fragment size")
    ax[row,col].set_ylabel("read number")
    ax[row,col].set_title(sample_name)
    col += 1
    
    if col == 3:
        col = 0
        row += 1
    
    if col == 1 and row == 1:
        col = 2
        
    
    elif col == 2 and row == 3:
        col = 0
        row = 4







fig.delaxes(ax[1,1])
fig.delaxes(ax[3,2])

plt.show()










"""
df_sizes = df_pseudo[["length", "est_counts"]]

df_sizes.sort_values(by = ["length"], inplace=True)
sizes = df_sizes["length"].unique()
data = list()
data_sizes = list()
for size in sizes:
    
    df_temp = df_sizes[df_sizes["length"] == size]
    if 1000 > df_temp["est_counts"].mean() > 1:
        data_sizes.append(size)
        data.append(df_temp["est_counts"].mean())
    
plt.plot(df_pseudo["length"],
         df_pseudo["est_counts"],
         "x")

plt.show()
"""


