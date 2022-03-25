from more_itertools import sample
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np


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

    df_chrom_13 = df_pseudo[df_pseudo["chrom"] == "ref|NC_001144|"] 
    df_chrom_13 = df_chrom_13[df_chrom_13["start"] > 440000 ]
    df_chrom_13 = df_chrom_13[df_chrom_13["stop"] < 480000 ]

    ax[row,col].bar(((df_chrom_13["stop"] - df_chrom_13["start"])/2) + df_chrom_13["start"],
                      height = df_chrom_13["est_counts"],
                      width = df_chrom_13["length"])
    
    
    sample_name = re.search("N_(.+?)_", file_pseudo).group(1)
    
    ax[row,col].set_title(f"{sample_name}\n"
                          f"fraction of the reads in the region : {round((df_chrom_13['est_counts'].sum() / df_pseudo['est_counts'].sum())*100)}%")
    ax[row,col].set_ylabel("reads")
    ax[row,col].set_xlabel("position on chromosome 13")
    
    col += 1
    
    if col == 3:
        col = 0
        row += 1
    
    if col == 1 and row == 1:
        col = 2
        
    
    elif col == 2 and row == 3:
        col = 0
        row = 4
    
    print(df_pseudo['est_counts'].sum())
    print(df_chrom_13["est_counts"].sum())
    print("")
    
    ratio_pseudo_list.append((df_chrom_13["est_counts"].sum() / df_pseudo['est_counts'].sum()) * 100 )
    names_pseudo_list.append(sample_name)






fig.delaxes(ax[1,1])
fig.delaxes(ax[3,2])


path = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/chrom_13/all_samples.eps"
plt.savefig(path, format = "eps")

plt.close()
plt.clf()


i = 0
width = 0.35

for row in df_mapping.itertuples():
    
    ratio_mapping = (row[5] / (row[2] + row[3] + row[4] + row[5] + row[6] + row[7])) * 100
    ratio_mapping_list.append(ratio_mapping)
    names_mapping_list.append(re.search("N_(.+?)_", row[1]).group(1))
    
    plt.bar(i - width/2,
            height = ratio_mapping,
            width=width,
            color = "blue")
    
    i += 1

plt.xticks(range(len(df_mapping)), names_mapping_list)


names_list = list()
ratio_list = list()

for name, ratio_pseudo in zip(names_pseudo_list, ratio_pseudo_list):
    plt.bar(names_mapping_list.index(name) + width/2,
            height = ratio_pseudo,
            color = "red",
            width = width)

plt.ylabel("%reads in the region  /  %multimapped reads")
    
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


