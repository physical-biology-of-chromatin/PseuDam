from operator import truediv
import pandas as pd
import matplotlib.pyplot as plt
import re
import random
import numpy as np

files = [
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4/abundance.tsv"
         ]

"""
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1D2_/abundance.tsv",
"""



"""
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_2RG1_/abundance.tsv",
         "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_2RG2_EKDL220001502-1a_H2C2YDSX3_L3/abundance.tsv",
"""

bin_size = 500

count_mean = list()

mapping_results = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/bowtie2_pe_plot.tsv"

file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/test/yes/sites.bed"

df_mapping = pd.read_csv(mapping_results,
                         sep = "\t",
                         header = 0)

region_start = 450000
region_stop = 470000

window = 40000

draw_number = 10000

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
    "ref|NC_001145|",
    "ref|NC_001146|",
    "ref|NC_001147|",
    "ref|NC_001148|",
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

ratio_file_list = list()

col = 0
row = 0

ratio_mean = list()

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
    df_pseudo = df_pseudo[df_pseudo["chrom"] != "ref|NC_001135|"]

    #df_pseudo = df_pseudo[df_pseudo["est_counts"] > 1]

    df_pseudo.reset_index(drop = True,
                          inplace=True)

    print(df_pseudo["est_counts"].mean())

    ratio_file = df_pseudo["est_counts"].mean() / df_pseudo["est_counts"].sum()

    print(ratio_file)

    ratio_file_list.append(ratio_file)


z_factor = ratio_file_list[0] / ratio_file_list[1]



df_signal = pd.read_csv(files[0],
                        sep = "\t",
                        header = 0)

df_noise = pd.read_csv(files[1],
                        sep = "\t",
                        header = 0)

df_sites = pd.read_csv(file_sites,
                    sep = "\t",
                    header = None)


df_signal.insert(0,
                column = "chrom",
                value = df_sites[0])

df_noise.insert(0,
                column = "chrom",
                value = df_sites[0])

df_signal.insert(1,
                column = "start",
                value = df_sites[1])

df_noise.insert(1,
                column = "start",
                value = df_sites[1])

df_signal.insert(2,
                column = "stop",
                value = df_sites[2])

df_noise.insert(2,
                column = "stop",
                value = df_sites[2])



df_noise = df_noise[df_noise["chrom"] == "ref|NC_001135|"]
df_signal = df_signal[df_signal["chrom"] == "ref|NC_001135|"]

df_noise = df_noise[df_noise["start"] > 92000 - window/2]
df_noise = df_noise[df_noise["stop"] < 92000 + window/2]

df_signal = df_signal[df_signal["start"] > 92000 - window/2]
df_signal = df_signal[df_signal["stop"] < 92000 + window/2]



bins = np.arange((92000 - window/2) - bin_size,
                 (92000 + window/2) + bin_size,
                 bin_size)

reads_bins_signal = list()
reads_bins_noise = list()

for bin in bins:
    df_bin_signal = df_signal[df_signal["stop"] > bin]
    df_bin_signal = df_signal[df_signal["start"] < bin + bin_size]
    
    reads_bins_signal.append(int(df_bin_signal["est_counts"].sum()) / len(df_bin_signal))
    
    
    df_bin_noise = df_noise[df_noise["stop"] > bin]
    df_bin_noise = df_noise[df_noise["start"] < bin + bin_size]
    
    reads_bins_noise.append(int(df_bin_noise["est_counts"].sum()) / len(df_bin_noise))

reads_ratio_bins = list(map(truediv, reads_bins_signal, reads_bins_noise))

print(z_factor)

middle = ((df_noise["stop"] - df_noise["start"])/2) + df_noise["start"]
ratio = df_signal["est_counts"] / df_noise["est_counts"]
width = df_noise["length"]

plt.figure("no_z")
plt.bar(middle,
        height=np.log2(ratio),
        width=width,
        color = ["red" if i < 0 else "blue" for i in np.log2(ratio)])

df_noise["est_counts"] = df_noise["est_counts"] / z_factor
df_signal["est_counts"] == df_signal["est_counts"] * z_factor

ratio = df_signal["est_counts"] / df_noise["est_counts"]

plt.figure("z")
plt.bar(middle,
        height=np.log2(ratio),
        width=width,
        color = ["red" if i < 0 else "blue" for i in np.log2(ratio)])











