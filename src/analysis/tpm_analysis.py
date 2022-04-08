import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file_pseudo = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3/abundance.tsv"

file_classical = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3_cov.bed"


file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/sites.bed"

df_pseudo = pd.read_csv(file_pseudo,
                        sep = "\t",
                        header = 0)

df_classical = pd.read_csv(file_classical,
                             sep = "\t",
                             header = None)


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

print(df_pseudo)


chrom_name = "ref|NC_001135|"

df_cur_pseudo = df_pseudo[df_pseudo["chrom"] == chrom_name]

df_cur_classical = df_classical[df_classical[0] == chrom_name]

print(df_classical[3].mean())
print(df_pseudo["est_counts"].mean())
print(np.mean(df_classical[3] - df_pseudo["est_counts"]))

plt.plot(((df_pseudo["stop"] - df_pseudo["start"])/2) + df_pseudo["start"],
         df_classical[3] - df_pseudo["est_counts"])

plt.show()

