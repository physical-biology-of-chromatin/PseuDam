import random
import re
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def bootstrap_percentile(df, percentile, b):
    
    bootstrap_means_list = list()
    
    if percentile >= 1:
        percentile = percentile / 100
    
    
    percentile_value = round((len(df) * (1 - percentile )/ 2))
    
    df = df.drop(df["est_counts"].nsmallest(percentile_value).index)
    
    df = df.drop(df["est_counts"].nlargest(percentile_value).index)
    
    sample_size = len(df)
    
    sample_reads = list(df["est_counts"])
    
    for i in range(b):
       
        bootstrap_sample = random.choices(sample_reads, k = sample_size)

        bootstrap_mean = np.mean(bootstrap_sample)
        
        bootstrap_means_list.append(bootstrap_mean)
        
    return(bootstrap_means_list)


def bootstrap_basic(df, b):
    
    df = df[df["est_counts"] >= 1]
    
    bootstraps_frags_dic = dict()
    
    data = list()
    
    i = 0
    
    for row in df.itertuples():
        
        bootstraps_frags_dic[row[4]] = list()
        
        for wheight in range(round(row[8])):
            
            data.append(row[4]) 
            
    
    sample_size = len(data)
    
    
    for i in range(b):
       
        bootstrap_sample = random.choices(data, k = sample_size)

        counted_bootstrap_sample = Counter(bootstrap_sample)
        
        for key in counted_bootstrap_sample.keys():
            
            bootstraps_frags_dic[key].append(counted_bootstrap_sample[key])
        

    return(bootstraps_frags_dic)


    
def bootstrap_analysis(dic, df):
    
    df["bootstrap_mean"] = np.NaN
    df["bootstrap_std"] = np.NaN
    df["bootstrap_var"] = np.NaN
    
    i = 0
    
    for key in dic:
        
        mean = np.mean(dic[key])
        std = np.std(dic[key])
        var = pow(std, 2)

        df.iloc[i,8] = mean 
        df.iloc[i,9] = std 
        df.iloc[i,10] = var
        
        i += 1
    
    df = df[(df["bootstrap_mean"] > 0)] 
    df = df[(df["bootstrap_var"] > 0)]
    
    return df
        


files_list = ["/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4/abundance.tsv",
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

bootstraps_number = 1000


fig, ax = plt.subplots(5, 3, 
                       sharex= True,
                       sharey= True,
                       figsize = (40, 40))

row = 0
col = 0

for file in files_list:


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

    df_pseudo.iloc[0,3] = "_0"


    # Removing the chromosome 12 and 2
    df_pseudo = df_pseudo[df_pseudo["chrom"] != "ref|NC_001144|"]

    print(df_pseudo["est_counts"].sum())
    
    print(len(df_pseudo[df_pseudo["chrom"] == "ref|NC_001135|"]) / len(df_pseudo) )
    
    print(df_pseudo.nlargest(10, columns = "est_counts"))
    
    df_pseudo = df_pseudo[df_pseudo["chrom"] != "ref|NC_001135|"]

    
    print(df_pseudo["est_counts"].sum())
    
    
    # Initial sampling of the data
    sample_size = len(df_pseudo)

    bootstrap_dic = bootstrap_basic(df_pseudo, bootstraps_number)

    df_bootstrap = bootstrap_analysis(bootstrap_dic, df_pseudo[df_pseudo["est_counts"] >= 1])

    # Associating a color scale to the density of the data
    x = np.array(df_bootstrap["bootstrap_mean"])
    y = np.array(df_bootstrap["bootstrap_var"])
    
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    
    ax[row,col].scatter(x,
                        y,
                        c=z,
                        s=50)

    ax[row,col].set_xlabel("mean")
    ax[row,col].set_ylabel("variance")
    
    sample_name = re.search("N_(.+?)_", file).group(1)
    
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
    
    print(sample_name)
    
fig.suptitle("scatterplots de la variance et de la moyenne des fragments GATC")

plt.show()
