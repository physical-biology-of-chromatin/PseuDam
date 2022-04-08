import pandas as pd
import matplotlib.pyplot as plt
import re
import random

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

region_start = 450000
region_stop = 470000

window = 20000

draw_number = 500

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
fig, ax = plt.subplots(5, 3, 
                       sharex= True,
                       sharey= True,
                       figsize = (40, 40))

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

    df_12 = df_pseudo[df_pseudo["chrom"] == "ref|NC_001144|"]
    
    df_12 = df_12[df_12["start"] > region_start]
    df_12 = df_12[df_12["stop"] < region_stop]
    
    count_12 = df_12["est_counts"].sum() / 150

    ratio_draws = list()
    count_draws = list()
    regions_draws = list()
    
    regions_chrom = list()
    count_chrom = list()
    ratio_chrom = list()
    
    for chrom in chrom_list:
        
        df_chrom = df_pseudo[df_pseudo["chrom"] == chrom]
        
        for draw in range(draw_number):

            position = random.randint(window/2, int(df_pseudo["stop"].max()) - window/2)
            
            region_draw = [position - window/2, position + window/2]

            df_temp = df_chrom[df_chrom["start"] > region_draw[0]]
            df_temp = df_temp[df_temp["stop"] < region_draw[1]]
            
            count_region = df_temp["est_counts"].sum()
            
            ratio_region = count_region / count_12
            
            ratio_draws.append(ratio_region)
            count_draws.append(count_region)
            regions_draws.append(region_draw)
        
        regions_chrom.append(regions_draws)
        count_chrom.append(count_draws)
        ratio_chrom.append(ratio_draws)
    
        ratio_draws = list()
        count_draws = list()
        regions_draws = list()
        
        
    ax[row,col].violinplot(ratio_chrom)
        
    sample_name = re.search("N_(.+?)_", file_pseudo).group(1)
    
    ax[row,col].set_xticks(range(1,len(ratio_chrom)+1))
    ax[row,col].set_xticklabels(chrom_names)
    ax[row,col].set_xlabel("chromosomes compared to the region")
    ax[row,col].set_ylabel("ratio of the total reads")
    
    ax[row,col].set_title(f"{sample_name}\n"
                          f"{round((df_12['est_counts'].sum() / df_pseudo['est_counts'].sum())*100)}% of the reads in region")
    col += 1
    
    if col == 3:
        col = 0
        row += 1
    
    if col == 1 and row == 1:
        col = 2
        
    
    elif col == 2 and row == 3:
        col = 0
        row = 4


"""plt.tight_layout()"""

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


