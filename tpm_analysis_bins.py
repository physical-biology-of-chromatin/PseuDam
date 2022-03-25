from operator import truediv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file_align_signal = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3/abundance.tsv"

file_align_noise = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4/abundance.tsv"


file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/sites.bed"

df_align_signal = pd.read_csv(file_align_signal,
                       sep = "\t",
                       header = 0)

df_align_noise = pd.read_csv(file_align_noise,
                       sep = "\t",
                       header = 0)


df_sites = pd.read_csv(file_sites,
                       sep = "\t",
                       header = None)


df_align_signal.insert(0,
                column = "chrom",
                value = df_sites[0])

df_align_signal.insert(1,
                column = "start",
                value = df_sites[1])

df_align_signal.insert(2,
                column = "stop",
                value = df_sites[2])


df_align_noise.insert(0,
                column = "chrom",
                value = df_sites[0])

df_align_noise.insert(1,
                column = "start",
                value = df_sites[1])

df_align_noise.insert(2,
                column = "stop",
                value = df_sites[2])

chrom_name = "ref|NC_001135|"

df_cur_signal = df_align_signal[df_align_signal["chrom"] == chrom_name]

df_cur_noise = df_align_noise[df_align_noise["chrom"] == chrom_name]

bin_size = 500

representation = "whole"
window_size = 100000

if representation == "whole":

# Creation of the bins
    regions = np.arange(0,
                        df_cur_signal["stop"].max()+bin_size,
                        bin_size)

else:
    regions = np.arange(92500 - (window_size/2) - bin_size,
                        92500 + (window_size/2) + bin_size,
                        bin_size)



i = 1
value_region = 0

regions_value_noise_list = list()
regions_pos_list = list()

gatc = 0
frag_number = 0

norm_choice = "norm"


for row_noise in df_cur_noise.itertuples():
    
    frag_number += 1
    frag_width = row_noise[3] - row_noise[2]
    
    if norm_choice == "norm":
        frag_reads = row_noise[7] / row_noise[5]
    elif norm_choice == "no_norm":
        frag_reads = row_noise[7]
    region = regions[i]
    
    # skips the first few bins without fragments
    if i == 1:
        i = 1
        if regions[i] < row_noise[2]:
            
            while not (regions[i-1] < row_noise[2] < regions[i]):
                
                
                regions_value_noise_list.append(0)
                i += 1
            
    # fragment fully inside the bin
    if (row_noise[3] < regions[i]
        and row_noise[2] > regions[i-1]):
        
        value_region += frag_reads
        
    # fragment with one end outside
    elif row_noise[3] > regions[i]:
        
        value_region += frag_reads
        
        #fragment over mulitple bins
        if row_noise[3] > regions[i+1]:
        
            while row_noise[3] > regions[i+1]:
                
                regions_value_noise_list.append(value_region / frag_number)
                
                value_region = 0
                
                i += 1
                
                value_region += frag_reads
                
                frag_number = 1
                
        regions_value_noise_list.append(value_region / frag_number)
        
        frag_number = 1
        value_region = 0
        i += 1
                
        value_region += frag_reads


i = 1
value_region = 0
regions_value_list = list()
regions_pos_list = list()
gatc = 0
frag_number = 0



for row in df_cur_signal.itertuples():
    
    frag_number += 1
    frag_width = row[5]
    
    if norm_choice == "norm":
        frag_reads = row[7] / row[5]
        
    elif norm_choice == "no_norm":
        frag_reads = row[7]
        
    region = regions[i]
    
    # skips the first few bins without fragments
    if i == 1:
        i = 1
        if regions[i] < row[2]:
            
            while not (regions[i-1] < row[2] < regions[i]):
                
                middle_region = (((regions[i] - regions[i-1]) / 2) + regions[i - 1])
                regions_pos_list.append(middle_region)
                
                regions_value_list.append(0)
                i += 1
            
    # fragment fully inside the bin
    if (row[3] < regions[i]
        and row[2] > regions[i-1]):
        
        value_region += frag_reads
        
    # fragment with one end outside
    elif row[3] > regions[i]:
        
        value_region += frag_reads
        
        #fragment over mulitple bins
        if row[3] > regions[i+1]:
        
            while row[3] > regions[i+1]:
                
                middle_region = (((regions[i] - regions[i-1]) / 2) + regions[i - 1])
                regions_pos_list.append(middle_region)
                
                regions_value_list.append(value_region / frag_number)
                
                value_region = 0
                
                i += 1
                
                value_region += frag_reads
                
                frag_number = 1
                
        middle_region = (((regions[i] - regions[i-1]) / 2) + regions[i - 1])
        regions_pos_list.append(middle_region)
        
        regions_value_list.append(value_region / frag_number)
        
        frag_number = 1
        value_region = 0
        i += 1
                
        value_region += frag_reads

regions_value_list = [i + 1 for i in regions_value_list]
regions_value_noise_list = [i + 1 for i in regions_value_noise_list]
ratio = list(map(truediv, regions_value_list, regions_value_noise_list))

plt.bar(regions_pos_list,
        height = np.log2(ratio),
        width = bin_size,
        color = ["red" if i < 0 else "blue" for i in np.log2(ratio)])

plt.show()