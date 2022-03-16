from operator import truediv
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


files_list = [
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4_cov.bed"
]


bin_sizes = [500]


representation_list = ["centered"]

window_list = [100000]

files_values = list()
for file in files_list:
    
    for representation, window_size in zip(representation_list, window_list):

        for bin_size in bin_sizes:

            
            df_dam = pd.read_csv(file, 
                                header = None, 
                                sep = '\t')



            chrom_name = "ref|NC_001135|"

            # Selection of the rows
            # Chosen chormosome / in the window
            
            if representation == "whole":
                
                df_currated = df_dam[df_dam[0] == chrom_name] 

            else:   
                df_currated = df_dam[df_dam[0] == chrom_name] 
                

                
                df_currated = df_currated[df_currated[1] > 92500 - (window_size/2)]
                df_currated = df_currated[df_currated[2] < 92500 + (window_size/2)] 

            
            
            # Addition of a column contianing the reads normalised by the size
            df_currated.insert(loc = 7,
                               column = "7",
                               value = df_currated[3] / df_currated[5])

            if representation == "whole":
                
                # Creation of the bins
                regions = np.arange(0,
                                    df_currated[2].max()+bin_size,
                                    bin_size)
                
            else:
                regions = np.arange(92500 - (window_size/2) - bin_size,
                                    92500 + (window_size/2) + bin_size,
                                    bin_size)

            

            i = 1
            value_region = 0

            regions_value_list = list()
            regions_pos_list = list()

            gatc = 0
            frag_number = 0
            
            
            # Transfert des reads des fragments dans les bins
            
            for row in df_currated.itertuples():
                
                frag_number += 1
                frag_width = row[3] - row[2]
                frag_reads = row[8]
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
    files_values.append(regions_value_list)



files_values_ratio = list(map(truediv, files_values[1], files_values[0]))
                          
            
plt.bar(regions_pos_list,
        height = files_values_ratio, 
        width = bin_size)
plt.ylabel("ratio")
plt.xlabel("positions")

plt.title("ratio of 3RG1D3 / 3RG1")
plt.suptitle(f"Mean number of reads per {bin_size}")
path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
        + "/norm/"
        + "analysis_norm_ratio"
        + "_"
        + representation
        + "_depth_bins_"
        + str(bin_size)
        + ".eps")
plt.savefig(path, format = "eps")
print(path)

