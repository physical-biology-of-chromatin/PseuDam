import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


files_list = [
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_2RG1_EKDL220001501-1a_H2LFNDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_2RG2_EKDL220001502-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG1D2_EKDL220001507-1a_H25VJDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3_cov.bed"]



bin_sizes = [100, 200, 500, 1000, 2000, 3000, 5000]


for file in files_list:
    
    for bin_size in bin_sizes:
        
        df_dam = pd.read_csv(file, 
                            header = None, 
                            sep = '\t')


        window_size = 80000
        chrom_name = "ref|NC_001135|"

        # Selection of the rows
        # Chosen chormosome / in the window
        df_currated = df_dam[df_dam[0] == chrom_name] 

        df_currated = df_currated[df_currated[1] > 90000 - (window_size/2)]
        df_currated = df_currated[df_currated[2] < 90000 + (window_size/2)] 
   

        
        # Addition of a column contianing the reads normalised by the size
        df_currated.insert(loc = 7,
                           column = "7",
                           value = df_currated[3] / df_currated[5])


        
        # Creation of the bins
        regions = np.arange(90000 - (window_size/2) - bin_size,
                            90000 + (window_size/2)+ bin_size,
                            bin_size)

        
        fig, ax = plt.subplots()

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
            
            if i == 1:
                i = 1
                if regions[i-1] > row[3]:
                    while regions[i-1] < row[3]:
                        
                        middle_region = (((regions[i] - regions[i-1]) / 2) + regions[i - 1])
                        regions_pos_list.append(middle_region)
                        
                        regions_value_list.append(0)
                        i += 1
                    
            
            if (row[3] < regions[i]
                and row[2] > regions[i-1]):
                
                value_region += frag_reads
                

            elif row[3] > regions[i]:
                
                value_region += frag_reads
                

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



        sample_name = re.search("N_(.+?)_EKDL", file).group(1)
        
        ax.bar(regions_pos_list,
                height = regions_value_list, 
                width = bin_size)

        ax.set_ylabel("reads / frag_length")
        ax.set_xlabel("positions")
        
        ax.set_title(sample_name)

        fig.suptitle(f"Mean number of reads per {bin_size} normalised by the fragment size")


        path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
                + sample_name
                + "/"
                + "analysis_focused_"
                + sample_name
                + "_depth_bins_"
                + str(bin_size))

        plt.savefig(path)
        plt.clf()
        plt.close()
