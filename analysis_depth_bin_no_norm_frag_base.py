import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


files_frag_list = [
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


files_base_list = [
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_2RG1_EKDL220001501-1a_H2LFNDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_2RG2_EKDL220001502-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_3RG1D2_EKDL220001507-1a_H25VJDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3_cov.bed"]

bin_sizes = [500]


representation_list = ["centered"]

window_list = [50000]


for file_frag, file_base in zip(files_frag_list, files_base_list):
    
    for representation, window_size in zip(representation_list, window_list):

        for bin_size in bin_sizes:
            
            
            fig, ax = plt.subplots(figsize = (20,10))
            
            df_dam_frag = pd.read_csv(file_frag, 
                                      header = None, 
                                      sep = '\t')

            df_dam_base = pd.read_csv(file_base, 
                                      header = None, 
                                      sep = '\t')

            
            


            chrom_name = "ref|NC_001135|"

            # Selection of the rows
            # Chosen chormosome / in the window
            
            if representation == "whole":
                
                df_currated = df_dam_frag[df_dam_frag[0] == chrom_name] 
                df_dam_base = df_dam_base[df_dam_base[0] == chrom_name]
            else:   
                df_currated = df_dam_frag[df_dam_frag[0] == chrom_name] 
                df_dam_base = df_dam_base[df_dam_base[0] == chrom_name]

                
                df_currated = df_currated[df_currated[1] > 92500 - (window_size/2)]
                df_currated = df_currated[df_currated[2] < 92500 + (window_size/2)]

                df_dam_base = df_dam_base[df_dam_base[1] > 92500 - (window_size/2)]
                df_dam_base = df_dam_base[df_dam_base[2] < 92500 + (window_size/2)] 



            # Transfert des reads des fragments dans les bins

            sample_name = re.search("N_(.+?)_EKDL", file_frag).group(1)
            
            df_dam_base[5] = df_dam_base[1] + df_dam_base[3]
            

            ax.bar(((df_currated[2] - df_currated[1]) / 2) + df_currated[1],
                    height = df_currated[3], 
                    width = df_currated[2] - df_currated[1],
                    color = ["green" if i > 500 else "blue" for i in df_currated[5]],
                    edgecolor = "black")

            ax.plot(df_dam_base[5], 
                     df_dam_base[4], 
                     color = "red")
            
            
            ax.set_ylabel("reads / frag_length")
            ax.set_xlabel("positions")
            
            
            ax.set_title(sample_name)

            fig.suptitle(f"Mean number of reads per {bin_size}")


            path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/base_analysis/"
                    + "analysis_base_frag"
                    + sample_name
                    + "_"
                    + representation
                    +".eps")

            
            plt.savefig(path, format = "eps")
            plt.clf()
            plt.close()
            print(path)


