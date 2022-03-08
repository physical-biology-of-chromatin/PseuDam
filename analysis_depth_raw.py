import re

import matplotlib.pyplot as plt
import pandas as pd

# Importing the data as a pandas dataframe
files_list = [
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_2RG1_EKDL220001501-1a_H2LFNDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_2RG2_EKDL220001502-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG1D2_EKDL220001507-1a_H25VJDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3_cov.bed"]

file_names_list = ["2RG1",
                   "2RG2",
                   "2RG3",
                   "3RG1",
                   "3RG1D1",
                   "3RG1D2",
                   "3RG1D3",
                   "3RG2",
                   "3RG2D1",
                   "3RG2D3",
                   "3RG3D1",
                   "3RG3D2",
                   "3RG3D3"]
for file, file_name in zip(files_list, file_names_list):
    

    df_dam = pd.read_csv(file, 
                        header = None, 
                        sep = '\t')
    
    chrom_name = "ref|NC_001135|"
    
    # Selection of the rows
    # Chosen chormosome 
    df_currated = df_dam[df_dam[0] == chrom_name]
    """
    df_currated = df_currated[df_currated[1] > 90000 - (window_size/2)]
    df_currated = df_currated[df_currated[2] < 90000 + (window_size/2)] 
    """      
    
    #file_name = re.search("N_(.+?)_EKDL", file).group(1)
    

    plt.bar(((df_currated[2] - df_currated[1]) / 2) + df_currated[1],
             height = df_currated[3], 
             width = df_currated[2] - df_currated[1])
    plt.ylabel("read/fragment")
    plt.xlabel("postitions")
    
    plt.ylim(0, 8000)
    
    plt.suptitle(file_name)
    
    path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
            + "raw/reduced_scale"
            + "/"
            + "analysis_raw_"
            + file_name
            + "_depth_frag")
    plt.savefig(path)
    plt.clf()
    plt.close()
