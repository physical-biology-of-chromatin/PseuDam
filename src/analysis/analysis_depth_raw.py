import re

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

# Importing the data as a pandas dataframe

window_size_list = [10000, 20000, 50000, 100000, 200000]

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


fig, ax = plt.subplots(nrows = 3,
                       ncols = 3,
                       sharex = True,
                       sharey = True)



for file in files_list:

    for window_size in window_size_list:
        
        df_dam = pd.read_csv(file, 
                                header = None, 
                                sep = '\t')
                
        chrom_name = "ref|NC_001135|"
        
        # Selection of the rows
        # Chosen chormosome 
        df_currated = df_dam[df_dam[0] == chrom_name]
        df_currated = df_currated[df_currated[1] > 90000 - (window_size/2)]
        df_currated = df_currated[df_currated[2] < 90000 + (window_size/2)] 
        
        
        file_name = re.search("N_(.+?)_EKDL", file).group(1)
        

        ax.bar(((df_currated[2] - df_currated[1]) / 2) + df_currated[1],
                height = df_currated[3], 
                width = df_currated[2] - df_currated[1])
        ax.set_ylabel("read/fragment")
        ax.set_xlabel("postitions")
        
        
        TetO = patches.Rectangle(xy = (92500,-100),
                                width = 1000,
                                height = 200,
                                color = "black", 
                                clip_on = False )

        ax.add_patch(TetO)
        
        
        ax.set_ylim(0, 8000)
        
        ax.set_title(file_name)
        
        
    path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
            + "raw/"
            + "analysis_raw_"
            + file_name
            + "_depth_frag_window_"
            + str(window_size))
    plt.savefig(path, format = "eps")
    plt.clf()
    plt.close()

