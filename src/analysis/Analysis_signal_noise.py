from operator import truediv
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


files_noise_list = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/coverage_per_frag/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3_cov.bed",
                    "/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/coverage_per_frag/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3_cov.bed",
                    "/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/coverage_per_frag/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4_cov.bed")



files_signal_list = [
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG1D2_EKDL220001507-1a_H25VJDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3_cov.bed"]


bin_sizes = [500]


representation_list = ["whole"]

window_list = [600000]

chrom_list = np.arange(33, 48, 1)


for chrom in chrom_list:
    
    for bin_size in bin_sizes:

        for representation, window_size in zip(representation_list, window_list):

            fig, ax = plt.subplots(3, 
                                    3,
                                    sharex=True,
                                    sharey=True,
                                    figsize = (25,25))

            col = 0
            rows = 0

            for file in files_signal_list:

                df_dam = pd.read_csv(file, 
                                    header = None, 
                                    sep = '\t')



                chrom_name = f"ref|NC_0011{chrom}|"

                # Selection of the rows
                # Chosen chormosome / in the window

                # Finds the value of G in the file name to use the right noise file
                df_noise = pd.read_csv(files_noise_list[int(re.search("N_3RG(.+?)D", file).group(1)) - 1],
                                        sep = "\t",
                                        header = None)

                if representation == "whole":

                    df_currated = df_dam[df_dam[0] == chrom_name] 
                    df_noise = df_noise[df_noise[0] == chrom_name]
                    
                else:   
                    df_currated = df_dam[df_dam[0] == chrom_name] 
                    df_noise = df_noise[df_noise[0] == chrom_name]

                    df_currated = df_currated[df_currated[1] > 92500 - (window_size/2)]
                    df_currated = df_currated[df_currated[2] < 92500 + (window_size/2)] 
                    df_noise = df_noise[df_noise[1] > 92500 - (window_size/2)]
                    df_noise = df_noise[df_noise[1] < 92500 + (window_size/2)]


                df_noise.insert(loc = 7,
                                column = 7,
                                value = df_noise[3] / df_noise[5])

                df_currated.insert(loc = 7,
                                    column = 7,
                                    value = df_currated[3] / df_currated[5])


                # Addition of a column contianing the reads normalised by the size


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

                regions_value_noise_list = list()
                regions_pos_list = list()

                gatc = 0
                frag_number = 0     
                            
                for row_noise in df_noise.itertuples():
        
                    frag_number += 1
                    frag_width = row_noise[3] - row_noise[2]
                    frag_reads = row_noise[8]
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


                # Transfert des reads des fragments dans les bins

                for row_signal in df_currated.itertuples():

                    frag_number += 1
                    frag_width = row_signal[3] - row_signal[2]
                    frag_reads = row_signal[8]
                    region = regions[i]

                    # skips the first few bins without fragments
                    if i == 1:
                        i = 1
                        if regions[i] < row_signal[2]:

                            while not (regions[i-1] < row_signal[2] < regions[i]):

                                middle_region = (((regions[i] - regions[i-1]) / 2) + regions[i - 1])
                                regions_pos_list.append(middle_region)

                                regions_value_list.append(0)
                                i += 1


                    # fragment fully inside the bin
                    if (row_signal[3] < regions[i]
                        and row_signal[2] > regions[i-1]):

                        value_region += frag_reads

                    # fragment with one end outside
                    elif row_signal[3] > regions[i]:

                        value_region += frag_reads

                        #fragment over mulitple bins
                        if row_signal[3] > regions[i+1]:
                        
                            while row_signal[3] > regions[i+1]:

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
                
                
                
                regions_ratio_list = [i / j for i, j in zip(regions_value_list, regions_value_noise_list)]
                
                regions_log2ratio_list = np.log2(regions_ratio_list)
                

                
                sample_name = re.search("N_(.+?)_EKDL", file).group(1)


                ax[rows,col].bar(regions_pos_list,
                                height = regions_log2ratio_list, 
                                width = bin_size)

                ax[rows,col].set_ylabel("log2(Dox+/Dox-)")
                ax[rows,col].set_xlabel("positions")


                col += 1

                if col == 3:
                    rows += 1
                    col = 0

                if col == 1 and rows == 1:
                    col += 1

            path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
                    + "figures/test/"
                    + "analysis_norm_log2ratio_end_chrom"
                    + str(chrom)
                    + "_"
                    + representation
                    + "_depth_bins_"
                    + str(bin_size))


            plt.savefig(path)
            plt.clf()
            plt.close()

            print(path)