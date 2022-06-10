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

            fig, ax = plt.subplots( 3, 
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

                if representation == "whole":

                    df_currated = df_dam[df_dam[0] == chrom_name] 

                else:   
                    df_currated = df_dam[df_dam[0] == chrom_name] 



                    df_currated = df_currated[df_currated[1] > 92500 - (window_size/2)]
                    df_currated = df_currated[df_currated[2] < 92500 + (window_size/2)] 

                # Finds the value of G in the file name to use the right noise file
                df_noise = pd.read_csv(files_noise_list[int(re.search("N_3RG(.+?)D", file).group(1)) - 1],
                                        sep = "\t",
                                        header = None)


                df_currated.insert(loc = 7,
                                    column = 7,
                                    value = df_currated[3] / df_noise[3])

                # Addition of a column contianing the reads normalised by the size
                df_currated.insert( loc = 8,
                                    column = 8,
                                    value = df_currated[7] / df_currated[5])

                

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
                    frag_reads = row[9]
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


                sample_name = re.search("N_(.+?)_EKDL", file).group(1)


                ax[rows,col].bar(regions_pos_list,
                                 height = np.log2(regions_value_list), 
                                 width = bin_size,
                                 color = ["red" if i < 0 else "blue" for i in regions])

                ax[rows,col].set_ylabel("reads / frag_length")
                ax[rows,col].set_xlabel("positions")


                col += 1

                if col == 3:
                    rows += 1
                    col = 0

                if col == 1 and rows == 1:
                    col += 1

            path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
                    + "figures/test/"
                    + "analysis_norm_chrom"
                    + str(chrom)
                    + "_"
                    + representation
                    + "_depth_bins_"
                    + str(bin_size))


            plt.savefig(path)
            plt.clf()
            plt.close()

            print(path)