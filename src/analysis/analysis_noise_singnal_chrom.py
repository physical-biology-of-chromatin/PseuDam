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

norm = "no_norm"

choice = "damc"

bin_sizes = [500]

representation_list = ["whole"]

window_list = [6000000]

files_values = list()






for chrom_number in range(33, 49):

    chrom_name = f"ref|NC_0011{str(chrom_number)}|"
    
    fig, ax = plt.subplots(3,
                           3,
                           sharex = True,
                           sharey = True,
                           figsize = (25,25))
    
    col = 0
    rows = 0
    
    for file_signal in files_signal_list:
        
        sample_noise_name = re.search("G(.+?)D", file_signal).group(1)
        
        sample_name = re.search("N_(.+?)_EKDL", file_signal).group(1)
        
        df_noise = pd.read_csv(files_noise_list[(int(re.search("N_3RG(.+?)D", file_signal).group(1)) - 1)],
                            sep = "\t",
                            header = None)

        
        for representation, window_size in zip(representation_list, window_list):

            for bin_size in bin_sizes:

                
                df_dam = pd.read_csv(file_signal, 
                                    header = None, 
                                    sep = '\t')





                # Selection of the rows
                # Chosen chormosome / in the window
                
                if representation == "whole":
                    
                    df_currated = df_dam[df_dam[0] == chrom_name] 
                    df_noise = df_noise[df_noise[0] == chrom_name]
                    
                else:   
                    df_currated = df_dam[df_dam[0] == chrom_name] 
                    df_noise = df_noise[df_noise[0] == chrom_name]

                    
                    df_currated = df_currated[df_currated[1] > 92500 - (window_size/2)]
                    df_currated = df_currated[df_currated[2] < 92500 + (window_size/2)] 

                    df_noise = df_noise[df_noise[1] > 92500 - (window_size/2)]
                    df_noise = df_noise[df_noise[2] < 92500 + (window_size/2)]
                
                
                # Addition of a column contianing the reads normalised by the size
                df_currated.insert(loc = 7,
                                column = 7,
                                value = df_currated[3] / df_currated[5])

                
                df_noise[7] = df_noise[3] / df_noise[5]
                
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

                ##########################################################################
                #                               Binning noise                            #
                ##########################################################################
                
                # Transfert des reads des fragments dans les bins
                for row_noise in df_noise.itertuples():
                    
                    frag_number += 1
                    frag_width = row_noise[3] - row_noise[2]
                    
                    if norm == "no_norm":
                        frag_reads = row_noise[4]
                    
                    elif norm == "norm":
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


                ##########################################################################
                #                               Binning signal                           #
                ##########################################################################

                i = 1
                value_region = 0

                regions_value_list = list()
                regions_pos_list = list()

                gatc = 0
                frag_number = 0


                for row in df_currated.itertuples():
                    
                    frag_number += 1
                    frag_width = row[3] - row[2]
                    
                    if norm == "no_norm":
                        frag_reads = row[4]
                    
                    elif norm == "norm":
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
        
        regions_value_list = [i+1 for i in regions_value_list]
        regions_value_noise_list = [i+1 for i in regions_value_noise_list]
        




        
        
                ##########################################################################
                #                               Plotting                                 #
                ##########################################################################
                
                
        if choice == "log2ratio":
            regions_values_ratio = list(map(truediv, 
                                            regions_value_list, 
                                            regions_value_noise_list))
            
            regions_values_log2ratio = np.log2(regions_values_ratio)
            
            ax[rows,col].bar(regions_pos_list,
                            height = regions_values_log2ratio,
                            width = bin_size,
                            color = ["red" if i < 0 else "blue" for i in regions_values_log2ratio])
        
        elif choice == "ratio":
            regions_values_ratio = list(map(truediv, 
                                            regions_value_list, 
                                            regions_value_noise_list))
            
            ax[rows,col].bar(regions_pos_list,
                            height = regions_values_ratio,
                            width = bin_size)
        
        elif choice == "square":
            regions_values_square = list(map(lambda list1, list2: np.sqrt(list1 / list2) - 1, 
                                         regions_value_list, 
                                         regions_value_noise_list))
            
            ax[rows, col].bar(regions_pos_list,
                              height = regions_values_square,
                              width = bin_size,
                              color = ["red" if i < 0 else "blue" for i in regions_values_square])
        
        elif choice == "damc":
            regions_values_damc = list(
                map(lambda signal, noise : (signal - noise) / signal, 
                                           regions_value_list, 
                                           regions_value_noise_list))
            
            ax[rows, col].bar(regions_pos_list,
                              height = regions_values_damc,
                              width = bin_size,
                              color = ["red" if i < 0 else "blue" for i in regions_values_damc])
            
        ax[rows,col].set_ylabel(choice)
        ax[rows,col].set_xlabel("positions")

        ax[rows,col].set_title(sample_name)

        col += 1
        if col == 3:
            rows += 1
            col = 0
        if col == 1 and rows == 1:
            col += 1
        
        
    path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/figures/"
            + "analysis_"
            + choice
            + "_"
            + norm
            + "_"
            + representation
            + "_depth_bins_"
            + str(bin_size)
            + "_chrom_"
            + str(chrom_number)
            +".eps")

    plt.savefig(path, format = "eps")
    plt.clf()
    plt.close()
    print(path)


