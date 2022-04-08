import re

import scipy.stats

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.signal import find_peaks

files_frag_list = [
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_2RG1_EKDL220001501-1a_H2LFNDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_2RG2_EKDL220001502-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG1D2_EKDL220001507-1a_H25VJDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_frag/old/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3_cov.bed"]


files_base_list = [
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_2RG1_EKDL220001501-1a_H2LFNDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_2RG2_EKDL220001502-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_3RG1D2_EKDL220001507-1a_H25VJDSX3_fused_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3_cov.bed",
"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage_base/old/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3_cov.bed"]

bin_sizes = [500]

center = 92500


representation_list = ["whole"]

window_list = [10000]

distances_peaks_samples = list()
distances_sides_samples = list()

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

                
                df_currated = df_currated[df_currated[1] > center - (window_size/2)]
                df_currated = df_currated[df_currated[2] < center + (window_size/2)]

                df_dam_base = df_dam_base[df_dam_base[1] > center - (window_size/2)]
                df_dam_base = df_dam_base[df_dam_base[2] < center + (window_size/2)] 
            
            distance_max_list = list()
            
            peaks_list = list()
            
            distances_peaks_list = list()
            distances_sides_list = list()
            
            for row in df_currated.itertuples():
                
                df_temp = df_dam_base[df_dam_base[1] == row[2]]
                df_temp = df_temp[df_temp[2] == row[3]]
                
                if row[3] != 0 and len(df_temp) / 3 >= 1:
                    
                    peaks, _ = find_peaks(df_temp[4], 
                                          distance = len(df_temp) / 3)
                    
                    heads = df_temp.head()
                    
                    if len(peaks) > 1:
                        
                        distances_peaks = peaks[1] - peaks[0]
                        
                        distances_peaks_list.append(distances_peaks)

                        distance_to_sides = [peaks[0], 
                                             row[3] - (row[2] + peaks[1])] 
                        
                        distances_sides_list.append(np.mean(distance_to_sides))
                    
                    elif len(peaks) > 0:
                        
                        distance_to_sides = peaks[0]
                        
                        if distance_to_sides < row[3] - (distance_to_sides + row[2]):
                            distances_sides_list.append(distance_to_sides)
                            
                        else:
                            distances_sides_list.append(row[3] - (distance_to_sides + row[2]))
                     
                    positions = [heads.index[0] + peak for peak in peaks]
                    
                    [peaks_list.append([row[2] + peaks[i], 
                                        df_temp[4][position]]) 
                     for position, i in zip(positions, range(0,len(peaks)))]


            
            """
                linemax = df_temp.loc[df_temp[4].idxmax()]
                
                distance_max = linemax[3]
                

                
                if linemax[4] == 0:
                    pass
                
                elif distance_max > row[3] - (distance_max + row[2]): 
                    
                    distance_max_list.append(row[3] - (distance_max + row[2]))
                
                else:
                     
                    distance_max_list.append(distance_max)
            """



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
            
            
            ax.plot([peak[0] for peak in peaks_list],
                    [peak[1] for peak in peaks_list],
                    "x")
            
            
            ax.set_ylabel("reads / frag_length")
            ax.set_xlabel("positions")
            
            
            ax.set_title(sample_name)

            fig.suptitle(f"Mean number of reads per {bin_size}\n"
                         f"Mean distance peak =  {np.mean(distances_peaks_list)}")


            path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/base_analysis/"
                    + "analysis_base_frag"
                    + sample_name
                    + "_"
                    + representation
                    +".eps")
            
            ax.set_xticks(np.arange(center - (window_size/2), center + (window_size/2), 100))
            
            plt.show()
            plt.clf()
            plt.close()
            
            distances_peaks_samples.append([sample_name, np.mean(distances_peaks_list)])
            distances_sides_samples.append([sample_name, np.mean(distances_sides_list)])
            
            
"""            
# T tests for the distances
print("list peaks : ")
print(distances_peaks_samples)

print("list sides : ")
print(distances_sides_samples)

control_peaks = list()
test_peaks = list()

control_sides = list()
test_sides = list()

for sample_peaks, sample_sides in zip(distances_peaks_samples, distances_sides_samples):
    if "D" in sample_peaks[0]:
        test_peaks.append(sample_peaks[1])
        test_sides.append(sample_sides[1])
    else:
        control_peaks.append(sample_peaks[1])
        control_sides.append(sample_sides[1])

print("t-test peaks")
print(scipy.stats.ttest_ind(test_peaks, control_peaks))

print("t-test sides")
print(scipy.stats.ttest_ind(test_sides, control_sides))          
"""

