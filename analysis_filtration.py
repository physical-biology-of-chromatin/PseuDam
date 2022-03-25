import pandas as pd
import matplotlib.pyplot as plt
import re


file_in_list = ["/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_2RG1_EKDL220001501-1a_H2LFNDSX3_L4.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_2RG1_EKDL220001501-1a_H25VJDSX3_L3.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_2RG2_EKDL220001502-1a_H2C2YDSX3_L3.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG1D2_EKDL220001507-1a_H2KMHDSX3_L4.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG1D2_EKDL220001507-1a_H25VJDSX3_L3.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3.bed",
                "/datas/nathan/vscode_nextflow/nextflow-nathan/results/raw_bed/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3.bed"]

file_out_list = ["/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_2RG1_EKDL220001501-1a_H2LFNDSX3_L4_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_2RG1_EKDL220001501-1a_H25VJDSX3_L3_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_2RG2_EKDL220001502-1a_H2C2YDSX3_L3_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG1D2_EKDL220001507-1a_H2KMHDSX3_L4_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG1D2_EKDL220001507-1a_H25VJDSX3_L3_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3_inter.bed",
                 "/datas/nathan/vscode_nextflow/nextflow-nathan/results/intersect/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3_inter.bed"]

fig, ax = plt.subplots(4,
                       4,
                       sharex=True,
                       sharey=True,
                       figsize = (30,30))

row = 0

col = 0

i = 0

for file_in, file_out in zip(file_in_list, file_out_list):

    df_in = pd.read_csv(file_in, 
                        header = None, 
                        sep = '\t')

    df_out = pd.read_csv(file_out,
                        header = None,
                        sep = "\t")

    sample_name = re.search("N_(.+?)_EKDL", file_in).group(1)
    
    print(df_in)

    df_in.insert(loc = 6, 
                 column = 6, 
                 value = (df_in[2] - df_in[1]))

    ax[row, col].hist(df_in[6],
                      color = "blue",
                      label = "before")

    print(df_out)

    df_out.insert(loc = 6, 
                  column = 6,
                  value = (df_out[2] - df_out[1]))

    
    ax[row, col].hist(df_out[6],
                color = "red",
                label = "after")

    ax[row, col].set_title(f"{sample_name} \n {round((1 - (df_out.shape[0] / df_in.shape[0]))*100)}% of the reads were filtered out ")    

    file_name = re.search("N_(.+?)_EKDL", file_in).group(1)
    
    
    mean = round((1 - (df_out.shape[0] / df_in.shape[0]))*100)
    
    fig.suptitle(f"{file_name}\n"
                 f"{mean}% of the reads were filtered out")
    
    col += 1   
    
    if col == 4:
        col = 0
        row +=1



path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
        +"filter_analysis"
        )

plt.legend()
plt.savefig(path)
plt.clf()
plt.close()