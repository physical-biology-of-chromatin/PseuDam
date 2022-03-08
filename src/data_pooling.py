import pandas as pd
import pybedtools

file1_list = ["/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_2RG1_EKDL220001501-1a_H2LFNDSX3_L4_cov.bed",
              "/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG1D2_EKDL220001507-1a_H2KMHDSX3_L4_cov.bed"]
file2_list = ["/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_2RG1_EKDL220001501-1a_H25VJDSX3_L3_cov.bed",
              "/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_3RG1D2_EKDL220001507-1a_H25VJDSX3_L3_cov.bed"]

new_names_list = ["N_2RG1_EKDL220001501-1a_H2LFNDSX3_fused_cov.bed","N_3RG1D2_EKDL220001507-1a_H25VJDSX3_fused_cov.bed"]

for file1, file2, file_name in zip(file1_list, file2_list, new_names_list):
    df1 = pd.read_csv(file1, header = None, sep = "\t")
    df2 = pd.read_csv(file2, header = None, sep = "\t")
    
    print(df1)
    print(df2)
    df3 = df1[[0, 1, 2]]
    df3.insert(loc = 3, column = "3", value = df1[3] + df2[3])
    df3.insert(loc = 4, column = "4", value = pd.NaT)
    df3.insert(loc = 5, column = "5", value = df1[5])
    df3.insert(loc = 6, column = "6", value = pd.NaT)
    
    bed = pybedtools.BedTool.from_dataframe(df3)
    path = f"/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/" + file_name 
    
    bed.saveas(path)
    
    