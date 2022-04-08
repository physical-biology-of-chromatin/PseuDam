from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import re

files_noise_list = {
    "3RG1" : "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3/abundance.tsv",
    "3RG2" : "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG2_EKDL220001505-1a_H2C2YDSX3_L3/abundance.tsv",
    "3RG3" : "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4/abundance.tsv"}



files_signal_list = [
    "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG2D1_EKDL220001509-1a_H2C2YDSX3_L3/abundance.tsv",
    "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D1_EKDL220001511-1a_H2C2YDSX3_L3/abundance.tsv",
    "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1D2_/abundance.tsv",
    "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D2_EKDL220001512-1a_H2C2YDSX3_L3/abundance.tsv",
    "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1D3_EKDL220001508-1a_H2KMHDSX3_L4/abundance.tsv",
    "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4/abundance.tsv",
    "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG2D3_EKDL220001510-1a_H2C2YDSX3_L3/abundance.tsv",
    "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG3D3_EKDL220001513-1a_H2C2YDSX3_L3/abundance.tsv"
]


chrom_names = {
    "ref|NC_001133|" : "chrom 1",
    "ref|NC_001134|" : "chrom 2",
    "ref|NC_001135|" : "chrom 3",
    "ref|NC_001136|" : "chrom 4",
    "ref|NC_001137|" : "chrom 5",
    "ref|NC_001138|" : "chrom 6",
    "ref|NC_001139|" : "chrom 7",
    "ref|NC_001140|" : "chrom 8",
    "ref|NC_001141|" : "chrom 9",
    "ref|NC_001142|" : "chrom 10",
    "ref|NC_001143|" : "chrom 11",
    "ref|NC_001144|" : "chrom 12",
    "ref|NC_001145|" : "chrom 13",
    "ref|NC_001146|" : "chrom 14",
    "ref|NC_001147|" : "chrom 15",
    "ref|NC_001148|" : "chrom 16",
    "ref|NC_001224|" : "mitoch",
}


file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/test/yes/sites.bed"

for file_signal in files_signal_list:

    df_signal = pd.read_csv(file_signal,
                            sep = "\t",
                            header = 0)


    df_sites = pd.read_csv(file_sites,
                        sep = "\t",
                        header = None)
    df_signal.insert(0,
                    column = "chrom",
                    value = df_sites[0])
    df_signal.insert(1,
                    column = "start",
                    value = df_sites[1])
    df_signal.insert(2,
                    column = "stop",
                    value = df_sites[2])
    df_signal.iloc[0,3] = "_0"

    df_signal = df_signal[df_signal["chrom"] != "ref|NC_001144|"]

    sample_name = re.search("N_(.+?)D", file_signal).group(1)

    file_noise = files_noise_list[sample_name]


    df_noise = pd.read_csv(file_noise,
                            sep = "\t",
                            header = 0)


    df_sites = pd.read_csv(file_sites,
                        sep = "\t",
                        header = None)
    df_noise.insert(0,
                    column = "chrom",
                    value = df_sites[0])
    df_noise.insert(1,
                    column = "start",
                    value = df_sites[1])
    df_noise.insert(2,
                    column = "stop",
                    value = df_sites[2])
    df_noise.iloc[0,3] = "_0"

    df_signal["ratio"] = np.log2((df_signal["est_counts"] + 1) / (df_noise["est_counts"] + 1))

    df_largest = df_signal.nlargest(10, columns = "ratio")


    fig, ax = plt.subplots(2,5,
                           sharey = True,
                           figsize = (50,30))


    print(df_largest)

    rows = 0
    col = 0

    for row in df_largest.itertuples():
        
        window = 20000
        
        df_signal_plot = df_signal[df_signal["chrom"] == row[1]]
        df_signal_plot = df_signal_plot[df_signal_plot["start"] >= row[2] - (window / 2)]
        df_signal_plot = df_signal_plot[df_signal_plot["stop"] <= row[3] + (window / 2)]


        df_noise_plot = df_noise[df_noise["chrom"] == row[1]]
        df_noise_plot = df_noise_plot[df_noise_plot["start"] >= row[2] - (window / 2)]
        df_noise_plot = df_noise_plot[df_noise_plot["stop"] <= row[3] + (window / 2)]
        
        
        ax[rows,col].bar(df_signal_plot["start"] + (df_signal_plot["length"] / 2),
                         height = np.log2((df_signal_plot["est_counts"] + 1)/( df_noise_plot["est_counts"] + 1)),
                         width = df_signal_plot["length"],
                         color = ["green" if ratio == row[9] else 
                                 "blue" if ratio > 0 else 
                                 "red" 
                                for ratio in df_signal_plot["ratio"]])
        
        ax[rows,col].set_ylabel("log2(signal/noise)")
        ax[rows,col].set_xlabel("positions")
        
        ax[rows,col].set_title(chrom_names[row[1]])
        
        col += 1
        
        if col == 5:
            col = 0
            rows += 1

    sample_name = re.search("N_(.+?)_", file_signal).group(1)
    
    fig.suptitle(f"regions surrounding the 10 largest signal/noise ratio peaks\n"
                 f"{sample_name}")

    path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/figures/highest_peak_ratio/"
            + sample_name
            + "_10_highest_peaks_ratio.eps")
    
    plt.savefig(path, format = "eps")
    
    