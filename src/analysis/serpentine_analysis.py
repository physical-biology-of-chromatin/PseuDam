import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from packages.Dam_ID_utilities import sleuth_norm_extraction



files_list = [
    "file_1",
    "file_n"
]

chrom_names_celegans = {
    "ENA|BX284601|BX284601.5" : "chromosome_I",
    "ENA|BX284602|BX284602.5" : "chromosome_II",
    "ENA|BX284603|BX284603.4" : "chromosome_III",
    "ENA|BX284604|BX284604.4" : "chromosome_IV",
    "ENA|BX284605|BX284605.5" : "chromosome_V",
    "ENA|BX284606|BX284606.5" : "chromosome_X"
    }


for file in files_list:
    
    df_norm = sleuth_norm_extraction(file, chrom_names_celegans)
    
    plt.hist(df_norm["length"],
             bins = 200)
    
    length = np.arange(300, df_norm["length"].max(), 50)
    
    
    mean_reads_range = list()
    median_reads_range = list()
    std_reads_range = list()
    positions = list()
    
    plt.clf()
    
    for i in range(len(length)-1):
        
        df_norm_filtered = df_norm[(df_norm["length"] < length[i])
                                   & df_norm["length"] > length[i+1]]
        
        mean_reads = df_norm_filtered["est_counts"].mean()
        median_reads = df_norm_filtered["est_counts"].median()
        std_reads = df_norm_filtered["est_counts"].std()
        
        mean_reads_range.append(mean_reads)
        median_reads_range.append(median_reads)
        std_reads_range.append(std_reads)
        
        positions.append(length[i] + 50)
        
        plt.boxplot(df_norm_filtered["est_counts"])
        
    plt.xticks(labels=[f"{position - 50} - {position + 50}" for position in positions])