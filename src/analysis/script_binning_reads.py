import os

from packages.binning import (abundance_file_conversion,
                              binning_application_kallisto,
                              binning_kallisto_reads)
from packages.Dam_ID_utilities import kallisto_abundance_reformating

epsilon_list_reads = [
    500,
    1000,
    2000,
    5000,
    10000
]


file_name = "abundance.tsv"


for condition_gal in range(1, 4):
    
    condition_list = [f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}D1/abundance.tsv",
                      f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}D2/abundance.tsv",
                      f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}D3/abundance.tsv",
                      f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}/abundance.tsv",
                      f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/2RG{condition_gal}/abundance.tsv"]


    df_control = kallisto_abundance_reformating(condition_list[4])

    for epsilon in epsilon_list_reads:
        
        df_binning_control = binning_kallisto_reads(df_control, epsilon = epsilon)
        
        for i, condition in enumerate(condition_list):       

            
            if not os.path.exists(condition):
                continue

            if i != 4:    
                
                df_condition = kallisto_abundance_reformating(condition)
                    
                df_binning_reads = binning_application_kallisto(df_condition, df_binning_control)
                

            
            if i < 3:
                path = f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}D{i+1}/bin_reads_{epsilon}.tsv"
                
                abundance_file_conversion(df_binning_reads, path)
            
            
            if i == 3:
                path = f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}/bin_reads_{epsilon}.tsv"

                abundance_file_conversion(df_binning_reads, path)

            
            if i == 4:
                path = f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/2RG{condition_gal}/bin_reads_{epsilon}.tsv"
                path_tx2gene = f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/tx2gene_bin_reads_{epsilon}.tsv"
                abundance_file_conversion(df_binning_control, path, path_tx2gene = path_tx2gene)
