import os

import pandas as pd

from packages.binning import (abundance_file_conversion,
                              binning_kallisto_gatc)
from packages.Dam_ID_utilities import kallisto_abundance_reformating

epsilon_list_gatc = [
    2,
    5,
    10
]

file_name = "abundance.tsv"

for condition_gal in range(1, 4):

    condition_list = [f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}D1/abundance.tsv",
                      f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}D2/abundance.tsv",
                      f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}D3/abundance.tsv",
                      f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}/abundance.tsv",
                      f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/2RG{condition_gal}/abundance.tsv"]
    
    for i, condition in enumerate(condition_list):
        
        if os.path.exists(condition) != True:
            continue
        
        df_condition = kallisto_abundance_reformating(condition)

        
        for epsilon in epsilon_list_gatc:
        
            df_binning_gatc = binning_kallisto_gatc(df_condition, epsilon = epsilon)
            
            
            if i < 3:
                path = f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}D{i+1}/bin_gatc_{epsilon}.tsv"
                abundance_file_conversion(df_binning_gatc, path)
                
                
            if i == 3:
                path = f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/3RG{condition_gal}/bin_gatc_{epsilon}.tsv"
                abundance_file_conversion(df_binning_gatc, path)
                
                
            if i == 4:
                path = f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/2RG{condition_gal}/bin_gatc_{epsilon}.tsv"
                path_tx2gene = f"/datas/nathan/Dam_ID_analysis/data/counts/G{condition_gal}/tx2gene_bin_gatc_{epsilon}.tsv"
                abundance_file_conversion(df_binning_gatc, path, path_tx2gene = path_tx2gene)



            