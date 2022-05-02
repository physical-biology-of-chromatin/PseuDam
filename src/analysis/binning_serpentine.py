from packages.binning import (binning_serpentine_kalisto, bed_fragments_generation, mean_replicates, sleuth_zeros_addition)
from packages.Dam_ID_utilities import (kallisto_abundance_extraction, sleuth_norm_extraction)
import pandas as pd
import numpy as np




def main():
    control_files_list = [
        "/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L2/counts_norm/srf-3i1-NLS-GFP_L2_rep2.csv",
        "/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L2/counts_norm/srf-3i1-NLS-GFP_L2_rep1.csv"
    ]

    test_files_list = [
        "/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L2/counts_norm/srf-3i1-rpb-6_L2_rep1.csv",
        "/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L2/counts_norm/srf-3i1-rpb-6_L2_rep2.csv"
    ]

    chrom_names_celegans = {
        "ENA|BX284601|BX284601.5" : "chromosome_I",
        "ENA|BX284602|BX284602.5" : "chromosome_II",
        "ENA|BX284603|BX284603.4" : "chromosome_III",
        "ENA|BX284604|BX284604.4" : "chromosome_IV",
        "ENA|BX284605|BX284605.5" : "chromosome_V",
        "ENA|BX284606|BX284606.5" : "chromosome_X"
        }


    bed_file = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/Dam_ID/c_elegans/sites/sites.bed"

    epsilon_list = [10, 50, 100, 200]

    teta_list = [100, 500, 1000, 2000]

    df_test_list = list()

    for test_file in test_files_list:
        
        df_test = sleuth_norm_extraction(test_file, 
                                        chrom_names_celegans)

        
        df_test_list.append(df_test)



    df_control_list = list()    

    for control_file in control_files_list:
        
        df_control = sleuth_norm_extraction(control_file, 
                                            chrom_names_celegans)
        
        df_control_list.append(df_control)


    df_control_mean = mean_replicates(df_control_list)
    df_test_mean = mean_replicates(df_test_list)

    df_test_mean = sleuth_zeros_addition(df_test_mean, bed_file)
    print(f"length df = {len(df_test_mean)}")
    df_control_mean = sleuth_zeros_addition(df_control_mean, bed_file)
    print(f"length df = {len(df_test_mean)}\n")
    
    
    for teta in teta_list:
        
        for epsilon in epsilon_list:
            
            if epsilon >= teta:
                continue
            
            df_bins_test, df_bins_control = binning_serpentine_kalisto(df_test_mean,
                                                                    df_control_mean,
                                                                    teta=teta,
                                                                    epsilon=epsilon)
            
            path_bed = f"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/binned_{epsilon}_{teta}.bed"
            
            path_tx = f"/datas/nathan/vscode_nextflow/nextflow-nathan/data/tx2gene/tx2gene_binned_{epsilon}_{teta}.txt"
            
            print(f"epsilon : {epsilon} \n"
                f"teta = {teta}\n")


        
            bed_fragments_generation(df_bins_test, 
                                    path_bed,
                                    path_tx2gene=path_tx)
        

main()