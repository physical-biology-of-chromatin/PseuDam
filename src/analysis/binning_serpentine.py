from packages.binning import (bed_fragments_generation,
                              binning_serpentine_kalisto, mean_replicates,
                              sleuth_zeros_addition)
from packages.Dam_ID_utilities import (sleuth_norm_extraction_cel)




def binning_celegans(epsilon, teta, m_length):
    control_files_list = [
    "/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L4/counts_norm/srf3_L4_control_1.csv",
    "/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L4/counts_norm/srf3_L4_control_2.csv"
    ]

    test_files_list = [
    "/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L4/counts_norm/srf3_L4_test_1.csv",
    "/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L4/counts_norm/srf3_L4_test_2.csv"
    ]

    chrom_names = {
        "ENA|BX284601|BX284601.5" : "chrI",
        "ENA|BX284602|BX284602.5" : "chrII",
        "ENA|BX284603|BX284603.4" : "chrIII",
        "ENA|BX284604|BX284604.4" : "chrIV",
        "ENA|BX284605|BX284605.5" : "chrV",
        "ENA|BX284606|BX284606.5" : "chrX"
        }
    

    bed_file = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/Dam_ID/c_elegans/sites/sites.bed"

    df_test_list = list()

    for test_file in test_files_list:
        
        df_test = sleuth_norm_extraction_cel(test_file, 
                                             chrom_names)

        
        df_test_list.append(df_test)



    df_control_list = list()    

    for control_file in control_files_list:
        
        df_control = sleuth_norm_extraction_cel(control_file, 
                                            chrom_names)
        
        df_control_list.append(df_control)


    df_control_mean = mean_replicates(df_control_list)
    df_test_mean = mean_replicates(df_test_list)

    df_test_mean = sleuth_zeros_addition(df_test_mean, bed_file, chrom_names)
    print(f"length df = {len(df_test_mean)}")
    df_control_mean = sleuth_zeros_addition(df_control_mean, bed_file, chrom_names)
    print(f"length df = {len(df_test_mean)}\n")
    

    df_bins_test= binning_serpentine_kalisto(df_test_mean,
                                                df_control_mean,
                                                teta=teta,
                                                epsilon=epsilon,
                                                m_length=m_length)
    
    path_bed = f"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/srf3_L4/sites_{epsilon}_{teta}_{m_length}.bed"
    
    
    print(f"epsilon : {epsilon} \n"
        f"teta = {teta}\n")



    bed_fragments_generation(df_bins_test, 
                            path_bed,
                            chrom_names)
        
        

        
def binning_Damc():
    control_files_list = [
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG2_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG1_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/2RG3_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/2RG2_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/2RG1_norm.csv"
    ]

    test_files_list = [
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG3D3_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG3D2_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG3D1_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG2D3_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG2D1_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG1D3_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG1D2_norm.csv",
    "/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG1D1_norm.csv"
    ]

    chrom_names = {
        "chrom_I" : "chrI",
        "chrom_II" : "chrII",
        "chrom_III" : "chrIII",
        "chrom_IV" : "chrIV",
        "chrom_V" : "chrV",
        "chrom_VI" : "chrVI",
        "chrom_VII" : "chrVII",
        "chrom_VIII" : "chrVIII",
        "chrom_IX" : "chrIX",
        "chrom_X" : "chrX",
        "chrom_XI" : "chrXI",
        "chrom_XII" : "chrXII",
        "chrom_XIII" : "chrXIII",
        "chrom_XIV" : "chrXIV",
        "chrom_XV" : "chrXV",
        "chrom_XVI" : "chrXVI",
        "chrom_circular" : "chrM",
        }
    
    print(f"epsilon : {epsilon} \n"
    f"teta : {teta}\n")

    bed_file = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/Dam_ID/DamC_split/sites_test/sites.bed"

    epsilon_list = [10]

    teta_list = [100]

    df_test_list = list()

    for test_file in test_files_list:
        
        df_test = sleuth_norm_extraction_cel(test_file, 
                                             chrom_names)

        
        df_test_list.append(df_test)



    df_control_list = list()    

    for control_file in control_files_list:
        
        df_control = sleuth_norm_extraction_cel(control_file, 
                                            chrom_names)
        
        df_control_list.append(df_control)


    df_control_mean = mean_replicates(df_control_list)
    df_test_mean = mean_replicates(df_test_list)

    df_test_mean = sleuth_zeros_addition(df_test_mean, bed_file, chrom_names)

    df_control_mean = sleuth_zeros_addition(df_control_mean, bed_file, chrom_names)

    
    
    for teta in teta_list:
        
        for epsilon in epsilon_list:
            
            if epsilon >= teta:
                continue
            
            df_bins_test= binning_serpentine_kalisto(df_test_mean,
                                                                    df_control_mean,
                                                                    teta=teta,
                                                                    epsilon=epsilon)
            
            path_bed = f"/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/dam_total/sites_10_100.bed"
            
            



        
            bed_fragments_generation(df_bins_test, 
                                    path_bed,
                                    chrom_names)
    
    print("Done")


binning_celegans(epsilon=10, teta=500, m_length=3000)

binning_celegans(epsilon=50, teta=100, m_length=3000)
