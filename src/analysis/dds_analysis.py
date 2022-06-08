import re
from os.path import exists

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.fftpack import diff
from scipy.stats import gaussian_kde
from sqlalchemy import null
from yaml import BlockSequenceStartToken

from packages.binning import (binning_application_deseq, binning_deseq_gatc,
                              binning_deseq_reads)
from packages.Dam_ID_utilities import (dam_id_window_filtering, kallisto_abundance_extraction,
                                       kallisto_abundance_reformating,
                                       kallisto_outup_reformating,
                                       position_extraction,
                                       sleuth_output_reformating)

#####################################################################
#                                                                   #
#                             Comparisons                           #
#                                                                   #
#####################################################################


def rlog_comparison(gal_condition):
    
    chrom_names = {
        "ref|NC_001133|" : "chrom_1",
        "ref|NC_001134|" : "chrom_2",
        "ref|NC_001135|" : "chrom_3",
        "ref|NC_001136|" : "chrom_4",
        "ref|NC_001137|" : "chrom_5",
        "ref|NC_001138|" : "chrom_6",
        "ref|NC_001139|" : "chrom_7",
        "ref|NC_001140|" : "chrom_8",
        "ref|NC_001141|" : "chrom_9",
        "ref|NC_001142|" : "chrom_10",
        "ref|NC_001143|" : "chrom_11",
        "ref|NC_001144|" : "chrom_12",
        "ref|NC_001145|" : "chrom_13",
        "ref|NC_001146|" : "chrom_14",
        "ref|NC_001147|" : "chrom_15",
        "ref|NC_001148|" : "chrom_16",
        "ref|NC_001224|" : "mitoch",
    }

    chrom = 3

    window = 0

    center = 0

    i = 0

    conditions_list = [f"3RG{gal_condition}D1",
                       f"3RG{gal_condition}D2", 
                       f"3RG{gal_condition}D3"]

    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures(3,1)
    
    file_rld = "/datas/nathan/Dam_ID_analysis/results/rld_total.csv"


    df_rld = pd.read_csv(file_rld, header = 0, sep = ",")


    df_rld = position_extraction(df_rld)

    df_rld = dam_id_window_filtering(df_rld, chrom, center = center, window = window)


    for condition in conditions_list:


        try:
            df_rld[condition]
        except KeyError:
            print(f"{condition} has not been found in the DF")
            continue
                

        subfigs[i].suptitle(condition)
        
        axs = subfigs[i].subplots(1, 2)


        axs[0].bar(df_rld["start"] + (df_rld["length"]/2),
                    height = df_rld[condition] - df_rld[f"2RG{gal_condition}"],
                    width = df_rld["stop"] - df_rld["start"],
                    color = ["red" if log_ratio < 0 
                            else "blue" 
                            for log_ratio in df_rld[condition] - df_rld[f"2RG{gal_condition}"]])

        axs[0].set_title(f"all\nmean = {round(df_rld[condition].mean(), 2)}")
        axs[0].set_ylabel("rlog(ratio)")
        axs[0].set_xlabel("postions")

        

        


        df_rld = pd.read_csv(file_rld, header = 0, sep = ",")


        df_rld = position_extraction(df_rld)

        df_rld = dam_id_window_filtering(df_rld, chrom, center = center, window = window)



        axs[1].bar(df_rld["start"] + (df_rld["length"]/2),
                    height = df_rld[condition] - df_rld[f"2RG{gal_condition}"],
                    width = df_rld["stop"] - df_rld["start"],
                    color = ["red" if log_ratio < 0 else "blue" for log_ratio in df_rld[condition] - df_rld[f"2RG{gal_condition}"]])

        axs[1].set_title(f"Gal{gal_condition} only\nmean = {round(df_rld[condition].mean(), 2)}")
        axs[1].set_ylabel("rlog(ratio)")
        axs[1].set_xlabel("postions")

        i += 1


    plt.show()



def dds_comparison(gal_condition):


    chrom_names = {
        "ref|NC_001133|" : "chrom_1",
        "ref|NC_001134|" : "chrom_2",
        "ref|NC_001135|" : "chrom_3",
        "ref|NC_001136|" : "chrom_4",
        "ref|NC_001137|" : "chrom_5",
        "ref|NC_001138|" : "chrom_6",
        "ref|NC_001139|" : "chrom_7",
        "ref|NC_001140|" : "chrom_8",
        "ref|NC_001141|" : "chrom_9",
        "ref|NC_001142|" : "chrom_10",
        "ref|NC_001143|" : "chrom_11",
        "ref|NC_001144|" : "chrom_12",
        "ref|NC_001145|" : "chrom_13",
        "ref|NC_001146|" : "chrom_14",
        "ref|NC_001147|" : "chrom_15",
        "ref|NC_001148|" : "chrom_16",
        "ref|NC_001224|" : "mitoch",
    }

    chrom = "mitoch"

    window = 0

    center = 0

    i = 0

    conditions_list = [f"3RG{gal_condition}D1",
                       f"3RG{gal_condition}D2", 
                       f"3RG{gal_condition}D3",
                       f"2RG{gal_condition}"]

    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures(4,1)
    
    file_dds = "/datas/nathan/Dam_ID_analysis/results/dds_total.csv"
    df_dds = pd.read_csv(file_dds, header = 0, sep = ",")   
    df_dds = position_extraction(df_dds)    
    df_dds = dam_id_window_filtering(df_dds, chrom, center = center, window = window)


    file_dds_gal = f"/datas/nathan/Dam_ID_analysis/results/dds_G{gal_condition}.csv"
    df_dds_gal = pd.read_csv(file_dds_gal, header = 0, sep = ",")
    df_dds_gal = position_extraction(df_dds_gal)
    df_dds_gal = dam_id_window_filtering(df_dds_gal, chrom, center = center, window = window)

    max_value = 0
    min_value = pow(10, 10)
    
    for condition in conditions_list:
        
        if df_dds[condition].max() > max_value:
            max_value = df_dds[condition].max()
        
        if df_dds_gal[condition].max() > max_value:
            max_value = df_dds_gal[condition].max()


    for condition in conditions_list:

        try:
            df_dds[condition]
        except KeyError:
            print(f"{condition} has not been found in the DF")
            continue
                

        
        
        axs = subfigs[i].subplots(1, 2)


        axs[0].bar(df_dds["start"] + (df_dds["length"]/2),
                    height = df_dds[condition],
                    width = df_dds["stop"] - df_dds["start"])

        axs[0].set_title(f"all\nmean = {round(df_dds[condition].mean(), 6)}")
        axs[0].set_ylabel("count")
        axs[0].set_xlabel("postions")
        axs[0].set_ylim(0,max_value*1.2)




        axs[1].bar(df_dds_gal["start"] + (df_dds_gal["length"]/2),
                    height = df_dds_gal[condition],
                    width = df_dds_gal["stop"] - df_dds_gal["start"])

        axs[1].set_title(f"Gal{gal_condition} only\nmean = {round(df_dds_gal[condition].mean(), 6)}")
        axs[1].set_ylabel("count")
        axs[1].set_xlabel("postions")
        axs[1].set_ylim(0,max_value*1.2)

        subfigs[i].suptitle(f"{condition}\nmeans ratio = {df_dds[condition].mean() / df_dds_gal[condition].mean()}")

        i += 1

    plt.show()



def raw_comparison(gal_condition):
    """Plots a comparison of the raw read counts after sizefactor for triple
    for each test sample in the given gal condition

    Args:
        gal_condition (_type_): _description_
    """
    
    chrom = "chrom_3"
    center = 0
    window = 0

    i = 0

    conditions_list = [f"3RG{gal_condition}D1",
                    f"3RG{gal_condition}D2",
                    f"3RG{gal_condition}D3",
                    f"2RG{gal_condition}"]
    
    
    
    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures(len(conditions_list),1)





    for condition in conditions_list:

        file_dds = "/datas/nathan/Dam_ID_analysis/results/dds_total.csv"



        df_dds = pd.read_csv(file_dds, header = 0, sep = ",")


        df_dds = position_extraction(df_dds)


        df_dds = dam_id_window_filtering(df_dds, chrom, center = center, window = window)


        try:
            df_dds[condition]
        except KeyError:
            print(f"{condition} has not been found in the DF")
            continue
        
        
        subfigs[i].suptitle(f"{condition}")
        
        axs = subfigs[i].subplots(1, 2,
                                sharex = True,
                                sharey = True)
        
        axs[0].bar(df_dds["start"] + (df_dds["length"]/2),
                    height = df_dds[condition] ,
                    width = df_dds["stop"] - df_dds["start"],
                    color = ["red" if log_ratio < 0 else "blue" for log_ratio in df_dds[condition] - df_dds[f"2RG{gal_condition}"]])

        axs[0].set_title(f"all\nmean = {round(df_dds[condition].mean(), 2)}")
        axs[0].set_ylabel("rlog(ratio)")
        axs[0].set_xlabel("postions")

        
        
        file_dds = f"/datas/nathan/Dam_ID_analysis/results/dds_G{gal_condition}.csv"


        df_dds = pd.read_csv(file_dds, header = 0, sep = ",")


        df_dds = position_extraction(df_dds)


        df_dds = dam_id_window_filtering(df_dds, chrom, center = center, window = window)


        axs[1].bar(df_dds["start"] + (df_dds["length"]/2),
                    height = df_dds[condition],
                    width = df_dds["stop"] - df_dds["start"],
                    color = ["red" 
                            if log_ratio < 0 
                            else "blue" 
                            for log_ratio in df_dds[condition] - df_dds[f"2RG{gal_condition}"]])

        axs[1].set_title(f"Gal{gal_condition} only\nmean = {round(df_dds[condition].mean(), 2)}")
        axs[1].set_ylabel("rlog(ratio)")
        axs[1].set_xlabel("postions")

        i += 1


        print((df_dds.nlargest(5, condition)))
    
    plt.show()



def size_factor_comparison(gal_condition):
    
    chrom_names = {
    "ref|NC_001133|" : "chrom_1",
    "ref|NC_001134|" : "chrom_2",
    "ref|NC_001135|" : "chrom_3",
    "ref|NC_001136|" : "chrom_4",
    "ref|NC_001137|" : "chrom_5",
    "ref|NC_001138|" : "chrom_6",
    "ref|NC_001139|" : "chrom_7",
    "ref|NC_001140|" : "chrom_8",
    "ref|NC_001141|" : "chrom_9",
    "ref|NC_001142|" : "chrom_10",
    "ref|NC_001143|" : "chrom_11",
    "ref|NC_001144|" : "chrom_12",
    "ref|NC_001145|" : "chrom_13",
    "ref|NC_001146|" : "chrom_14",
    "ref|NC_001147|" : "chrom_15",
    "ref|NC_001148|" : "chrom_16",
    "ref|NC_001224|" : "mitoch",
    }



    chrom = "chrom_3"
    center = 0
    window = 0

    i = 0

    conditions_list = [f"3RG{gal_condition}D1",
                    f"3RG{gal_condition}D2",
                    f"3RG{gal_condition}D3",
                    f"2RG{gal_condition}"]
    
    
    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures(len(conditions_list),1)

    file_dds = "/datas/nathan/Dam_ID_analysis/results/dds_total.csv"



    df_dds = pd.read_csv(file_dds, header = 0, sep = ",")


    df_dds = position_extraction(df_dds)

    print(df_dds.nlargest(1,"3RG3D2"))

    df_dds = dam_id_window_filtering(df_dds, chrom, center = center, window = window)

    max_value = 0

    for count, condition in enumerate(conditions_list):
        
        path_file_condition = f"/datas/nathan/Dam_ID_analysis/data/yes/{condition}/abundance.tsv"
        
        path_file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/test/yes/sites.bed"
        
        df_raw = kallisto_outup_reformating(path_file_condition, path_file_sites, chrom_names = chrom_names)
        
        if condition == "3RG3D2":
            print(df_raw["est_counts"].mean())
        
        df_raw = dam_id_window_filtering(df_raw, chrom = chrom, window = window, center = center)

        if df_dds[condition].max() > max_value:
            max_value = df_dds[condition].max()
               
        if df_raw["est_counts"].max() > max_value:
            max_value = df_raw["est_counts"].max()




    for count, condition in enumerate(conditions_list):
          
        try:
            df_dds[condition]
        except KeyError:
            print(f"{condition} has not been found in the DF")
            continue

        
        axs = subfigs[count].subplots(1, 2)
        
        axs[0].bar(df_dds["start"] + (df_dds["length"]/2),
                    height = df_dds[condition] ,
                    width = df_dds["stop"] - df_dds["start"],
                    color = ["red" if log_ratio < 0 else "blue" for log_ratio in df_dds[condition] - df_dds[f"2RG{gal_condition}"]])

        axs[0].set_title(f"size factor\nmean = {round(df_dds[condition].mean(), 2)}")
        axs[0].set_ylabel("rlog(ratio)")
        axs[0].set_xlabel("postions")

        axs[0].set_ylim(-20, max_value*1.3)

        path_file_condition = f"/datas/nathan/Dam_ID_analysis/data/yes/{condition}/abundance.tsv"
        
        path_file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/test/yes/sites.bed"
        
        df_raw = kallisto_outup_reformating(path_file_condition, path_file_sites, chrom_names = chrom_names)
        
        df_raw = dam_id_window_filtering(df_raw, chrom = chrom, window = window, center = center)
        
        axs[1].bar(df_raw["start"] + (df_raw["length"]/2),
                height = df_raw["est_counts"],
                width = df_raw["length"])
        
        axs[1].set_title(f"raw\nmean = {round(df_raw['est_counts'].mean(), 2)}")
        axs[1].set_ylabel("read count")
        axs[1].set_xlabel("postions")
        
        axs[1].set_ylim(-20, max_value*1.3)
        
        subfigs[count].suptitle(f"{condition}\n{round(df_dds[condition].mean()/df_raw['est_counts'].mean(), 5)}")
          
    plt.show()


#####################################################################
#                                                                   #
#                           Scatterplots                            #
#                                                                   #
#####################################################################


def scatter_rld(gal_condition):

    chrom = "chrom_3"
    center = 0
    window = 0

    i = 0

    conditions_list = [f"3RG{gal_condition}D1",
                       f"3RG{gal_condition}D2",
                       f"3RG{gal_condition}D3"]

    
    control = f"2RG{gal_condition}"


    file_rld_total = "/datas/nathan/Dam_ID_analysis/results/rld_total.csv"

    df_rld_total = pd.read_csv(file_rld_total, header = 0, sep = ",")
    df_rld_total = position_extraction(df_rld_total)
    df_rld_total = dam_id_window_filtering(df_rld_total, chrom, center = center, window = window)
    
    
    
    file_rld_g = f"/datas/nathan/Dam_ID_analysis/results/rld_G{gal_condition}.csv"
    
    df_rld_g = pd.read_csv(file_rld_g, header = 0, sep = ",")
    df_rld_g = position_extraction(df_rld_g)
    df_rld_g = dam_id_window_filtering(df_rld_g, chrom, center = center, window = window)    


    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures(len(conditions_list),1)


    for count, condition in enumerate(conditions_list):
        
        print(condition)



        try:
            df_rld_total[condition]
        except KeyError:
            print(f"{condition} has not been found in the DF")
            continue
                

        subfigs[count].suptitle(condition)
        
        axs = subfigs[count].subplots(1, 1)
     

        # Associating a color scale to the density of the data
        x = np.array(df_rld_total[control])
        y = np.array(df_rld_total[condition])
        
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        
        axs.scatter(x,
                       y,
                       c=z,
                       s=50)
        

    plt.show()



def scatter_dds_cel(gal_condition):
    

    chrom_names = {
    "ref|NC_001133|" : "chrom_1",
    "ref|NC_001134|" : "chrom_2",
    "ref|NC_001135|" : "chrom_3",
    "ref|NC_001136|" : "chrom_4",
    "ref|NC_001137|" : "chrom_5",
    "ref|NC_001138|" : "chrom_6",
    "ref|NC_001139|" : "chrom_7",
    "ref|NC_001140|" : "chrom_8",
    "ref|NC_001141|" : "chrom_9",
    "ref|NC_001142|" : "chrom_10",
    "ref|NC_001143|" : "chrom_11",
    "ref|NC_001144|" : "chrom_12",
    "ref|NC_001145|" : "chrom_13",
    "ref|NC_001146|" : "chrom_14",
    "ref|NC_001147|" : "chrom_15",
    "ref|NC_001148|" : "chrom_16",
    "ref|NC_001224|" : "mitoch",
    }

    

    chrom = "mitoch"


    control = ["srf-3i1-NLS-GFP_L2_rep1",
               "srf-3i1-NLS-GFP_L2_rep2"]

    
    conditions_list = ["srf-3i1-rpb-6_L2_rep1",
                       "srf-3i1-rpb-6_L2_rep2"]


    file_dds_total = f"/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L2/dds_condition.csv"

    df_dds_total = pd.read_csv(file_dds_total, header = 0, sep = ",")
    df_dds_total = position_extraction(df_dds_total)
    #df_dds_total = dam_id_window_filtering(df_dds_total, chrom, center = center, window = window)
    print(len(df_dds_total))



    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures()



    axs = subfigs.subplots(1, 1)
    

    """
    df_dds_total = df_dds_total[(df_dds_total[conditions_list[0]] != 0) 
                                & (df_dds_total[conditions_list[1]] != 0)]
    
    df_dds_total = df_dds_total[(df_dds_total[control[0]] != 0)
                                & (df_dds_total[control[1]] != 0)]
    """

    
    
    # Associating a color scale to the density of the data

    x_line = np.arange(0,df_dds_total[conditions_list[0]].max(), 2)
    y_line = x_line
    
    x = [np.mean([control_1, control_2]) 
        for control_1, control_2 
        in zip(df_dds_total[control[0]],
                df_dds_total[control[1]])]
    
    y = [np.mean([condition_1, condition_2]) 
        for condition_1, condition_2 
        in zip(df_dds_total[conditions_list[0]],
                df_dds_total[conditions_list[1]])]
    
    x_new = list()
    y_new = list()
    
    for mean_x, mean_y in zip(x, y):
        if mean_x >= 1 and mean_y >= 1:
            x_new.append(mean_x)
            y_new.append(mean_y)
        
    x = np.array(x_new)
    y = np.array(y_new)
    
    
    print(np.mean(x))
    print(np.mean(y))

    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    
    axs.scatter(x,
                    y,
                    c=z,
                    s=50)
    
    print(len(x_new))
    
    
    axs.plot(x_line, y_line, color = "black")
    

    """
    axs.scatter(df_dds_total[control][(df_dds_total["chrom"] != "chrom_XII") 
                                    & (df_dds_total["chrom"] != "chrom_III")
                                    & (df_dds_total["chrom"] != "chrom_circular")],
                df_dds_total[condition][(df_dds_total["chrom"] != "chrom_XII") 
                                        & (df_dds_total["chrom"] != "chrom_III")
                                        & (df_dds_total["chrom"] != "chrom_circular")],
                color = "blue",
                label = "other")
    
    

    axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_XII"],
                df_dds_total[condition][df_dds_total["chrom"] == "chrom_XII"],
                color = "green",
                label = "chrom 12")        
    

    axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_III"],
                df_dds_total[condition][df_dds_total["chrom"] == "chrom_III"],
                color = "red",
                label = "chrom 3")
    
    axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_circulaire"],
                df_dds_total[condition][df_dds_total["chrom"] == "chorm_circulaire"],
                color = "yellow",
                label = "Mitochondrial DNA")        
    """
    
    axs.set_ylabel(conditions_list[0])
    axs.set_xlabel(control[0])
    axs.set_yscale("log")
    axs.set_xscale("log")

    plt.legend()
    plt.show()

    return fig



def scatter_raw_cel(gal_condition):
    


    file_control_1 = "/datas/nathan/Dam_ID_analysis/data/counts/C_elegans/srf3/L2/srf-3i1-NLS-GFP_L2_rep1/abundance.tsv"
    file_control_2 = "/datas/nathan/Dam_ID_analysis/data/counts/C_elegans/srf3/L2/srf-3i1-NLS-GFP_L2_rep2/abundance.tsv"
    
    file_test_1 = "/datas/nathan/Dam_ID_analysis/data/counts/C_elegans/srf3/L2/srf-3i1-rpb-6_L2_rep1/abundance.tsv"
    file_test_2 = "/datas/nathan/Dam_ID_analysis/data/counts/C_elegans/srf3/L2/srf-3i1-rpb-6_L2_rep2/abundance.tsv"


    chrom = "mitoch"



    df_control_1 = kallisto_abundance_reformating(file_control_1)
    df_control_2 = kallisto_abundance_reformating(file_control_2)

    print(len(df_control_1))
    
    df_test_1 = kallisto_abundance_reformating(file_test_1)
    df_test_2 = kallisto_abundance_reformating(file_test_2)

    df_control = df_control_1
    
    df_test = df_test_1

    df_control["est_counts"] = (df_control_1["est_counts"] + df_control_2["est_counts"]) / 2

    df_test["est_counts"] = (df_test_1["est_counts"] + df_test_2["est_counts"]) / 2

    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures()



    axs = subfigs.subplots(1, 1)

    df_tot = pd.DataFrame()

    df_tot["est_counts_test"] = df_test["est_counts"]
    df_tot["est_counts_control"] = df_control["est_counts"]

    

    df_tot = df_tot[(df_tot["est_counts_test"] > 1)
                     & (df_tot["est_counts_control"] > 1)]

    print(len(df_tot))
    
    # Associating a color scale to the density of the data

    x_line = np.arange(0,df_test["est_counts"].max(), 2)
    y_line = x_line


    x = np.array(df_tot["est_counts_control"])
    y = np.array(df_tot["est_counts_test"])
    
    print(np.mean(x))
    print(np.mean(y))
    
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    
    axs.scatter(x,
                    y,
                    c=z,
                    s=50)

    
    axs.plot(x_line, y_line, color = "black")
    
    """
    axs.scatter(df_dds_total[control][(df_dds_total["chrom"] != "chrom_XII") 
                                    & (df_dds_total["chrom"] != "chrom_III")
                                    & (df_dds_total["chrom"] != "chrom_circular")],
                df_dds_total[condition][(df_dds_total["chrom"] != "chrom_XII") 
                                        & (df_dds_total["chrom"] != "chrom_III")
                                        & (df_dds_total["chrom"] != "chrom_circular")],
                color = "blue",
                label = "other")
    
    

    axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_XII"],
                df_dds_total[condition][df_dds_total["chrom"] == "chrom_XII"],
                color = "green",
                label = "chrom 12")        
    

    axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_III"],
                df_dds_total[condition][df_dds_total["chrom"] == "chrom_III"],
                color = "red",
                label = "chrom 3")
    
    axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_circulaire"],
                df_dds_total[condition][df_dds_total["chrom"] == "chorm_circulaire"],
                color = "yellow",
                label = "Mitochondrial DNA")        
    """
    
    axs.set_ylabel("test")
    axs.set_xlabel("control")
    axs.set_yscale("log")
    axs.set_xscale("log")

    plt.legend()
    plt.show()

    
    
    return fig



def scatter_dds_sacc(gal_condition):
    

    chrom_names = {
    "ref|NC_001133|" : "chrom_1",
    "ref|NC_001134|" : "chrom_2",
    "ref|NC_001135|" : "chrom_3",
    "ref|NC_001136|" : "chrom_4",
    "ref|NC_001137|" : "chrom_5",
    "ref|NC_001138|" : "chrom_6",
    "ref|NC_001139|" : "chrom_7",
    "ref|NC_001140|" : "chrom_8",
    "ref|NC_001141|" : "chrom_9",
    "ref|NC_001142|" : "chrom_10",
    "ref|NC_001143|" : "chrom_11",
    "ref|NC_001144|" : "chrom_12",
    "ref|NC_001145|" : "chrom_13",
    "ref|NC_001146|" : "chrom_14",
    "ref|NC_001147|" : "chrom_15",
    "ref|NC_001148|" : "chrom_16",
    "ref|NC_001224|" : "mitoch",
    }

    

    chrom = "mitoch"


    control = f"2RG{gal_condition}"

    
    conditions_list = [f"3RG{gal_condition}D1",
                       f"3RG{gal_condition}D2",
                       f"3RG{gal_condition}D3",
                       f"3RG{gal_condition}",]


    file_dds_total = f"/datas/nathan/Dam_ID_analysis/results/dds_G{gal_condition}.csv"

    df_dds_total = pd.read_csv(file_dds_total, header = 0, sep = ",")
    df_dds_total = position_extraction(df_dds_total)
    #df_dds_total = dam_id_window_filtering(df_dds_total, chrom, center = center, window = window)
    print(len(df_dds_total))



    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures(len(conditions_list))

    i = 0

    for condition in conditions_list:
        
        try:
            x = df_dds_total[condition]
        except KeyError:
            continue 

        subfigs[i].suptitle(condition)

        axs = subfigs[i].subplots(1, 1)
        """
        df_dds_total = df_dds_total[(df_dds_total[conditions_list[0]] != 0) 
                                    & (df_dds_total[conditions_list[1]] != 0)]
        
        df_dds_total = df_dds_total[(df_dds_total[control] != 0)
                                    & (df_dds_total[control] != 0)]
        """

        
        
        # Associating a color scale to the density of the data

        x_line = np.arange(0,1e6, 2)
        y_line = x_line
        
        x = np.array(df_dds_total[control])
        y = np.array(df_dds_total[condition])
   
        
        print(np.mean(x))
        print(np.mean(y))

        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        
        axs.scatter(x,
                        y,
                        c=z,
                        s=50)
        
        
        axs.plot(x_line, y_line, color = "black")
        

        """
        axs.scatter(df_dds_total[control][(df_dds_total["chrom"] != "chrom_XII") 
                                        & (df_dds_total["chrom"] != "chrom_III")
                                        & (df_dds_total["chrom"] != "chrom_circular")],
                    df_dds_total[condition][(df_dds_total["chrom"] != "chrom_XII") 
                                            & (df_dds_total["chrom"] != "chrom_III")
                                            & (df_dds_total["chrom"] != "chrom_circular")],
                    color = "blue",
                    label = "other")
        
        

        axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_XII"],
                    df_dds_total[condition][df_dds_total["chrom"] == "chrom_XII"],
                    color = "green",
                    label = "chrom 12")        
        

        axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_III"],
                    df_dds_total[condition][df_dds_total["chrom"] == "chrom_III"],
                    color = "red",
                    label = "chrom 3")
        
        axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_circulaire"],
                    df_dds_total[condition][df_dds_total["chrom"] == "chorm_circulaire"],
                    color = "yellow",
                    label = "Mitochondrial DNA")        
        """
        
        axs.set_ylabel(conditions_list[0])
        axs.set_xlabel(control[0])
        axs.set_yscale("log")
        axs.set_xscale("log")

        i += 1
        
    plt.legend()
    plt.show()

    return fig



def scatter_raw_sacc(gal_condition):
    

    chrom_names = {
    "ref|NC_001133|" : "chrom_1",
    "ref|NC_001134|" : "chrom_2",
    "ref|NC_001135|" : "chrom_3",
    "ref|NC_001136|" : "chrom_4",
    "ref|NC_001137|" : "chrom_5",
    "ref|NC_001138|" : "chrom_6",
    "ref|NC_001139|" : "chrom_7",
    "ref|NC_001140|" : "chrom_8",
    "ref|NC_001141|" : "chrom_9",
    "ref|NC_001142|" : "chrom_10",
    "ref|NC_001143|" : "chrom_11",
    "ref|NC_001144|" : "chrom_12",
    "ref|NC_001145|" : "chrom_13",
    "ref|NC_001146|" : "chrom_14",
    "ref|NC_001147|" : "chrom_15",
    "ref|NC_001148|" : "chrom_16",
    "ref|NC_001224|" : "mitoch",
    }

    file_control = f"/datas/nathan/Dam_ID_analysis/data/counts/G{gal_condition}/2RG{gal_condition}/abundance.tsv"


    files_condition_list = [f"/datas/nathan/Dam_ID_analysis/data/counts/G{gal_condition}/3RG{gal_condition}D1/abundance.tsv",
                            f"/datas/nathan/Dam_ID_analysis/data/counts/G{gal_condition}/3RG{gal_condition}D2/abundance.tsv",
                            f"/datas/nathan/Dam_ID_analysis/data/counts/G{gal_condition}/3RG{gal_condition}D3/abundance.tsv",
                            f"/datas/nathan/Dam_ID_analysis/data/counts/G{gal_condition}/3RG{gal_condition}/abundance.tsv"]


    conditions_names_list = [f"3RG{gal_condition}D1",
                             f"3RG{gal_condition}D2",
                             f"3RG{gal_condition}D3",
                             f"3RG{gal_condition}"]


    df_control = kallisto_abundance_reformating(file_control)
    

    fig = plt.figure(constrained_layout = True)


    subfigs = fig.subfigures(len(files_condition_list))

    i = 0
    
    for file_condition, condition in zip(files_condition_list, conditions_names_list):

    

        if exists(file_condition) == False:
            continue

        
        df_condition = kallisto_abundance_reformating(file_condition)

        axs = subfigs[i].subplots(1, 1)

        subfigs[i].suptitle(condition)
        
        # Associating a color scale to the density of the data

        x_line = np.arange(0, 1e6, 2)
        y_line = x_line

        df_fus = pd.DataFrame()
        
        df_fus["control"] = df_control["est_counts"]
        df_fus["condition"] = df_condition["est_counts"]
        
        df_fus = df_fus[(df_fus["condition"] >= 1)
                        & (df_fus["control"] >= 1)]
        
        df_fus = df_fus

        x = np.array(df_fus["control"])
        y = np.array(df_fus["condition"])

        
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        
        axs.scatter(x,
                        y,
                        c=z,
                        s=50)

        
        axs.plot(x_line, y_line, color = "black")
        
        """
        axs.scatter(df_dds_total[control][(df_dds_total["chrom"] != "chrom_XII") 
                                        & (df_dds_total["chrom"] != "chrom_III")
                                        & (df_dds_total["chrom"] != "chrom_circular")],
                    df_dds_total[condition][(df_dds_total["chrom"] != "chrom_XII") 
                                            & (df_dds_total["chrom"] != "chrom_III")
                                            & (df_dds_total["chrom"] != "chrom_circular")],
                    color = "blue",
                    label = "other")
        
        

        axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_XII"],
                    df_dds_total[condition][df_dds_total["chrom"] == "chrom_XII"],
                    color = "green",
                    label = "chrom 12")        
        

        axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_III"],
                    df_dds_total[condition][df_dds_total["chrom"] == "chrom_III"],
                    color = "red",
                    label = "chrom 3")
        
        axs.scatter(df_dds_total[control][df_dds_total["chrom"] == "chrom_circulaire"],
                    df_dds_total[condition][df_dds_total["chrom"] == "chorm_circulaire"],
                    color = "yellow",
                    label = "Mitochondrial DNA")        
        """
        
        axs.set_ylabel(f"test")
        axs.set_xlabel("control")
        axs.set_yscale("log")
        axs.set_xscale("log")

        i += 1


    plt.legend()
    plt.show()

    
    
    return fig



def scatter_dds_binned_reads(gal_condition, epsilon):

    chrom = "chrom_3"


    conditions_list = [f"3RG{gal_condition}D1",
                       f"3RG{gal_condition}D2",
                       f"3RG{gal_condition}D3",
                       f"3RG{gal_condition}",]

    
    control = f"2RG{gal_condition}"


    file_dds_total = f"/datas/nathan/Dam_ID_analysis/results/dds_G{gal_condition}.csv"

    df_dds_total = pd.read_csv(file_dds_total, header = 0, sep = ",")
    df_dds_total = position_extraction(df_dds_total)
    #df_dds_total = dam_id_window_filtering(df_dds_total, chrom, center = center, window = window)

    df_control = binning_deseq_reads(df_dds_total, control, epsilon, chrom = "all")   


    fig = plt.figure(constrained_layout = True,
                     figsize=[50,50])
    fig.suptitle(chrom)

    subfigs = fig.subfigures(len(conditions_list),1)


    for count, condition in enumerate(conditions_list):
        

        try:
            df_dds_total[condition]
        except KeyError:
            print(f"{condition} has not been found in the DF")
            continue
                

        df_condition = binning_application_deseq(df_dds_total, df_control, condition)


        subfigs[count].suptitle(condition)
        
        axs = subfigs[count].subplots(1, 1)



        # Associating a color scale to the density of the data


        
        x_line = np.arange(0,df_condition["value"].max(), 2)
        y_line = x_line
        
        axs.scatter(df_control["value"][(df_control["chrom"] != "chrom_12") 
                                           & (df_control["chrom"] != "chrom_3")],
                    df_condition["value"][(df_condition["chrom"] != "chrom_12") 
                                         & (df_condition["chrom"] != "chrom_3")],
                    color = "blue",
                    label = "other")
        
        axs.plot(x_line, y_line, color = "black")

        axs.scatter(df_control["value"][df_control["chrom"] == "chrom_12"],
                    df_condition["value"][df_condition["chrom"] == "chrom_12"],
                    color = "green",
                    label = "chrom 12")        
        

        axs.scatter(df_control["value"][df_control["chrom"] == "chrom_3"],
                    df_condition["value"][df_condition["chrom"] == "chrom_3"],
                    color = "red",
                    label = "chrom 3")
        
        
        axs.set_ylabel(condition)
        axs.set_xlabel(control)
        axs.set_yscale("log")
        axs.set_xscale("log")

    plt.legend()


    return fig  



def scatter_dds_binned_gatc(gal_condition, epsilon):
    
    chrom = "chrom_3"


    conditions_list = [f"3RG{gal_condition}D1",
                       f"3RG{gal_condition}D2",
                       f"3RG{gal_condition}D3",
                       f"3RG{gal_condition}",]

    
    control = f"2RG{gal_condition}"


    file_dds_total = f"/datas/nathan/Dam_ID_analysis/results/dds_G{gal_condition}.csv"

    df_dds_total = pd.read_csv(file_dds_total, header = 0, sep = ",")
    df_dds_total = position_extraction(df_dds_total)
    #df_dds_total = dam_id_window_filtering(df_dds_total, chrom, center = center, window = window)

    df_control = binning_deseq_gatc(df_dds_total, control, epsilon, chrom = "all")   


    fig = plt.figure(constrained_layout = True,
                     figsize=[50,50])
    fig.suptitle(chrom)

    subfigs = fig.subfigures(len(conditions_list),1)


    for count, condition in enumerate(conditions_list):
        

        try:
            df_dds_total[condition]
        except KeyError:
            print(f"{condition} has not been found in the DF")
            continue
                

        df_condition = binning_deseq_gatc(df_dds_total, condition, epsilon, chrom = "all")


        subfigs[count].suptitle(condition)
        
        axs = subfigs[count].subplots(1, 1)



        # Associating a color scale to the density of the data


        
        x_line = np.arange(0,df_condition["value"].max(), 2)
        y_line = x_line
        
        axs.scatter(df_control["value"][(df_control["chrom"] != "chrom_12") 
                                           & (df_control["chrom"] != "chrom_3")],
                    df_condition["value"][(df_condition["chrom"] != "chrom_12") 
                                         & (df_condition["chrom"] != "chrom_3")],
                    color = "blue",
                    label = "other")
        
        axs.plot(x_line, y_line, color = "black")

        axs.scatter(df_control["value"][df_control["chrom"] == "chrom_12"],
                    df_condition["value"][df_condition["chrom"] == "chrom_12"],
                    color = "green",
                    label = "chrom 12")        
        

        axs.scatter(df_control["value"][df_control["chrom"] == "chrom_3"],
                    df_condition["value"][df_condition["chrom"] == "chrom_3"],
                    color = "red",
                    label = "chrom 3")
        
        
        axs.set_ylabel(condition)
        axs.set_xlabel(control)
        axs.set_yscale("log")
        axs.set_xscale("log")

    plt.legend()
    plt.show()

    return fig


def raw_deseq_comp():
    
    chrom_names_celegans = {
    "ENA|BX284601|BX284601.5" : "chromosome_I",
    "ENA|BX284602|BX284602.5" : "chromosome_II",
    "ENA|BX284603|BX284603.4" : "chromosome_III",
    "ENA|BX284604|BX284604.4" : "chromosome_IV",
    "ENA|BX284605|BX284605.5" : "chromosome_V",
    "ENA|BX284606|BX284606.5" : "chromosome_X"
    }




    control = ["srf-3i1-NLS-GFP_L2_rep1",
               "srf-3i1-NLS-GFP_L2_rep2"]

    
    conditions_list = ["srf-3i1-rpb-6_L2_rep1",
                       "srf-3i1-rpb-6_L2_rep2"]


    file_dds_total = f"/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L4/dds_condition.csv"
    df_dds_total = pd.read_csv(file_dds_total, header = 0, sep = ",")
    df_dds_total = position_extraction(df_dds_total)




    file_raw = "/datas/nathan/Dam_ID_analysis/data/counts/C_elegans/srf3/L4/srf-3i1-NLS-GFP_L4_rep1/abundance.tsv"
    df_raw = kallisto_abundance_extraction(file_raw, chrom_names_celegans)

    
    file_results = "/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L4/results_condition.csv"
    df_results = pd.read_csv(file_results, header = 0, sep = ",")
    df_results = position_extraction(df_results)

    
    

    
    df_results = df_results.sort_values(["length", "start"])
    df_raw = df_raw.sort_values(["length", "start"])
    df_dds_total = df_dds_total.sort_values(["length", "start"])


    df_results.reset_index(drop=True, inplace=True)
    df_raw.reset_index(drop=True, inplace=True)
    df_dds_total.reset_index(drop=True, inplace=True)
    
    print(df_raw)
    print(df_dds_total)
    
    
    x = df_raw["est_counts"]
    y = df_dds_total["srf-3i1-NLS-GFP_L4_rep1"]
    c = ["red" if pval < 0.05
         else "blue"
         for pval in df_results["pvalue"]]

    x_line = np.arange(0, x.max(), 2)
    
    y_line = x_line
    
    x_y = pd.DataFrame()
    x_y["raw"] = df_raw["est_counts"]
    x_y["DESe2"] = df_dds_total["srf-3i1-NLS-GFP_L4_rep1"]
    x_y["diff"] = x_y["DESe2"] - x_y["raw"]
    print(x_y["diff"].nsmallest(10))
    

    plt.scatter(x, 
                y,
                color = c)

    plt.plot(x_line, y_line, color = "black")



    plt.show()


#####################################################################
#                                                                   #
#                           Other                                   #
#                                                                   #
#####################################################################


def log2_dds(gal_condition):


    chrom = "chrom_3"
    center = 0
    window = 0

    i = 0

    conditions_list = [f"3RG{gal_condition}D1",
                    f"3RG{gal_condition}D2",
                    f"3RG{gal_condition}D3"]
    
    control =  f"2RG{gal_condition}"
    
    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures(len(conditions_list),1)




    for condition in conditions_list:


        file_dds = "/datas/nathan/Dam_ID_analysis/results/dds_total.csv"

        df_dds = pd.read_csv(file_dds, header = 0, sep = ",")


        df_dds = position_extraction(df_dds)


        df_dds = dam_id_window_filtering(df_dds, chrom, center = center, window = window)


        try:
            df_dds[condition]
        except KeyError:
            print(f"{condition} has not been found in the DF")
            continue
        
        
        subfigs[i].suptitle(f"{condition}")
        
        axs = subfigs[i].subplots(1, 2,
                                sharex = True,
                                sharey = True)
        
        axs[0].bar(df_dds["start"] + (df_dds["length"]/2),
                    height = np.log2(df_dds[condition] / df_dds[control]) ,
                    width = df_dds["stop"] - df_dds["start"],
                    color = ["red" 
                             if log_ratio < 0 
                             else "blue" 
                             for log_ratio in np.log2(df_dds[condition] / df_dds[control])])

        axs[0].set_title(f"all\nmean = {round(df_dds[condition].mean(), 2)}")
        axs[0].set_ylabel("rlog(ratio)")
        axs[0].set_xlabel("postions")

        axs[0]
        
        file_dds = f"/datas/nathan/Dam_ID_analysis/results/dds_G{gal_condition}.csv"


        df_dds = pd.read_csv(file_dds, header = 0, sep = ",")


        df_dds = position_extraction(df_dds)


        df_dds = dam_id_window_filtering(df_dds, chrom, center = center, window = window)


        axs[1].bar(df_dds["start"] + (df_dds["length"]/2),
                    height = np.log2(df_dds[condition] / df_dds[control]),
                    width = df_dds["stop"] - df_dds["start"],
                    color = ["red" 
                            if log_ratio < 0 
                            else "blue" 
                            for log_ratio in np.log2(df_dds[condition] / df_dds[control])])

        axs[1].set_title(f"Gal{gal_condition} only\nmean = {round(df_dds[condition].mean(), 2)}")
        axs[1].set_ylabel("rlog(ratio)")
        axs[1].set_xlabel("postions")

        i += 1


        print((df_dds.nlargest(5, condition)))
    
    plt.show()

    
    
def multi_conditions_analsyis(gal_condition):

    chrom = "chrom_3"

    window = 0

    center = 0

    i = 0

    conditions_list = [f"3RG{gal_condition}D1",
                       f"3RG{gal_condition}D2", 
                       f"3RG{gal_condition}D3",
                       f"2RG{gal_condition}"]


    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures(3,1)
    
    file_rld_total = "/datas/nathan/Dam_ID_analysis/results/rld_total.csv"
    df_rld_total = pd.read_csv(file_rld_total, header = 0, sep = ",")
    df_rld_total = position_extraction(df_rld_total)
    df_rld_total = dam_id_window_filtering(df_rld_total, 
                                           chrom, 
                                           center = center, 
                                           window = window)


    file_rld_conditions = "/datas/nathan/Dam_ID_analysis/results/rld_conditions.csv"
    df_rld_conditions = pd.read_csv(file_rld_conditions, header = 0, sep = ",")
    df_rld_conditions = position_extraction(df_rld_conditions)
    df_rld_conditions = dam_id_window_filtering(df_rld_conditions, 
                                                chrom, 
                                                center = center, 
                                                window = window)


    print(df_rld_conditions[(df_rld_conditions["start"] > 90000) & (df_rld_conditions["stop"] < 95000)])
    
    max_value = 0
    min_value = pow(10,100)
    for condition in conditions_list:
        
        if df_rld_conditions[condition].max() > max_value:
            max_value = df_rld_conditions[condition].max()
            
        if df_rld_total[condition].max() > max_value:
            max_value = df_rld_total[condition].max()
        
        """
        df_rld_conditions = df_rld_conditions.dropna()
        df_rld_total = df_rld_total.dropna()
        
        if (df_rld_conditions[condition] - df_rld_conditions[control]).max() > max_value:
            max_value = (df_rld_conditions[condition] - df_rld_conditions[control]).max() 

        
        if (df_rld_total[condition] - df_rld_total[control]).max() > max_value:
            max_value = (df_rld_total[condition] - df_rld_total[control]).max()
    

        if (df_rld_conditions[condition] - df_rld_conditions[control]).min() < min_value:
            min_value = (df_rld_conditions[condition] - df_rld_conditions[control]).min()
            
        if (df_rld_total[condition] - df_rld_total[control]).min() < min_value:
            min_value = (df_rld_total[condition] - df_rld_total[control]).min()        
        """
    for condition in conditions_list:
        
                
        try:
            df_rld_conditions[condition]
        except KeyError:
            print(f"{condition} has not been found in the DF")
            continue
                

        subfigs[i].suptitle(condition)
        
        axs = subfigs[i].subplots(1, 2)


        axs[0].bar(df_rld_conditions["start"] + (df_rld_conditions["length"]/2),
                    height = df_rld_conditions[condition] ,
                    width = df_rld_conditions["stop"] - df_rld_conditions["start"],
                    color = ["red" if log_ratio < 0 
                             else "blue" 
                             for log_ratio in df_rld_conditions[condition] - df_rld_conditions[f"2RG{gal_condition}"]])

        axs[0].set_title(f"conditions\nmean = {round(df_rld_conditions[condition].mean(), 2)}")
        axs[0].set_ylabel("rlog(ratio)")
        axs[0].set_xlabel("postions")

        axs[0].set_ylim(0, max_value*1.2)

        axs[1].bar(df_rld_total["start"] + (df_rld_total["length"]/2),
                    height = df_rld_total[condition] ,
                    width = df_rld_total["stop"] - df_rld_total["start"],
                    color = ["red" 
                             if log_ratio < 0 
                             else "blue" 
                             for log_ratio in df_rld_total[condition] - df_rld_total[f"2RG{gal_condition}"]])

        axs[1].set_title(f"all\nmean = {round(df_rld_total[condition].mean(), 2)}")
        axs[1].set_ylabel("rlog(ratio)")
        axs[1].set_xlabel("postions")
        
        axs[1].set_ylim(min_value*1.2, max_value*1.2)

        i += 1


    plt.show()
    
    

def reads_distribution_analysis(gal_condition, chrom = None):
    
    chrom_names = {
    "ref|NC_001133|" : "chrom_1",
    "ref|NC_001134|" : "chrom_2",
    "ref|NC_001135|" : "chrom_3",
    "ref|NC_001136|" : "chrom_4",
    "ref|NC_001137|" : "chrom_5",
    "ref|NC_001138|" : "chrom_6",
    "ref|NC_001139|" : "chrom_7",
    "ref|NC_001140|" : "chrom_8",
    "ref|NC_001141|" : "chrom_9",
    "ref|NC_001142|" : "chrom_10",
    "ref|NC_001143|" : "chrom_11",
    "ref|NC_001144|" : "chrom_12",
    "ref|NC_001145|" : "chrom_13",
    "ref|NC_001146|" : "chrom_14",
    "ref|NC_001147|" : "chrom_15",
    "ref|NC_001148|" : "chrom_16",
    "ref|NC_001224|" : "mitoch",
    }



    center = 0
    window = 0

    i = 0

    conditions_list = [f"3RG{gal_condition}D1",
                       f"3RG{gal_condition}D2",
                       f"3RG{gal_condition}D3",
                       f"2RG{gal_condition}"]
    
    
    fig = plt.figure(constrained_layout = True)
    fig.suptitle("full genome without 0s")

    subfigs = fig.subfigures(len(conditions_list), 1,)

    file_dds_condition = "/datas/nathan/Dam_ID_analysis/results/dds_conditions.csv"

    df_dds_condition = pd.read_csv(file_dds_condition, header = 0, sep = ",")
    df_dds_condition = position_extraction(df_dds_condition)
    
    
    
    file_dds_all = "/datas/nathan/Dam_ID_analysis/results/dds_total.csv"

    df_dds_all = pd.read_csv(file_dds_all, header = 0, sep = ",")
    df_dds_all = position_extraction(df_dds_all)    
    
    

    
    
    
    """
    max_value = 0
    
    for condition in conditions_list:
        
        path_file_condition = f"/datas/nathan/Dam_ID_analysis/data/yes/{condition}/abundance.tsv"
        
        path_file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/test/yes/sites.bed"
        
        df_raw = kallisto_outup_reformating(path_file_condition, path_file_sites, chrom_names = chrom_names)
        
        if df_dds[condition].max() > max_value:
            max_value = df_dds[condition].max()
        
        if df_raw["est_counts"].max() > max_value:
            max_value = df_raw["est_counts"].max()
    """
        
        

    for count, condition in enumerate(conditions_list):
        
        try:
            df_dds_condition[condition]
        except KeyError:
            print(f"{condition} has not been found in the DF")
            continue



            
        axs = subfigs[count].subplots(1, 3,
                                      sharex = True,
                                      sharey = True)
        
        
        path_file_condition = f"/datas/nathan/Dam_ID_analysis/data/yes/{condition}/abundance.tsv"
        
        path_file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/test/yes/sites.bed"
        
        df_raw = kallisto_outup_reformating(path_file_condition, path_file_sites, chrom_names = chrom_names)
        
        if chrom != None:
            df_dds_condition = dam_id_window_filtering(df_dds_condition, chrom = chrom, window = window, center = center)
            df_dds_all = dam_id_window_filtering(df_dds_all, chrom = chrom, window = window, center = center)
            df_raw = dam_id_window_filtering(df_raw, chrom = chrom, window = window, center = center)

        df_dds_all = df_dds_all[df_dds_all[condition] != 0]
        df_dds_condition = df_dds_condition[df_dds_condition[condition] !=0]
        df_raw = df_raw[df_raw["est_counts"] != 0]


        axs[0].hist(df_dds_condition[condition],
                    bins = 5000)
        
        axs[1].hist(df_dds_all[condition],
                    bins = 5000)
        
        axs[0].set_xlabel("count/frag")
        axs[0].set_ylabel("population")
        axs[0].set_title(f"conditions\nmean={round(df_dds_condition[condition].mean(), 2)}")
        #axs[0].set_ylim(0,1.2*max_value)
        axs[0].set_xlim(0, 2000)
        
        axs[1].set_xlabel("count/frag")
        axs[1].set_ylabel("population")
        axs[1].set_title(f"all\nmean={round(df_dds_all[condition].mean(), 2)}")
        #axs[1].set_ylim(0,1.2*max_value)
        axs[1].set_xlim(0, 2000)
        
        
        
        axs[2].hist(df_raw["est_counts"],
                    bins = 5000)
        axs[2].set_xlabel("count/frag")
        axs[2].set_ylabel("population")
        axs[2].set_title(f"raw\nmean={round(df_raw['est_counts'].mean(), 2)}")
        #axs[2].set_ylim(0,1.2*max_value)
        axs[2].set_xlim(0, 2000)
        
        
        subfigs[count].suptitle(f"{condition}")
        
    plt.show()


def lfc_deseq_plotting(tissue, stage, ax):


    file_results = f"/datas/nathan/Dam_ID_analysis/results/C_elegans/{tissue}/{stage}/results_condition.csv"
    
    df_results = pd.read_csv(file_results, header = 0, sep = ",")
    
    df_results = position_extraction(df_results)

    df_plot = dam_id_window_filtering(df_results, "chromosome_IV")
    
    df_plot = df_plot[df_plot["log2FoldChange"] == df_plot["log2FoldChange"]]
    
    ax.bar(df_plot["start"] + (df_plot["length"]/2),
            height = df_plot["log2FoldChange"],
            width = df_plot["stop"] - df_plot["start"],
            color = ["red"
                        if lfc < -1
                        else "blue"
                     if lfc > 1
                        else "grey"
                        for lfc in df_plot["log2FoldChange"]])


    ax.set_title(f"{tissue} {stage}")
    ax.set_ylabel("LFC test/control")


    return ax


#####################################################################
#                                                                   #
#                              Pvalue                               #
#                                                                   #
#####################################################################


def volcano_plot_sleuth():

    
    file_sleuth = "/datas/nathan/Dam_ID_analysis/results/sleuth_G3.csv"
    
    df_sleuth = sleuth_output_reformating(file_sleuth)
    
    print(df_sleuth)



    df_circles = df_sleuth[(df_sleuth["chrom"] != "chrom_circular")
                           & (df_sleuth["chrom"] != "chrom_III")
                           & (df_sleuth["chrom"] != "chrom_XII")]    
    
    plt.scatter(df_circles["b"],
                -np.log10(df_circles["pval"]),
                color = ["red" 
                         if b < -1 and pvalue > -np.log10(0.05)
                         else "green"
                         if b > 1 and pvalue > -np.log10(0.05)
                         else "grey"
                         for b, pvalue in zip(df_circles["b"], -np.log10(df_circles["pval"]))],
                label = "other")
    
    
    
    df_star = df_sleuth[df_sleuth["chrom"] == "chrom_III"]
    
    plt.scatter(df_star["b"],
                -np.log10(df_star["pval"]),
                color = ["red" 
                         if b < -1 and pvalue > -np.log10(0.05)
                         else "green"
                         if b > 1 and pvalue > -np.log10(0.05)
                         else "grey"
                         for b, pvalue in zip(df_star["b"],
                                              -np.log10(df_star["pval"]))],
                marker = "*",
                edgecolors = "black",
                label = "chrom III")
    
    
    
    df_triangle = df_sleuth[df_sleuth["chrom"] == "chrom_circular"]
    
    plt.scatter(df_triangle["b"],
                -np.log10(df_triangle["pval"]),
                color = ["red" 
                         if b < -1 and pvalue > -np.log10(0.05)
                         else "green"
                         if b > 1 and pvalue > -np.log10(0.05)
                         else "grey"
                         for b, pvalue in zip(df_triangle["b"],
                                              -np.log10(df_triangle["pval"]))],
                marker = "^",
                edgecolors = "black",
                label = "ADN mitochondiral")
    
    
    
    df_square = df_sleuth[df_sleuth["chrom"] == "chrom_XII"]
    
    plt.scatter(df_square["b"],
                -np.log10(df_square["pval"]),
                color = ["red" 
                         if b < -1 and pvalue > -np.log10(0.05)
                         else "green"
                         if b > 1 and pvalue > -np.log10(0.05)
                         else "grey"
                         for b, pvalue in zip(df_square["b"],
                                              -np.log10(df_square["pval"]))],
                marker = "s",
                edgecolors = "black",
                label = "chrom XII")
             
             
                

                
    plt.xlabel("log2(ratio)")
    plt.ylabel("-log10(pvalue)")
    
    plt.legend()
    plt.show()



def pvalue_plotting_nlargest():


    file_sleuth = "/datas/nathan/Dam_ID_analysis/results/sleuth_G3.csv"
    
    df_sleuth = sleuth_output_reformating(file_sleuth)
    
    print(df_sleuth.nsmallest(5, "pval"))
    
    fig = plt.figure(constrained_layout = True)
    
    subfigs = fig.subfigures(4,1)
    
    
    for index, row in df_sleuth.nsmallest(4, "pval").iterrows():
        
        axs = subfigs[index].subplots()
        
        center = row["start"] + (row["length"]/2)
        
        window = 50000
        
        df_plot = dam_id_window_filtering(df_sleuth,
                                          row["chrom"],
                                          window = window,
                                          center = center)
        
        axs.bar(df_plot["start"] + (df_plot["length"]/2),
                height = -np.log10(df_plot["pval"]),
                width = df_plot["stop"] - df_plot["start"],
                color = ["red"
                         if b < -1
                         else "blue"
                         if b > 1
                         else "grey"
                         for b in df_plot["b"]])

        x = np.arange(center - window/2,
                      center + window/2,
                      1)
        
        y = [-np.log10(0.05) for i in x]
        
        plt.plot(x,
                 y,
                 "--",
                 color = "black",
                 label = "pvalue = 0.05")
        
        axs.set_title(row["chrom"])
        axs.set_ylabel("-log10(pvalue)")
        axs.set_xlabel("postions")

    plt.show()



def qvalue_plotting_sleuth(file_sleuth):
    
    
    df_sleuth = sleuth_output_reformating(file_sleuth)

    df_plot = dam_id_window_filtering(df_sleuth, 
                                      "chrom_III",
                                      window = 20000,
                                      center = 9.62e6)
    
    plt.bar(df_plot["start"] + (df_plot["length"]/2),
            height = -np.log10(df_plot["qval"]),
            width = df_plot["stop"] - df_plot["start"],
            color = ["red"
                        if b < -1
                        else "blue"
                     if b > 1
                        else "grey"
                        for b in df_plot["b"]])

    x = np.arange(0,
                  df_plot["stop"].max(),
                  1)
    
    y = [-np.log10(0.05) for i in x]
    
    plt.plot(x,
                y,
                "--",
                color = "black",
                label = "pvalue = 0.05")
    
    plt.title("chrom_III")
    plt.ylabel("-log10(pvalue)")
    plt.xlabel("postions")
    plt.xlim(0, df_plot["stop"].max()+1000)


    plt.show()


def lfc_plotting_sleuth(axis, file_sleuth):
    
    
    chrom_names_celegans = {
    "ENA|BX284601|BX284601.5" : "chromosome_I",
    "ENA|BX284602|BX284602.5" : "chromosome_II",
    "ENA|BX284603|BX284603.4" : "chromosome_III",
    "ENA|BX284604|BX284604.4" : "chromosome_IV",
    "ENA|BX284605|BX284605.5" : "chromosome_V",
    "ENA|BX284606|BX284606.5" : "chromosome_X"
    }

    
    if file_sleuth == "/datas/nathan/Dam_ID_analysis/results/C_elegans/dpy7/L2/sleuth_results.csv":
        df_sleuth = sleuth_output_reformating(file_sleuth, chrom_names_celegans)

    else:
        df_sleuth = sleuth_output_reformating(file_sleuth)
    
    df_plot = dam_id_window_filtering(df_sleuth,
                                      "chromosome_IV",
                                      window = 60000,
                                      center = 9.62e6)
    
    df_plot["b"] = [0
                    if lfc == np.nan
                    else lfc
                    for lfc in df_plot["b"]]
    
    axis.bar(df_plot["start"] + (df_plot["length"]/2),
            height = df_plot["b"],
            width = df_plot["stop"] - df_plot["start"],
            color = ["red"
                     if b < -1
                     else "blue"
                     if b > 1
                     else "grey"
                     for b in df_plot["b"]])


    
    axis.set_xlim(df_plot["start"].min(), df_plot["stop"].max()+1000)
    axis.set_ylabel("lfc")
    axis.set_xlabel("postions")


    return axis



def pvalue_plotting_deseq2(tissue, stage, ax):
    

    file_results = f"/datas/nathan/Dam_ID_analysis/results/C_elegans/{tissue}/{stage}/results_condition.csv"
    
    df_results = pd.read_csv(file_results, header = 0, sep = ",")
    
    df_results = position_extraction(df_results)

    center = 6861000
    
    window = 20000

    df_plot = dam_id_window_filtering(df_results, "chromosome_II", center = 6861000, window = 20000)
    
    ax.bar(df_plot["start"] + (df_plot["length"]/2),
            height = [-np.log10(padj)
                      if lfc > 0
                      else np.log10(padj)
                      if lfc < 0
                      else 0
                      for padj, lfc 
                      in zip(df_plot["padj"], df_plot["log2FoldChange"])],
            width = df_plot["stop"] - df_plot["start"],
            color = ["red"
                        if lfc < -1
                        else "blue"
                     if lfc > 1
                        else "grey"
                        for lfc in df_plot["log2FoldChange"]])

    x_top = np.arange(center - window/2,
                  center + window/2,
                  1)
    
    y_top = [-np.log10(0.05) for i in x_top]
    
    ax.plot(x_top,
             y_top,
             "--",
             color = "black",
             label = "pvalue = 0.05")

    x_bot = np.arange(center - window/2,
                  center + window/2,
                  1)
    
    y_bot = [np.log10(0.05) for i in x_top]
    
    ax.plot(x_bot,
             y_bot,
             "--",
             color = "black")

    x_middle = np.arange(center - window/2,
                         center + window/2,
                         1)
    
    y_middle = [0 for i in x_middle]
    
    ax.plot(x_middle, y_middle)
    
    ax.set_title(f"{tissue} {stage}")
    ax.set_ylabel("-log10(pvalue)")
    ax.set_xlabel("postions")
    ax.set_xlim(center - window/2,
                  center + window/2,)



    return ax



def qval_plotting_sleuth(axis, file_sleuth):
    
    
    chrom_names_celegans = {
    "ENA|BX284601|BX284601.5" : "chromosome_I",
    "ENA|BX284602|BX284602.5" : "chromosome_II",
    "ENA|BX284603|BX284603.4" : "chromosome_III",
    "ENA|BX284604|BX284604.4" : "chromosome_IV",
    "ENA|BX284605|BX284605.5" : "chromosome_V",
    "ENA|BX284606|BX284606.5" : "chromosome_X"
    }



    center = 9.62e6
    
    window = 60000
    
    if file_sleuth == "/datas/nathan/Dam_ID_analysis/results/C_elegans/dpy7/L2/sleuth_results.csv":
        df_sleuth = sleuth_output_reformating(file_sleuth, chrom_names_celegans)

    else:
        df_sleuth = sleuth_output_reformating(file_sleuth)
    
    df_plot = dam_id_window_filtering(df_sleuth,
                                      "chromosome_IV",
                                      window = 60000,
                                      center = 9.62e6)
    
    df_plot["qval"] = [1
                       if qval == np.nan 
                        or qval == 0
                       else qval
                       for qval in df_plot["qval"]]
    
    axis.bar(df_plot["start"] + (df_plot["length"]/2),
            height = -np.log10(df_plot["qval"]),
            width = df_plot["stop"] - df_plot["start"],
            color = ["red"
                     if b < -1
                     else "blue"
                     if b > 1
                     else "grey"
                     for b in df_plot["b"]])



    x_top = np.arange(center - window/2,
                  center + window/2,
                  1)
    
    y_top = [-np.log10(0.05) for i in x_top]
    
    axis.plot(x_top,
             y_top,
             "--",
             color = "black",
             label = "pvalue = 0.05")

    x_bot = np.arange(center - window/2,
                  center + window/2,
                  1)
    
    y_bot = [np.log10(0.05) for i in x_top]
    
    axis.plot(x_bot,
             y_bot,
             "--",
             color = "black")
  
    axis.set_xlim(df_plot["start"].min(), df_plot["stop"].max()+1000)
    axis.set_ylabel("lfc")
    axis.set_xlabel("postions")


    return axis



files_list = [
    "/datas/nathan/Dam_ID_analysis/results/C_elegans/dpy7/L2/sleuth_results.csv",
    "/datas/nathan/Dam_ID_analysis/results/C_elegans/10_100/results_sleuth_dpy7L2.csv",
    "/datas/nathan/Dam_ID_analysis/results/C_elegans/100_2000/results_sleuth_dpy7L2.csv"
]

fig, axis = plt.subplots(3,
                         1,
                         sharex=True,
                         sharey=True)

for i, file in enumerate(files_list):
    axis[i] = qval_plotting_sleuth(axis[i], file)
    


axis[0].set_title("norm unbinned")
axis[1].set_title("binned 10_100")
axis[2].set_title("binned 100_2000")
fig.suptitle("elt1")

plt.show()











"""


fig, ax = plt.subplots(4,1,
                       sharex= True,
                       sharey= True,
                       figsize=(30,30))


tissue_list = ["srf3", "dpy7"]

stage_list = ["L2", "L4"]

i = 0


for tissue in tissue_list:
    
    for stage in stage_list:

        
        ax[i] = pvalue_plotting_deseq2(tissue, stage, ax[i])
        
        i += 1

fig.suptitle("aff-1")

plt.show()

"""







"""
epsilon = 5
gal_condition = 3
scatter_dds_binned_gatc(gal_condition, epsilon)

for i in range(1,4):
    
    gal_condition = i
    
    fig = scatter_dds_binned_gatc(gal_condition, epsilon)
    path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/figures/DESeq2/scatters/"
            +"scatter_reads_"
            + str(epsilon)
            + "_bins_G"
            + str(gal_condition)
            + ".eps")
    plt.savefig(path, format = "eps")
    plt.close()
"""
