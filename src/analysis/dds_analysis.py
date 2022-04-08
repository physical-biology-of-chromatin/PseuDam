import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

from packages.Dam_ID_utilities import position_extraction, kallisto_outup_reformating, dam_id_window_filtering


def rlog_comparison():
    
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

    conditions_list = ["3RG3D1", "3RG3D2", "3RG3D3"]

    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures(3,1)

    for condition in conditions_list:


        subfigs[i].suptitle(condition)
        
        axs = subfigs[i].subplots(1, 2)



        file_rld = "/datas/nathan/Dam_ID_analysis/rld_total.csv"



        df_rld = pd.read_csv(file_rld, header = 0, sep = ",")



        df_rld = position_extraction(df_rld)

        df_rld = dam_id_window_filtering(df_rld, chrom, center = center, window = window)


        axs[0].bar(df_rld["start"] + (df_rld["length"]/2),
                    height = df_rld[condition] - df_rld["2RG3"],
                    width = df_rld["stop"] - df_rld["start"],
                    color = ["red" if log_ratio < 0 else "blue" for log_ratio in df_rld[condition] - df_rld["2RG3"]])

        axs[0].set_title(f"all\nmean = {round(df_rld[condition].mean(), 2)}")
        axs[0].set_ylabel("rlog(ratio)")
        axs[0].set_xlabel("postions")

        
        
        file_dds = "/datas/nathan/Dam_ID_analysis/dds_G3.csv"

        file_rld = "/datas/nathan/Dam_ID_analysis/rld_G3.csv"

        df_dds = pd.read_csv(file_dds, header = 0, sep = ",")

        df_rld = pd.read_csv(file_rld, header = 0, sep = ",")


        df_dds = position_extraction(df_dds)
        df_rld = position_extraction(df_rld)

        df_rld = dam_id_window_filtering(df_rld, chrom, center = center, window = window)

        df_dds = dam_id_window_filtering(df_dds, chrom, center = center, window = window)


        axs[1].bar(df_rld["start"] + (df_rld["length"]/2),
                    height = df_rld[condition] - df_rld["2RG3"],
                    width = df_rld["stop"] - df_rld["start"],
                    color = ["red" if log_ratio < 0 else "blue" for log_ratio in df_rld[condition] - df_rld["2RG3"]])

        axs[1].set_title(f"Gal3 only\nmean = {round(df_rld[condition].mean(), 2)}")
        axs[1].set_ylabel("rlog(ratio)")
        axs[1].set_xlabel("postions")

        i += 1


    plt.show()




def raw_comparison():
    
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

    conditions_list = ["3RG3D1", "3RG3D2", "3RG3D3", "2RG3"]

    fig = plt.figure(constrained_layout = True)
    fig.suptitle(chrom)

    subfigs = fig.subfigures(len(conditions_list),1)

    for condition in conditions_list:


        subfigs[i].suptitle(f"{condition}")
        
        axs = subfigs[i].subplots(1, 2,
                                sharex = True,
                                sharey = True)

        file_dds = "/datas/nathan/Dam_ID_analysis/dds_total.csv"



        df_dds = pd.read_csv(file_dds, header = 0, sep = ",")




        df_dds = position_extraction(df_dds)


        df_dds = dam_id_window_filtering(df_dds, chrom, center = center, window = window)


        axs[0].bar(df_dds["start"] + (df_dds["length"]/2),
                    height = df_dds[condition] ,
                    width = df_dds["stop"] - df_dds["start"],
                    color = ["red" if log_ratio < 0 else "blue" for log_ratio in df_dds[condition] - df_dds["2RG3"]])

        axs[0].set_title(f"all\nmean = {round(df_dds[condition].mean(), 2)}")
        axs[0].set_ylabel("rlog(ratio)")
        axs[0].set_xlabel("postions")

        
        
        file_dds = "/datas/nathan/Dam_ID_analysis/dds_G3.csv"


        df_dds = pd.read_csv(file_dds, header = 0, sep = ",")




        df_dds = position_extraction(df_dds)


        df_dds = dam_id_window_filtering(df_dds, chrom, center = center, window = window)


        axs[1].bar(df_dds["start"] + (df_dds["length"]/2),
                    height = df_dds[condition],
                    width = df_dds["stop"] - df_dds["start"],
                    color = ["red" if log_ratio < 0 else "blue" for log_ratio in df_dds[condition] - df_dds["2RG3"]])

        axs[1].set_title(f"Gal3 only\nmean = {(df_dds[condition].mean(), 2)}")
        axs[1].set_ylabel("rlog(ratio)")
        axs[1].set_xlabel("postions")

        i += 1


        print((df_dds.nlargest(5, condition)))
    
    plt.show()


raw_comparison()