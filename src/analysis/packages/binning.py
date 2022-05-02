import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



#####################################################################
#                                                                   #
#                           size binning                            #
#                                                                   #
#####################################################################


def binning_deseq_fixed_size(df_in, 
                              condition, 
                              bin_size,
                              chrom = "chrom_3", 
                              region_start = 0, 
                              region_end = 330000, 
                              percent = True):
    """Bins the given gatc fragments from a pd dataframe
    into bins of a given size 

    Args:
        df_in (pd dataframe): DESeq2 output pd df
        
        condition (_type_): Condition to pool
        
        bin_size (_type_): size of the bins
        
        chrom (str, optional): Chromosome to pool on. Defaults to "chrom_3".
        
        region_start (int, optional): Start of the pooled region. Defaults to 0.
        
        region_end (int, optional): End of teh pooled region. Defaults to 330000.
        
        percent (bool, optional): gives the bin either the percentage of the read
        count of the overlapping fragment or the whole count but divides the bin 
        by the number of fragments. Defaults to True.

    Returns:
        pd dataframe : binned dataframe
    """

    
    bins = np.arange(region_start,
                     region_end + (region_end % bin_size),
                     bin_size)
    

    df_in = df_in[(df_in["chrom"] == chrom) &
                  (df_in["start"] > region_start) &
                  (df_in["stop"] < region_end)]
    
    
    
    df_bins = pd.DataFrame(bins, columns = ["start"])
    
    df_bins["stop"] = df_bins["start"] + 500
    df_bins["start"] = df_bins["start"] + 1
    
    df_bins[condition] = np.zeros(len(df_bins))
    
    
    df_in.sort_values("start", inplace = True)
    df_in.reset_index(drop = True, inplace = True)
    
    df_in = df_in[["chrom", "start", "stop", "length", condition]]
    
    bin_number = 0
    value_bin = 0
    frag_number = 0
    
    for index, row in df_in.iterrows():
        
        start = row["start"]
        stop = row["stop"]
        length = row["length"]
        value = row[condition]
        
        
        if bin_number == 0:
            if not df_bins["start"][bin_number] > start:
                while not (df_bins["start"][bin_number] 
                           < start 
                           < df_bins["stop"][bin_number]):
                    
                    bin_number += 1
                
        if ((df_bins["start"][bin_number] < start)
            and df_bins["stop"][bin_number] > stop):
            
            value_bin += value
            
            if percent != True:
                
                frag_number += 1
            
            
        elif stop > df_bins["stop"][bin_number]:
            
            if percent == True:
                value_bin += (value 
                              * ((df_bins["stop"][bin_number] 
                                  - start) 
                                 / length))
                
            else:
                value_bin += value
                frag_number += 1
            
            while stop > df_bins["stop"][bin_number]:
                
                if percent == True:
                    df_bins[condition][bin_number] = value_bin
                else:
                    df_bins[condition][bin_number] = value_bin / frag_number
                    frag_number = 1
                    
                value_bin = 0
                
                bin_number += 1

                if percent == True:
                    
                    if ((df_bins["stop"][bin_number] - start) 
                        > length):
                        
                        value_bin += value * (bin_size / length)
                        
                    else:
                        
                        value_bin += (value 
                                      * ((df_bins["stop"][bin_number] 
                                          - start) 
                                         / length))
                        
                else:
                    
                    value_bin += value
                    frag_number += 1
        
        if percent == True:
            df_bins[condition][bin_number] = value_bin
        else:
            df_bins[condition][bin_number] = value_bin / frag_number
            frag_number = 1
        
    
    return df_bins



def binning_kallisto_fixed_size(df_in, 
                                bin_size,
                                chrom = "chrom_3", 
                                region_start = 0, 
                                region_end = 330000, 
                                percent = True):
    
    """Bins the given gatc fragments from a pd dataframe
    into bins of a given size 

    Args:
        df_in (pd dataframe): Kallisto output pd df
        
        condition (_type_): Condition to pool
        
        bin_size (_type_): size of the bins
        
        chrom (str, optional): Chromosome to pool on. Defaults to "chrom_3".
        
        region_start (int, optional): Start of the pooled region. Defaults to 0.
        
        region_end (int, optional): End of teh pooled region. Defaults to 330000.
        
        percent (bool, optional): gives the bin either the percentage of the read
        count of the overlapping fragment or the whole count but divides the bin 
        by the number of fragments. Defaults to True.

    Returns:
        pd dataframe : binned dataframe
    """


    bins = np.arange(region_start,
                     region_end + (region_end % bin_size),
                     bin_size)
    

    df_in = df_in[(df_in["chrom"] == chrom) &
                  (df_in["start"] > region_start) &
                  (df_in["stop"] < region_end)]
    
    
    
    df_bins = pd.DataFrame(bins, columns = ["start"])
    
    df_bins["stop"] = df_bins["start"] + 500
    df_bins["start"] = df_bins["start"] + 1
    
    df_bins["est_counts"] = np.zeros(len(df_bins))
    
    
    df_in.sort_values("start", inplace = True)
    df_in.reset_index(drop = True, inplace = True)
    
    df_in = df_in[["chrom", "start", "stop", "length", "est_counts"]]
    
    bin_number = 0
    value_bin = 0
    frag_number = 0
    
    for index, row in df_in.iterrows():
        
        start = row["start"]
        stop = row["stop"]
        length = row["length"]
        value = row["est_counts"]
        
        
        if bin_number == 0:
            if not df_bins["start"][bin_number] > start:
                while not (df_bins["start"][bin_number] 
                           < start 
                           < df_bins["stop"][bin_number]):
                    
                    bin_number += 1
                
        if ((df_bins["start"][bin_number] < start)
            and df_bins["stop"][bin_number] > stop):
            
            value_bin += value
            
            if percent != True:
                
                frag_number += 1
            
            
        elif stop > df_bins["stop"][bin_number]:
            
            if percent == True:
                value_bin += (value 
                              * ((df_bins["stop"][bin_number] 
                                  - start) 
                                 / length))
                
            else:
                value_bin += value
                frag_number += 1
            
            while stop > df_bins["stop"][bin_number]:
                
                if percent == True:
                    df_bins["est_counts"][bin_number] = value_bin
                else:
                    df_bins["est_counts"][bin_number] = value_bin / frag_number
                    frag_number = 1
                    
                value_bin = 0
                
                bin_number += 1

                if percent == True:
                    
                    if ((df_bins["stop"][bin_number] - start) 
                        > length):
                        
                        value_bin += value * (bin_size / length)
                        
                    else:
                        
                        value_bin += (value 
                                      * ((df_bins["stop"][bin_number] 
                                          - start) 
                                         / length))
                        
                else:
                    
                    value_bin += value
                    frag_number += 1
        
        if percent == True:
            df_bins["est_counts"][bin_number] = value_bin
            
        else:
            df_bins["est_counts"][bin_number] = value_bin / frag_number
            frag_number = 1
        
    
    return df_bins


#####################################################################
#                                                                   #
#                           reads binning                           #
#                                                                   #
#####################################################################


def binning_deseq_reads(df_in, 
                       condition, 
                       epsilon, 
                       chrom = "all",
                       region_start = 0,
                       region_end = 0):
    """Creates bins containing epsilon reads over the given region from a DESeq2 output pd df

    Args:
        df_in (pd dataframe): DESeq2 output dataframe 
        (formated by position_extraction) 
        
        condition (stirng): name of the column containing the counts
        
        epsilon (int): number of reads in each bin
        
        chrom (str, optional): chromosome to bin on. If = all bins all the chromsomomes in 
        the dataset. Defaults to "chrom_3".
        
        region_start (int, optional): start of the binning region. Defaults to 0.
        
        region_end (int, optional): end of the binning region. 
        Defaults to the length of the chromosome if 0.

    Returns:
        pd dataframe: binned dataset (contains only the binned region)
    """

    #TODO decide if the last bin should be included even if not full
    #     (included for now)
    
    if type(df_in) != pd.DataFrame:
        raise TypeError("df_in must be a pandas dataframe")
        
        
    chrom_list = list()
    
    if chrom == "all": 
        chrom_list = df_in["chrom"].unique()
        
        
        
    elif (region_end and region_start) == 0: 
        df_in = df_in[(df_in["chrom"] == chrom)]
        chrom_list.append(chrom)
        
        
    elif region_start >= region_end:
        raise ValueError("Region_end can't be equal or inferior to region start") 
    
    else:
        df_in = df_in[(df_in["chrom"] == chrom)
                      & (df_in["start"] > region_start)
                      & (df_in["stop"] < region_end)]        
        chrom_list.append(chrom)
    
    df_bins = pd.DataFrame()
    df_bins["start"] = np.zeros(len(df_in))
    df_bins["stop"] = np.zeros(len(df_in))
    df_bins["value"] = np.zeros(len(df_in))
    df_bins["chrom"] = np.zeros(len(df_in))
    df_bins["gatc"] = np.zeros(len(df_in))
    
    
    bin_start = 0
    bin_end = 0
    bin_value = 0
    
    bin_number = 0

    gatc_number = 0

    

    for chromosome in chrom_list:
        
        df_chrom = df_in
        
        bin_start = 0
        bin_end = 0
        bin_value = 0
        gatc_number = 0
        
        if chrom == "all":
            df_chrom = df_in[(df_in["chrom"] == chromosome)]
            
        df_chrom.sort_values("start", inplace = True)
        
        for index, row in df_chrom.iterrows():
            
            gatc_number += 1
            
            if bin_end == 0:
                bin_start = row["start"]
                gatc_number += 1
            
            bin_value += row[condition]
            bin_end = row["stop"]
            
            if bin_value >= epsilon:
                
                df_bins["value"][bin_number] = bin_value
                df_bins["start"][bin_number] = bin_start
                df_bins["stop"][bin_number] = bin_end
                df_bins["gatc"][bin_number] = gatc_number
                df_bins["chrom"][bin_number] = chromosome
                
                bin_end = 0
                bin_start = 0
                bin_value = 0            
                gatc_number = 0
                
                bin_number += 1 

        df_bins["value"][bin_number] = bin_value
        df_bins["start"][bin_number] = bin_start
        df_bins["stop"][bin_number] = bin_end
        df_bins["gatc"][bin_number] = gatc_number
        bin_number += 1
        
        
    df_bins = df_bins[df_bins["stop"] != 0]
        
    return df_bins



def binning_kallisto_reads(df_in, 
                       epsilon, 
                       chrom = "all",
                       region_start = 0,
                       region_end = 0):
    """Creates bins containing epsilon reads over the given region from a kallisto abundance pd df

    Args:
        df_in (pd dataframe): DESeq2 output dataframe 
        (formated by kallisto_out_reformating)
        
        epsilon (int): number of reads in each bin
        
        chrom (str, optional): chromosome to bin on. Defaults to "all".
        
        region_start (int, optional): start of the binning region. Defaults to 0.
        
        region_end (int, optional): end of the binning region. Defaults to the length of the chromosome if 0.

    Returns:
        pd dataframe: binned dataset (contains only the binned region)
    """

    #TODO decide if the last bin should be included even if not full
    #     (included for now)
    
    if type(df_in) != pd.DataFrame:
        raise TypeError("df_in must be a pandas dataframe")
        
        
    chrom_list = list()
    
    if chrom == "all": 
        chrom_list = df_in["chrom"].unique()
        
        
        
    elif (region_end and region_start) == 0: 
        df_in = df_in[(df_in["chrom"] == chrom)]
        chrom_list.append(chrom)
        
        
    elif region_start >= region_end:
        raise ValueError("Region_end can't be equal or inferior to region start") 
    
    else:
        df_in = df_in[(df_in["chrom"] == chrom)
                      & (df_in["start"] > region_start)
                      & (df_in["stop"] < region_end)]        
        chrom_list.append(chrom)
    
    df_bins = pd.DataFrame()
    df_bins["start"] = np.zeros(len(df_in))
    df_bins["stop"] = np.zeros(len(df_in))
    df_bins["value"] = np.zeros(len(df_in))
    df_bins["chrom"] = np.zeros(len(df_in))
    df_bins["gatc"] = np.zeros(len(df_in))
    
    
    bin_start = 0
    bin_end = 0
    bin_value = 0
    
    bin_number = 0

    gatc_number = 0

    

    for chromosome in chrom_list:
        
        df_chrom = df_in
        
        bin_start = 0
        bin_end = 0
        bin_value = 0
        gatc_number = 0
        
        if chrom == "all":
            df_chrom = df_in[(df_in["chrom"] == chromosome)]
            
        df_chrom.sort_values("start", inplace = True)
        
        for index, row in df_chrom.iterrows():
            
            gatc_number += 1
            
            if bin_end == 0:
                bin_start = row["start"]
                gatc_number += 1
            
            bin_value += row["est_counts"]
            bin_end = row["stop"]
            
            if bin_value >= epsilon:
                
                df_bins["value"][bin_number] = bin_value
                df_bins["start"][bin_number] = bin_start
                df_bins["stop"][bin_number] = bin_end
                df_bins["gatc"][bin_number] = gatc_number
                df_bins["chrom"][bin_number] = chromosome
                
                bin_end = 0
                bin_start = 0
                bin_value = 0            
                gatc_number = 0
                
                bin_number += 1 

        df_bins["value"][bin_number] = bin_value
        df_bins["start"][bin_number] = bin_start
        df_bins["stop"][bin_number] = bin_end
        df_bins["gatc"][bin_number] = gatc_number
        df_bins["chrom"][bin_number] = chromosome
        bin_number += 1
        
        
    df_bins = df_bins[df_bins["stop"] != 0]
        
    return df_bins


#####################################################################
#                                                                   #
#                           GATC binning                            #
#                                                                   #
#####################################################################


def binning_deseq_gatc(df_in, 
                       condition, 
                       epsilon, 
                       chrom = "all",
                       region_start = 0,
                       region_end = 330000):
    """Creates bins containing epsilon gatc_sites (epsilon+1 gatc fragments) 
    over the given region from a DESeq2 output pd df

    Args:
        df_in (pd dataframe): DESeq2 output dataframe 
        (formated by position_extraction) 
        
        condition (stirng): name of the column containing the counts
        
        epsilon (int): number of reads in each bin
        
        chrom (str, optional): chromosome to bin on. Defaults to "chrom_3".
        
        region_start (int, optional): start of the binning region. Defaults to 0.
        
        region_end (int, optional): end of the binning region. 
        Defaults to the length of the chromosome if 0.

    Returns:
        pd dataframe: binned dataset (contains only the binned region)
    """



    if type(df_in) != pd.DataFrame:
        raise TypeError("df_in must be a pandas dataframe")
        
        
    chrom_list = list()
    
    if chrom == "all": 
        chrom_list = df_in["chrom"].unique()
        
    
    elif (region_end and region_start) == 0: 
        df_in = df_in[(df_in["chrom"] == chrom)]
        chrom_list.append(chrom)
        
        
    elif region_start >= region_end:
        raise ValueError("Region_end can't be equal or inferior to region start") 
    
    else:
        df_in = df_in[(df_in["chrom"] == chrom)
                      & (df_in["start"] > region_start)
                      & (df_in["stop"] < region_end)]        
        chrom_list.append(chrom)
    
    df_bins = pd.DataFrame()
    df_bins["start"] = np.zeros(len(df_in))
    df_bins["stop"] = np.zeros(len(df_in))
    df_bins["value"] = np.zeros(len(df_in))
    df_bins["chrom"] = np.zeros(len(df_in))
    
    bin_start = 0
    bin_end = 0
    bin_value = 0
    
    bin_number = 0

    gatc_number = 0

    for chromosome in chrom_list:
        df_chrom = df_in
        
        bin_start = 0
        bin_end = 0
        bin_value = 0
        gatc_number = 0

        if chrom == "all":
            df_chrom = df_in[(df_in["chrom"] == chromosome)]
            
        df_chrom.sort_values("start", inplace = True)

        for index, row in df_chrom.iterrows():
            
            gatc_number += 1
            
            if bin_end == 0:
                bin_start = row["start"]
                
            bin_end = row["stop"]
            bin_value += row[condition]
            
            if gatc_number >= epsilon:
                
                df_bins["start"][bin_number] = bin_start
                df_bins["stop"][bin_number] = bin_end
                df_bins["value"][bin_number] = bin_value
                df_bins["chrom"][bin_number] = chromosome
                
                bin_start = 0
                bin_end = 0
                bin_value = 0
                gatc_number = 0
                
                bin_number += 1
            
        df_bins["start"][bin_number] = bin_start
        df_bins["stop"][bin_number] = bin_end
        df_bins["value"][bin_number] = bin_value
        df_bins["chrom"][bin_number] = chromosome
        
        bin_number += 1

    df_bins = df_bins[df_bins["stop"] != 0]
    
    
    return df_bins



def binning_kallisto_gatc(df_in, 
                       epsilon, 
                       chrom = "all",
                       region_start = 0,
                       region_end = 0):
    """Creates bins containing epsilon gatc_sites (epsilon-1 gatc fragments) over the 
    given region from a kallisto abundance pd df

    Args:
        df_in (pd dataframe): DESeq2 output dataframe 
        (formated by kallisto_out_reformating)
        
        epsilon (int): number of reads in each bin
        
        chrom (str, optional): chromosome to bin on. Defaults to "chrom_3".
        
        region_start (int, optional): start of the binning region. Defaults to 0.
        
        region_end (int, optional): end of the binning region. Defaults to the length 
        of the chromosome if 0.

    Returns:
        pd dataframe: binned dataset (contains only the binned region)
    """
    
    
    
    if type(df_in) != pd.DataFrame:
        raise TypeError("df_in must be a pandas dataframe")
        
        
    chrom_list = list()
    
    if chrom == "all": 
        chrom_list = df_in["chrom"].unique()
        
    
    elif (region_end and region_start) == 0: 
        df_in = df_in[(df_in["chrom"] == chrom)]
        chrom_list.append(chrom)
        
        
    elif region_start >= region_end:
        raise ValueError("Region_end can't be equal or inferior to region start") 
    
    else:
        df_in = df_in[(df_in["chrom"] == chrom)
                      & (df_in["start"] > region_start)
                      & (df_in["stop"] < region_end)]        
        chrom_list.append(chrom)
    
    df_bins = pd.DataFrame()
    df_bins["start"] = np.zeros(len(df_in))
    df_bins["stop"] = np.zeros(len(df_in))
    df_bins["value"] = np.zeros(len(df_in))
    df_bins["chrom"] = np.zeros(len(df_in))
    
    bin_start = 0
    bin_end = 0
    bin_value = 0
    
    bin_number = 0

    gatc_number = 0

    for chromosome in chrom_list:
        df_chrom = df_in
        
        bin_start = 0
        bin_end = 0
        bin_value = 0
        gatc_number = 0

        if chrom == "all":
            df_chrom = df_in[(df_in["chrom"] == chromosome)]
            
        df_chrom.sort_values("start", inplace = True)

        for index, row in df_chrom.iterrows():
            
            gatc_number += 1
            
            if bin_end == 0:
                bin_start = row["start"]
                
            bin_end = row["stop"]
            bin_value += row["est_counts"]
            
            if gatc_number >= epsilon:
                
                df_bins["start"][bin_number] = bin_start
                df_bins["stop"][bin_number] = bin_end
                df_bins["value"][bin_number] = bin_value
                df_bins["chrom"][bin_number] = chromosome
                
                bin_start = 0
                bin_end = 0
                bin_value = 0
                gatc_number = 0
                
                bin_number += 1
            
        df_bins["start"][bin_number] = bin_start
        df_bins["stop"][bin_number] = bin_end
        df_bins["value"][bin_number] = bin_value
        df_bins["chrom"][bin_number] = chromosome
        
        bin_number += 1

    df_bins = df_bins[df_bins["stop"] != 0]
    
    
    return df_bins



#####################################################################
#                                                                   #
#                       Binning Serpentine                          #
#                                                                   #
#####################################################################



def binning_serpentine_kalisto(df_test_mean,
                               df_control_mean,
                               epsilon,
                               teta):
    """Bins the given test and control datasets from kallisto abundance output
    using a 1D version of the serpentine algortithm from L. Baudry et al. 
    genome analysis 2020

    Args:
        df_test (list of dataframes): list of df containing the test condition datas
        df_control (list of dataframes): list of df containing the control condition datas
        epsilon (int): lower threshold value
        teta (int): higher threshold value

    Returns:
        pd dataframe: df containing the binned control condition
        pd dataframe: df containing the binned test condition
    """

    df_bins_test = pd.DataFrame()
    df_bins_control = pd.DataFrame()


    chrom_list = list()
    



    chrom_list = df_control_mean["chrom"].unique()
        



    df_bins_test["chrom"] = np.zeros(len(df_test_mean), dtype = str)
    df_bins_test["est_counts"] = np.zeros(len(df_test_mean))
    df_bins_test["start"] = np.zeros(len(df_test_mean))
    df_bins_test["stop"] = np.zeros(len(df_test_mean))
    df_bins_test["length"] = np.zeros(len(df_test_mean))


    df_bins_control["chrom"] = np.zeros(len(df_test_mean), dtype = str)
    df_bins_control["est_counts"] = np.zeros(len(df_test_mean))
    df_bins_control["start"] = np.zeros(len(df_test_mean))
    df_bins_control["stop"] = np.zeros(len(df_test_mean))
    df_bins_control["length"] = np.zeros(len(df_test_mean))

        
    

    i = 0
    
    for chrom in chrom_list:
        
        value_test = 0
        value_control = 0
    
        stop_pos = 0
        
        
        df_control_chrom_mean = df_control_mean[df_control_mean["chrom"] == chrom]
        df_control_chrom_mean = df_control_chrom_mean.sort_values(["start"])
        
        
        df_test_chrom_mean = df_test_mean[df_test_mean["chrom"] == chrom]
        df_test_chrom_mean = df_test_chrom_mean.sort_values(["start"])
        
        for ((index_test,
              row_test),
             (index_control,
              row_control)) in zip(df_test_chrom_mean.iterrows(),
                                   df_control_chrom_mean.iterrows()):

            value_test += row_test["est_counts"]
            value_control += row_control["est_counts"]
            #chrom_pos = row_test["chrom"]
            
            
            if stop_pos == 0:
                start_pos = row_test["start"]
            
            stop_pos = row_test["stop"]

            if ((value_control > teta 
                or value_test > teta) 
                and (value_control > epsilon 
                    and value_test > epsilon)):
                
                df_bins_test.at[i, "chrom"] = chrom
                df_bins_test.at[i, "est_counts"] = value_test
                df_bins_test.at[i, "start"] = start_pos
                df_bins_test.at[i, "stop"] = stop_pos
                df_bins_test.at[i, "length"] = stop_pos - start_pos
                
                
                df_bins_control.at[i, "chrom"] = chrom
                df_bins_control.at[i, "est_counts"] = value_control
                df_bins_control.at[i, "start"] = start_pos
                df_bins_control.at[i, "stop"] = stop_pos
                df_bins_control.at[i, "length"] = stop_pos - start_pos
        
                value_test = 0
                value_control = 0
                start_pos = 0
                stop_pos = 0
                
                i += 1

        df_bins_control = df_bins_control[df_bins_control["stop"] != 0]
        df_bins_test = df_bins_test[df_bins_test["stop"] != 0]

        df_bins_control["start"] = df_bins_control["start"].astype(int)
        df_bins_test["start"] = df_bins_test["start"].astype(int)
        
        df_bins_control["stop"] = df_bins_control["stop"].astype(int)
        df_bins_test["stop"] = df_bins_test["stop"].astype(int)
        
        

    return (df_bins_control, df_bins_test)


def mean_replicates(replicates_list):
    
    df_list_mean = pd.DataFrame()
    
    for replicate in replicates_list:
        
        try: 
            test = df_list_mean["est_counts"]
        except KeyError:
            df_list_mean["est_counts"] = np.zeros(len(replicate))
            df_list_mean["chrom"] = replicate["chrom"]
            df_list_mean["start"] = replicate["start"]
            df_list_mean["stop"] = replicate["stop"]
            df_list_mean["length"] = replicate["length"]
            df_list_mean["target_id"] = replicate["target_id"]

        df_list_mean["est_counts"] += replicate["est_counts"]

        df_list_mean["est_counts"] = df_list_mean["est_counts"]/len(replicates_list)

    
    return df_list_mean
    
    
def binning_serpentine_deseq(df_test,
                               df_control,
                               epsilon,
                               teta,
                               condition):
    """Bins the given test and control datasets from DESeq2 outup
    using a 1D version of the serpentine algortithm from L. Baudry et al. 
    genome analysis 2020

    Args:
        df_test (pd dataframe): df containing the test condition datas
        df_control (pd dataframe): df containing the control condition datas
        epsilon (int): lower threshold value
        teta (int): higher threshold value
        condition (string) : condition in the DESeq2 file to bin from


    Raises:
        TypeError: _description_
        ValueError: _description_

    Returns:
        pd dataframe: df containing the binned control condition
        pd dataframe: df containing the binned test condition
    """

    df_bins_test = pd.DataFrame()
    df_bins_control = pd.DataFrame()

    if type(df_control) != pd.DataFrame or type(df_test) != pd.DataFrame:
        raise TypeError("df_control and df_test must be pandas dataframes")
    
    

    chrom_list = list()
    
    
    chrom_list = df_control["chrom"].unique()
        

    df_bins_test[["value", "start", "stop", "length"]] = np.zeros(len(df_test))
    df_bins_control["value", "start", "stop", "length"] = np.zeros(len(df_test))

        
    
    value_test = 0
    value_control = 0
    
    stop_pos = 0
    
    i = 0
    
    
    for chrom in chrom_list:
        
        
        for ((index_test,
            row_test),
            (index_control,
            row_control)) in zip(df_test.iterrows(),
                                df_control.iterrows()):
        
            value_test += row_test[condition]
            value_control += row_control[condition]
            
            if stop_pos == 0:
                start_pos = row_test["start"]
            
            stop_pos = row_test["stop"]

            if ((value_control > teta 
                or value_test > teta) 
                and (value_control > epsilon 
                    and value_test > epsilon)):
                
                df_bins_test["value"][i] = value_test
                df_bins_test["start"][i] = start_pos
                df_bins_test["stop"][i] = stop_pos
                df_bins_test["length"][i] = stop_pos - start_pos
                
                
                df_bins_control["value"][i] = value_control
                df_bins_control["start"][i] = start_pos
                df_bins_control["stop"][i] = stop_pos
                df_bins_control["length"][i] = stop_pos - start_pos
        
                value_test = 0
                value_control = 0
                start_pos = 0
                stop_pos = 0
                
                i += 1

        df_bins_control = df_bins_control[df_bins_control["stop"] != 0]
        df_bins_test = df_bins_test[df_bins_test["stop"] != 0]
    
    return (df_bins_control, df_bins_test)


#####################################################################
#                                                                   #
#                             Utilities                             #
#                                                                   #
#####################################################################


def binning_application_deseq(df_in, df_control, condition):
    """ Applies a binning to other conditions

    Args:
        df_in (pd dataframe): Dataframe to bin
        
        df_control (pd dataframe): Dataframe to get the bins from

        condition (string): Name of the column containing the values
        
    Returns:
        pd dataframe: Binned dataframe
    """
    
    df_bins = pd.DataFrame(df_control[["start", "stop", "chrom"]])
    
    df_bins["value"] = np.zeros([len(df_bins)])
    
    
    for index, row in df_bins.iterrows():
        
        
        df_temp = df_in[(df_in["start"] >= row["start"])
                        & (df_in["stop"] <= row["stop"])
                        & (df_in["chrom"] == row["chrom"])]
        
        df_bins["value"][index] = df_temp[condition].sum()
        

    return df_bins



def tx2gene_geneid_generator_df(df_in, path_out):
    """Creates a geneid dataframe for DESeq2 from an existing dataframe
    and saves it if a path if given

    Args:
        df_in (pd dataframe): _description_
        path_out (string, optional): Path to save the geneIDs to. Defaults to None.

    Returns:
        pd dataframe: dataframe containing the geneIDs and the txnames
    """
    
    
    df_tx2gene = pd.DataFrame()

    sample_name = [chrom_name for chrom_name in list(df_in["chrom"])]

    start = [str(start) for start in list(df_in["start"])]

    stop = [str(stop) for stop in list(df_in["stop"])]

    gene_names = list(map(lambda sample_name, start, stop: sample_name + "_" + start + "_" + stop,
                        sample_name, 
                        start, 
                        stop))

    df_tx2gene["TXNAME"] = [("_" + str(i)) for i in range(len(gene_names))]

    df_tx2gene["GENEID"] = gene_names

    if path_out == None:
        pass
    elif type(path_out) == str:
        df_tx2gene.to_csv(path_out, sep = "\t", index = False)
    
    return df_tx2gene



def abundance_file_conversion(df_in, path, path_tx2gene = None):
    """Transforms a binned dataframe into an abundance file
    usable as an input for DESeq2

    Args:
        df_in (pd dataframe): binned dataframe
        
        path (string): path to the location to save the file
        
        path_tx2gene(string): path to the location to save the tx2gene file. Defaults to None
    """
    
    df_abundance = pd.DataFrame()

    df_abundance["target_id"] = [chrom+"_"+str(int(start))+"_"+str(int(end)) 
                                 for chrom, start, end 
                                 in zip(df_in["chrom"], 
                                        df_in["start"], 
                                        df_in["stop"])]

    df_abundance["length"] = [int(stop - start) for stop, start in zip(df_in["stop"], df_in["start"])]
    
    df_abundance["eff_length"] = np.zeros(len(df_in))
    
    df_abundance["est_counts"] = df_in["value"]
    
    df_abundance["tpm"] = np.zeros(len(df_in))
    
    
    

    
    df_abundance.to_csv(path, sep = "\t", header = True, index = False)


    if path_tx2gene != None:
        
        df_tx2gene = pd.DataFrame()
        
        df_tx2gene["TXNAME"] = df_abundance["target_id"]
        df_tx2gene["GENEID"] = df_abundance["target_id"]

        df_tx2gene.to_csv(path_tx2gene, sep = "\t", header = True, index = False)
        


def bed_fragments_generation(df_in, path, path_tx2gene = None):
    """Generates a .bed file using a binned dataframe

    Args:
        df_in (pd dataframe): Binned dataframe to create the bed file from
        path (string): Path to save the bed file at
        path_tx2gene (string, optional): Path to save the tx2gene conversion file at . Defaults to None.
    """
    
    
    df_bed = pd.DataFrame()
    
    chrom_names_celegans = {
    "chromosome_I" : "ENA|BX284601|BX284601.5",
    "chromosome_II" : "ENA|BX284602|BX284602.5",
    "chromosome_III" : "ENA|BX284603|BX284603.4",
    "chromosome_IV" : "ENA|BX284604|BX284604.4",
    "chromosome_V" : "ENA|BX284605|BX284605.5",
    "chromosome_X" : "ENA|BX284606|BX284606.5"
    }
    
    
    df_bed["chrom"] = [chrom_names_celegans[name] 
                       for name 
                       in df_in["chrom"]]
    df_bed["start"] = df_in["start"]
    df_bed["stop"] = df_in["stop"]
    df_bed["name"] = [chrom+"_"+str(int(start))+"_"+str(int(end)) 
                      for chrom, start, end 
                      in zip(df_in["chrom"], 
                             df_in["start"], 
                             df_in["stop"])]
    
    df_bed.to_csv(path, header = False, sep = "\t", index = False)
    
    
    if path_tx2gene is not None:
        
        df_tx2gene = pd.DataFrame()
        
        df_tx2gene["TXNAME"] = df_bed["name"]
        df_tx2gene["GENEID"] = df_bed["name"]

        df_tx2gene.to_csv(path_tx2gene, sep = "\t", header = True, index = False)


        
def binning_application_kallisto(df_in, df_control):
    """ Applies a binning to other conditions

    Args:
        df_in (pd dataframe): Dataframe to bin
        
        df_control (pd dataframe): Dataframe to get the bins from

    Returns:
        pd dataframe: Binned dataframe
    """
    
    df_bins = pd.DataFrame(df_control[["start", "stop", "chrom"]])
    
    df_bins["value"] = np.zeros([len(df_bins)])
    
    
    for index, row in df_bins.iterrows():
        
        
        df_temp = df_in[(df_in["start"] >= row["start"])
                        & (df_in["stop"] <= row["stop"])
                        & (df_in["chrom"] == row["chrom"])]
        
        df_bins["value"][index] = df_temp["est_counts"].sum()
        
    

    return df_bins



def sleuth_zeros_addition(df_in, bed_file):

    chrom_names_celegans = {
    "ENA|BX284601|BX284601.5" : "chromosome_I",
    "ENA|BX284602|BX284602.5" : "chromosome_II",
    "ENA|BX284603|BX284603.4" : "chromosome_III",
    "ENA|BX284604|BX284604.4" : "chromosome_IV",
    "ENA|BX284605|BX284605.5" : "chromosome_V",
    "ENA|BX284606|BX284606.5" : "chromosome_X"
    }



    bed = pd.read_csv(bed_file, header = None, sep = "\t")
    
    bed = bed.rename(columns = {0 : "chrom",
                                1 : "start",
                                2 : "stop",
                                3 : "target_id"})
    
    bed["chrom"] = [chrom_names_celegans[name] 
                    for name 
                    in bed["chrom"]]
    bed["length"] = [stop - start for start, stop in zip(bed["start"], bed["stop"])]
    bed["est_counts"] = np.zeros(len(bed))

    # Makes a list of the targets to replace by in the null df
    ref = df_in["target_id"].to_numpy()

    # Removes all the ids that are in the list
    df_out = bed[~bed["target_id"].isin(ref)]
    
    # Merge to replace the removed rows
    df_out = pd.merge(df_out, df_in, how = "outer")



    """
    i = 0
    ref = df_in["target_id"].to_numpy()
    
    
    for row in df_in.itertuples(index = False):
        
        i += 1
        
        index = bed[bed["target_id"] == row[5]].index.values[0]
        
        row_dic = {"chrom" : row[1],
                   "start" : row[2],
                   "stop" : row[3],
                   "target_id" : row[5],
                   "length" : row[4],
                   "est_counts" : row[0]}
        
        bed.loc[index] = row_dic
    
        if i%1000 == 0:
            print(i)
    """  
              
    return df_out


#####################################################################
#                                                                   #
#                          Binning Analysis                         #
#                                                                   #
#####################################################################


def binning_stats(df_in, show = True):
    
    try:
        test = df_in["length"]
    except KeyError:
        print("Warning : your df doesn't have a length column")
        df_in["length"] = df_in["stop"] - df_in["start"]
    
    fig, ax = plt.subplots()    
    
    ax.hist(df_in["length"],
             bins = 200)
    ax.set_title(f"mean = {df_in['length'].mean()}\n"
                 f"median = {df_in['length'].median()}\n"
                 f"std = {df_in['length'].std()}")
    ax.set_ylabel("Number of fragments")
    ax.set_xlabel("Fragment length")
    
    if show is True:
        plt.show()
    
    return fig
    
    