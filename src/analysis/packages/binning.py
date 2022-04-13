import pandas as pd
import numpy as np


#####################################################################
#                                                                   #
#                           size binning                            #
#                                                                   #
#####################################################################


def binning_deseq2_fixed_size(df_in, 
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
                       chrom = "chrom_3",
                       region_start = 0,
                       region_end = 0):
    """Creates bins containing epsilon reads over the given region from a DESeq2 output pd df

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

    #TODO decide if the last bin should be included even if not full
    #     (included for now)
    
    
    if (region_end and region_start) == 0: 
        df_in = df_in[(df_in["chrom"] == chrom)]

    
    elif region_start >= region_end:
        raise ValueError("Region_end can't be equal or inferior to region start") 
    
    else:
        df_in = df_in[(df_in["chrom"] == chrom)
                      & (df_in["start"] > region_start)
                      & (df_in["stop"] < region_end)]        

    
    df_bins = pd.DataFrame()
    
    bin_start = 0
    bin_end = 0
    bin_value = 0
    
    bin_number = 0

    gatc_number = 0

    for index, row in df_in.iterrows():
        
        gatc_number += 1
        
        if bin_end == 0:
            bin_start = df_in["start"]
            gatc_number += 1
        
        bin_value += row[condition]
        bin_end = df_in["stop"]
        
        if bin_value >= epsilon:
            
            df_bins["value"][bin_number] = bin_value
            df_bins["start"][bin_number] = bin_start
            df_bins["stop"][bin_number] = bin_end
            df_bins["gatc"][bin_number] = gatc_number
            
            
            bin_end = 0
            bin_start = 0
            bin_value = 0            
            gatc_number = 0
            
            bin_number += 1 

        df_bins["value"][bin_number] = bin_value
        df_bins["start"][bin_number] = bin_start
        df_bins["stop"][bin_number] = bin_end
        df_bins["gatc"][bin_number] = gatc_number
    
                      
    return df_bins



def binning_kallisto_reads(df_in, 
                       epsilon, 
                       chrom = "chrom_3",
                       region_start = 0,
                       region_end = 0):
    """Creates bins containing epsilon reads over the given region from a kallisto abundance pd df

    Args:
        df_in (pd dataframe): DESeq2 output dataframe 
        (formated by kallisto_out_reformating)
        
        epsilon (int): number of reads in each bin
        
        chrom (str, optional): chromosome to bin on. Defaults to "chrom_3".
        
        region_start (int, optional): start of the binning region. Defaults to 0.
        
        region_end (int, optional): end of the binning region. Defaults to the length of the chromosome if 0.

    Returns:
        pd dataframe: binned dataset (contains only the binned region)
    """

    #TODO decide if the last bin should be included even if not full
    #     (included for now)
    
    if (region_end and region_start) == 0: 
        df_in = df_in[(df_in["chrom"] == chrom)]

    
    elif region_start >= region_end:
        raise ValueError("Region_end can't be equal or inferior to region start") 
    
    else:
        df_in = df_in[(df_in["chrom"] == chrom)
                      & (df_in["start"] > region_start)
                      & (df_in["stop"] < region_end)]        
    
    df_bins = pd.DataFrame()
    
    bin_start = 0
    bin_end = 0
    bin_value = 0
    
    bin_number = 0

    gatc_number = 0

    for index, row in df_in.iterrows():
        
        gatc_number += 1
        
        if bin_end == 0:
            bin_start = df_in["start"]
            gatc_number += 1
        
        bin_value += row["est_counts"]
        bin_end = df_in["stop"]
        
        if bin_value >= epsilon:
            
            df_bins["value"][bin_number] = bin_value
            df_bins["start"][bin_number] = bin_start
            df_bins["stop"][bin_number] = bin_end
            df_bins["gatc"][bin_number] = gatc_number
            
            
            bin_end = 0
            bin_start = 0
            bin_value = 0            
            gatc_number = 0
            
            bin_number += 1 

        df_bins["value"][bin_number] = bin_value
        df_bins["start"][bin_number] = bin_start
        df_bins["stop"][bin_number] = bin_end
        df_bins["gatc"][bin_number] = gatc_number
    
                     
    return df_bins



#####################################################################
#                                                                   #
#                           GATC binning                            #
#                                                                   #
#####################################################################


def binning_deseq_reads(df_in, 
                       condition, 
                       epsilon, 
                       chrom = "chrom_3",
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



    if (region_end and region_start) == 0: 
        df_in = df_in[(df_in["chrom"] == chrom)]

    
    elif region_start >= region_end:
        raise ValueError("Region_end can't be equal or inferior to region start") 
    
    else:
        df_in = df_in[(df_in["chrom"] == chrom)
                      & (df_in["start"] > region_start)
                      & (df_in["stop"] < region_end)]  

    df_bins = pd.DataFrame()
    
    bin_start = 0
    bin_end = 0
    bin_value = 0
    
    bin_number = 0

    gatc_number = 0


    for index, row in df_in.iterrows():
        
        gatc_number += 1
        
        if bin_end == 0:
            bin_start = row["start"]
            
        bin_end = row["stop"]
        bin_value += row[condition]
        
        if gatc_number >= epsilon:
            
            df_bins["start"][bin_number] = bin_end
            df_bins["stop"][bin_number] = bin_start
            df_bins["value"][bin_number] = bin_value
            
            bin_start = 0
            bin_end = 0
            bin_value = 0
            gatc_number = 0
            
            bin_number += 1
        
    df_bins["start"][bin_number] = bin_end
    df_bins["stop"][bin_number] = bin_start
    df_bins["value"][bin_number] = bin_value

    return df_bins



def binning_kallisto_gatc(df_in, 
                       epsilon, 
                       chrom = "chrom_3",
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
    
    
    
    if (region_end and region_start) == 0: 
        df_in = df_in[(df_in["chrom"] == chrom)]

    
    elif region_start >= region_end:
        raise ValueError("Region_end can't be equal or inferior to region start") 
    
    else:
        df_in = df_in[(df_in["chrom"] == chrom)
                      & (df_in["start"] > region_start)
                      & (df_in["stop"] < region_end)]  
        
        
    df_bins = pd.DataFrame()
    
    bin_start = 0
    bin_end = 0
    bin_value = 0
    
    bin_number = 0

    gatc_number = 0


    for index, row in df_in.iterrows():
        
        gatc_number += 1
        
        if bin_end == 0:
            bin_start = row["start"]
            gatc_number += 1
            
            
        bin_end = row["stop"]
        bin_value += row["est_counts"]
        
        if gatc_number >= epsilon:
            
            df_bins["start"][bin_number] = bin_end
            df_bins["stop"][bin_number] = bin_start
            df_bins["value"][bin_number] = bin_value
            
            bin_start = 0
            bin_end = 0
            bin_value = 0
            gatc_number = 0
            
            bin_number += 1
        
    df_bins["start"][bin_number] = bin_end
    df_bins["stop"][bin_number] = bin_start
    df_bins["value"][bin_number] = bin_value

    return df_bins