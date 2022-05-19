import pandas as pd
import numpy as np
import re


def position_extraction(df):
    """Extracts the fragments' position from the id of DESeq2 outputs

    Args:
        df (pd dataframe): your DESeq2 extracted dataframe

    Returns:
        pandas dataframe: your df with positions (chrom, start, stop)
    """


    try: 
        x = df["id"]
    except KeyError:
        df.rename(columns={"Unnamed: 0" : "id"}, inplace = True)
    
    # If the first letter is c takes everything until the second _ if it starts with an m takes until the first _
    df["chrom"] = [re.search("^([^_]*_[^_]*)_.*$", id).group(1) 
                    if id[0] == "c" 
                    else re.search("^(.+?)_", id).group(1) 
                    if id[0] == "m" 
                    else "wth" for id in df["id"]]


    # Reverses everything then gets eveything between the 2 first _
    df["start"] = [int((re.search("_(.+?)_", id[::-1]).group(1))[::-1]) for id in df["id"]]

    # Reverses the list then takes until the first _
    df["stop"] = [int((re.search("^(.+?)_", id[::-1]).group(1))[::-1]) for id in df["id"]]    

    df["length"] = df["stop"] - df["start"]

    return df


def kallisto_abundance_reformating(file_abundance):
    """Extracts the fragments' position from an abudance.tsv file from kallisto

    Args:
        df (pd dataframe): your DESeq2 extracted dataframe

    Returns:
        pandas dataframe: your df with positions (chrom, start, stop)
    """

    df = pd.read_csv(file_abundance,
                     sep = "\t",
                     header = 0)
    
    # If the first letter is c takes everything until the second _ if it starts with an m takes until the first _
    df["chrom"] = [re.search("^([^_]*_[^_]*)_.*$", id).group(1) 
                    if id[0] == "c" 
                    else re.search("^(.+?)_", id).group(1) 
                    if id[0] == "m" 
                    else "wth" for id in df["target_id"]]


    # Reverses everything then gets eveything between the 2 first _
    df["start"] = [int((re.search("_(.+?)_", id[::-1]).group(1))[::-1]) for id in df["target_id"]]

    # Reverses the list then takes until the first _
    df["stop"] = [int((re.search("^(.+?)_", id[::-1]).group(1))[::-1]) for id in df["target_id"]]    

    df["length"] = df["stop"] - df["start"]

    return df

def kallisto_abundance_extraction(file_abundance, chrom_dic):
    """Extracts the fragments' position from an abudance.tsv file from kallisto

    Args:
        df (pd dataframe): your DESeq2 extracted dataframe

    Returns:
        pandas dataframe: your df with positions (chrom, start, stop)
    """

    df = pd.read_csv(file_abundance,
                     sep = "\t",
                     header = 0,
                     index_col = False,
                     )
    
    # If the first letter is c takes everything until the second _ if it starts with an m takes until the first _
    df["chrom"] = [chrom_dic[re.search("^(.+?)_", id).group(1)] 
                   for id in df["target_id"]]


    # Reverses everything then gets eveything between the 2 first _
    df["start"] = [int((re.search("_(.+?)_", id[::-1]).group(1))[::-1]) for id in df["target_id"]]

    # Reverses the list then takes until the first _
    df["stop"] = [int((re.search("^(.+?)_", id[::-1]).group(1))[::-1]) for id in df["target_id"]]    

    df["length"] = df["stop"] - df["start"]

    return df


def sleuth_norm_extraction_cel(file_abundance, chrom_dic):
    """Extracts the fragments' position from an abudance.tsv file from kallisto

    Args:
        df (pd dataframe): your DESeq2 extracted dataframe

    Returns:
        pandas dataframe: your df with positions (chrom, start, stop)
    """

    df = pd.read_csv(file_abundance,
                     sep = ",",
                     header = 0,
                     index_col = False,
                     )
    
    # If the first letter is c takes everything until the second _ if it starts with an m takes until the first _
    df["chrom"] = [chrom_dic[re.search("^(.+?)_", id).group(1)] 
                   for id in df["target_id"]]


    # Reverses everything then gets eveything between the 2 first _
    df["start"] = [int((re.search("_(.+?)_", id[::-1]).group(1))[::-1]) for id in df["target_id"]]

    # Reverses the list then takes until the first _
    df["stop"] = [int((re.search("^(.+?)_", id[::-1]).group(1))[::-1]) for id in df["target_id"]]    

    df["length"] = df["stop"] - df["start"]

    return df



def sleuth_norm_extraction_sac(file_abundance, chrom_dic):
    """Extracts the fragments' position from an abudance.tsv file from kallisto

    Args:
        df (pd dataframe): your DESeq2 extracted dataframe

    Returns:
        pandas dataframe: your df with positions (chrom, start, stop)
    """

    df = pd.read_csv(file_abundance,
                     sep = ",",
                     header = 0,
                     index_col = False,
                     )
    
    # If the first letter is c takes everything until the second _ if it starts with an m takes until the first _
    df["chrom"] = [chrom_dic[re.search("^([^_]*_[^_]*)_.*$", id).group(1)] 
                   for id 
                   in df["target_id"]]



    # Reverses everything then gets eveything between the 2 first _
    start_stop = [(re.search("^([^_]*_[^_]*)_.*$", id[::-1]).group(1))[::-1] for id in df["target_id"]]

    # Reverses everything then gets eveything between the 2 first _
    df["start"] = [int(re.search("^(.+?)_", id).group(1)) for id in start_stop]

    # Reverses the list then takes until the first _
    df["stop"] = [int((re.search("^(.+?)_", id[::-1]).group(1))[::-1]) for id in start_stop]    

    df["length"] = df["stop"] - df["start"]

    return df



def sleuth_output_reformating(file_sleuth, chrom_dic = None):
    """Extracts the fragments' position from  sleuth analsyis output file 

    Args:
        df (pd dataframe): your DESeq2 extracted dataframe
        chrom_dic (dic): dictionnary to convert chrom id to chrom names
    Returns:
        pandas dataframe: your df with positions (chrom, start, stop)
    """

    df = pd.read_csv(file_sleuth,
                     sep = ",",
                     header = 0)
    
    if chrom_dic is None:
        # If the first letter is c takes everything until the second _ if it starts with an m takes until the first _
        df["chrom"] = [re.search("^([^_]*_[^_]*)_.*$", id).group(1) 
                    if id[0] == "c" 
                    else re.search("^(.+?)_", id).group(1) 
                    if id[0] == "m" 
                    else "wth" for id in df["target_id"]]

    else:
        df["chrom"] = [chrom_dic[re.search("^(.+?)_", id).group(1)]
                       for id in df["target_id"]]
        

    # Reverses everything then gets eveything between the 2 first _
    df["start"] = [int((re.search("_(.+?)_", id[::-1]).group(1))[::-1]) for id in df["target_id"]]

    # Reverses the list then takes until the first _
    df["stop"] = [int((re.search("^(.+?)_", id[::-1]).group(1))[::-1]) for id in df["target_id"]]    

    df["length"] = df["stop"] - df["start"]

    return df


def kallisto_outup_reformating(file_df, file_sites, chrom_names = {}):
    """Adds chrom, start and stop columns to the kallisto pipeline output

    Args:
        file_df (path): kallisto anbundance.tsv output path

        file_sites (path): path to the output .bed of GATC_finder
        
        chrom_names (dict, optional): dictionnary of chrom names to convert from IDs. Defaults to {}.

    Returns:
        pandas dataframe: abundance dataframe with fragment position
    """
    
    
    df_sites = pd.read_csv(file_sites,
                           sep = "\t",
                           header = None)

    df = pd.read_csv(file_df,
                     sep = "\t",
                     header = 0)


    df.insert(0,
              column = "chrom",
              value = df_sites[0])

    df.insert(1,
              column = "start",
              value = df_sites[1])

    df.insert(2,
              column = "stop",
              value = df_sites[2])

    if len(chrom_names) != 0:
        
        df["chrom"] = [chrom_names[chrom] for chrom in df["chrom"]]

    return df


def dam_id_window_filtering(df, chrom, window = 0, center = 0):
    """Filters a pandas df to contain only the desired region in the desired chromosome

    Args:
        df (pandas dataframe): GATC fragments
        chrom (str): desired chromosome to filter by
        window (int, optional): window around the center of the region. Defaults to 0.
        center (int, optional): center of the region. Defaults to 0.

    Raises:
        ValueError: Raised if you don't enter both a window and a center


    Returns:
        pandas dataframe: filtered GATC fragments
    """
    
    if window == 0 and center != 0:
        raise ValueError("Please enter both a window and a center")
    elif window != 0 and center == 0:
        raise ValueError("Please enter both a window and a center")
    
    df = df[df["chrom"] == chrom]
    
    if df["stop"].max() < (center + window/2):
        print("Warning: the upper limit of your window is exceeding the gatc fragments span")
    if center - window/2 < 0:
        print("Warning: your lower window is is negative")
        
    
    if window != 0:
        df = df[df["start"] > center - window/2]
        df = df[df["stop"] < center + window/2]

    return df
        


def chrom_rename(path, chrom_dic):
    """Renames the chrom names of the chrom column of a bedgraph

    Args:
        path (string): path to the bedgraph
        chrom_dic (dic): dic to rename the chroms

    Returns:
        pd dataframe: df with renamed chroms
    """
    
    
    df = pd.read_csv(path,
                     sep = "\t",
                     header = None,
                     skiprows=1)
    
    df.columns = ["chrom", "start", "stop", "value"]

    
    df["chrom"] = [chrom_dic[name]
                    for name in df["chrom"]]
    
    first_line = str(open(path).readline())
    
    df.to_csv(path,
              sep = "\t",
              index = False,
              header = None)
    
    with open(path, 'r+') as file:
    
        content = file.read()
        
        file.seek(0)
        
        file.write(first_line + content)
        
        
def gff_to_abundance(file_bed, path):
    
    df_bed = pd.read_csv(file_bed, sep = "\t", header = None)

    df_abundance = pd.DataFrame()
    
    df_abundance["target_id"] = [f"chr{chrom}"+"_"+str(int(start))+"_"+str(int(end)) 
                                 for chrom, start, end 
                                 in zip(df_bed[0], 
                                        df_bed[3], 
                                        df_bed[4])]



    df_abundance["length"] = [int(stop - start) for stop, start in zip(df_bed[4], df_bed[3])]
    
    df_abundance["eff_length"] = [length / 150 
                                  if (length / 150) >= 1
                                  else 1
                                  for length in df_abundance["length"]]      
    

    df_abundance["est_counts"] = df_bed[9]
    

    
    rpk_list = df_abundance["est_counts"] / 10000000
    
    df_abundance["tpm"] = [rpk / rpk_list.sum() for rpk in rpk_list]
    
    
    
    df_abundance.to_csv(path, sep = "\t", header = True, index = False)


def bed_sites_chrom_reformating(path_bed, path_out):
    
    df_bed = pd.read_csv(path_bed, header=None, sep="\t")
    
    
    chrom_names = {
        "ENA|BX284601|BX284601.5" : "I",
        "ENA|BX284602|BX284602.5" : "II",
        "ENA|BX284603|BX284603.4" : "III",
        "ENA|BX284604|BX284604.4" : "IV",
        "ENA|BX284605|BX284605.5" : "V",
        "ENA|BX284606|BX284606.5" : "X"
        }
    
    df_bed[0] = [chrom_names[name] for name in df_bed[0]]
    
    df_bed.to_csv(path_out, header = None, index = False, sep = "\t")
    
    
    
bed_sites_chrom_reformating("/datas/nathan/test_pipeline/test/sites.bed", "/datas/nathan/test_pipeline/test/sites_reform.bed")

    
#gff_to_abundance("/datas/nathan/test_pipeline/dpy7_L2_rep1.gff", "/datas/nathan/test_pipeline/abundance.tsv")