
#!/usr/bin/env python

"""
Description :

Gives the mean length of the fragments from a given bed file  or the std of the length

Input :
--bed : path to the .bed file containing the fragments
--mean : if given, prints the mean of the fragments length
--std : if given, prints the std of the fragments length distribution
    
Output :
mean : mean of the fragments' length
std : std of the fragments' length

"""


import argparse

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("--bed", action="store",
                    help = "<path> Path to the bed file containing the fragments")

parser.add_argument("--mean", action="store_true")

parser.add_argument("--std", action="store_true")



args = parser.parse_args()

bed_file = args.bed
mean_length = args.mean
std_length = args.std

df_bed = pd.read_csv(bed_file, header = None, sep = "\t")

if mean_length != False and std_length != False:
    raise ValueError("Both std and mean can't be specified at once")

if mean_length != False:
    
    mean = np.mean(df_bed[2] - df_bed[1])
    
    print(mean)
    
if std_length != False:
    
    std = np.std(df_bed[2] - df_bed[1])
    
    print(mean)
