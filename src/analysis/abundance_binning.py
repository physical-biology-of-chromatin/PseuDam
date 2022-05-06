import pandas as pd
import numpy as np

from packages.Dam_ID_utilities import kallisto_outup_reformating

from packages.binning import (abundance_file_conversion, binning_application_kallisto,
                              binning_kallisto_gatc,
                              binning_kallisto_fixed_size)

file_sites = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/Dam_ID/sites/sites.bed"

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

epsilon = 5

for gal_condition in range(1,4):
    
    conditions_list = [
        f"3RG{gal_condition}G1",
        f"3RG{gal_condition}G2",
        f"3RG{gal_condition}G3",
        f"3RG{gal_condition}",
        f"2RG{gal_condition}",
    ]

    for condition in conditions_list:
        
        
        file_abundance = f"/datas/nathan/Dam_ID_analysis/data/counts/{condition}/abundance.tsv"
        
        
        df_condition = kallisto_outup_reformating(file_abundance, file_sites, chrom_names)
        
        df_condition_binned = binning_kallisto_gatc(df_condition, epsilon, chrom = "all")
        
        path = f"/datas/nathan/Dam_ID_analysis/data/counts/{condition}/abundance_binned_gatc_{epsilon}.tsv"
        
        abundance_file_conversion(df_condition_binned, path)