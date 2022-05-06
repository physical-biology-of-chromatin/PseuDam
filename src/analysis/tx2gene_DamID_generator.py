import pandas as pd
import re

def remove_substring(string, substring):
    
    new_string = string.removeprefix(substring)
    
    return new_string


sites_file = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/Dam_ID/c_elegans/sites/sites.bed"

chrom_names_yeast = {
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


chrom_names_celegans = {
    "ENA|BX284601|BX284601.5" : "chromosome_I",
    "ENA|BX284602|BX284602.5" : "chromosome_II",
    "ENA|BX284603|BX284603.4" : "chromosome_III",
    "ENA|BX284604|BX284604.4" : "chromosome_IV",
    "ENA|BX284605|BX284605.5" : "chromosome_V",
    "ENA|BX284606|BX284606.5" : "chromosome_X"
}



df_sites = pd.read_csv(sites_file, sep = "\t", header = None)

df_tx2gene = pd.DataFrame()

list_ids_line = [chrom_names_celegans[re.search("\A(.+?)_", id).group(1)] for id in df_sites[3]]

list_no_id = list(map(remove_substring, list(df_sites[3]), list(df_sites[0])))

gene_names_chrom = list(map(lambda chrom_name, no_id: chrom_name + no_id, list_ids_line, list_no_id))

gene_names_id = df_sites[3]

"""
sample_name = [chrom_names[chrom_id] for chrom_id in list(df_sites[0])]

start = [str(start) for start in list(df_sites[1])]

stop = [str(stop) for stop in list(df_sites[2])]

gene_names = list(map(lambda sample_name, start, stop: sample_name + "_" + start + "_" + stop,
                      sample_name, 
                      start, 
                      stop))
"""

df_tx2gene["TXNAME"] = gene_names_id

df_tx2gene["GENEID"] = gene_names_chrom

path = "/datas/nathan/Dam_ID_analysis/data/counts/C_elegans/dpy7/L2/tx2gene_ids.cvs"

df_tx2gene.to_csv(path, sep = "\t", index = False)