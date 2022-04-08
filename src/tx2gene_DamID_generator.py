import pandas as pd
sites_file = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites/test/yes/sites.bed"

pseudo_file = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/pseudo/test/yes/N_3RG1_EKDL220001504-1a_H2C2YDSX3_L3/abundance.tsv"

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

df_sites = pd.read_csv(sites_file, sep = "\t", header = None)

df_pseudo = pd.read_csv(pseudo_file, sep = "\t", header = 0)

df_pseudo.iloc[0,0] = "_0"

df_tx2gene = pd.DataFrame()

sample_name = [chrom_names[chrom_id] for chrom_id in list(df_sites[0])]

start = [str(start) for start in list(df_sites[1])]

stop = [str(stop) for stop in list(df_sites[2])]

gene_names = list(map(lambda sample_name, start, stop: sample_name + "_" + start + "_" + stop,
                      sample_name, 
                      start, 
                      stop))

df_tx2gene["TXNAME"] = df_pseudo["target_id"]

df_tx2gene["GENEID"] = gene_names

path = "/datas/nathan/Dam_ID_analysis/geneid_dam_id.csv"

df_tx2gene.to_csv(path, sep = "\t", index = False)