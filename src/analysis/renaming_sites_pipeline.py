import pandas as pd

chrom_names = {
        "ENA|BX284601|BX284601.5" : "I",
        "ENA|BX284602|BX284602.5" : "II",
        "ENA|BX284603|BX284603.4" : "III",
        "ENA|BX284604|BX284604.4" : "IV",
        "ENA|BX284605|BX284605.5" : "V",
        "ENA|BX284606|BX284606.5" : "X"
        }


path_sites = "/datas/nathan/genome_test_pipeline/Celegans.gff"

file_sites = pd.read_csv(path_sites, sep = "\t", header = None)

files_sites_mod = file_sites

files_sites_mod[0] = [chrom_names[name] for name in file_sites[0]]

files_sites_mod.to_csv("/datas/nathan/genome_test_pipeline/Celegans_mod.gff",
                       sep = "\t",
                       header = None,
                       index = None)