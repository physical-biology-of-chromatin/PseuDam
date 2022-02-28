import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/coverage.bed", 
                 header = None, 
                 sep = '\t')

sites_list = list()
values_list = list()

for row in df.itertuples():
    site = list()
    
    if row[1] == "NC_004354.4" and row[2] > 1000000 and row[3] < 15000000:
        
        site = row[2]

        
        value = (row[7])
        
        sites_list.append(site)
        values_list.append(value)


values = list()

for site, value in zip(sites_list, values_list):
    
    for i in range(round(value * 100)):
        values.append(site + 3)

plt.hist(values, bins=sites_list)
plt.show()