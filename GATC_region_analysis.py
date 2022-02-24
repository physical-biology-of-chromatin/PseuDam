import pandas as pd
import matplotlib.pyplot as plt





df = pd.read_csv("/datas/nathan/vscode_nextflow/nextflow-nathan/results/GATC/sites_yeast_new.bed", 
                 header = None, 
                 sep = '\t')


sites = list()

region_size = 20000
at = 0.3085
gc = 0.1915
number_sites = region_size * gc * gc * at * at

window_sizes = [2500, 5000, 7500, 10000, 12500, 15000, 25000]


for window in window_sizes:

    sites = list()

    for row in df.itertuples():
            
        if (row[1] == "ref|NC_001135|" 
            and row[3] > (91324 - window) 
            and row[3] < (92418 + window)):
            
            gatc = row[3]
            
            if row[3] >= 91871:
                gatc += 1824
        
            sites.append(gatc)


    bin_list = list()
    values_list = list()
    j = 0
    w = 0

    for i in range (1,len(sites)):
        
        bin = ((sites[i] - sites[i-1]) / 2) + sites[i - 1]  
        
        
        
        if sites[i-1] >= 91871:
            j += 1
            value = w - j
            
        else:
            w +=1
            value = w
            

        bin_list.append(bin)
        
        values_list.append(value)
        
    data = list()

    for bin, value in zip(bin_list, values_list):
        
        for i in range(value):
            if (bin > 91871 
                and bin < 91871 + 1824):
                pass
            else:
                data.append(bin)


    plt.hist(data, bins = sites, edgecolor = "black")

    path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/sites analysis/"
            "region_analysis_window_size="
            + str(window*2))

    plt.savefig(path)
    
    plt.clf()