from more_itertools import only
import pandas as pd
import matplotlib.pyplot as plt

import re
from os import listdir
from os.path import isfile, join


mypath = "/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites"

onlyfiles = [mypath + "/" + f for f in listdir(mypath) if isfile(join(mypath, f))]


list_length = list()
list_labels = list()

for i, file in enumerate(onlyfiles):
    df = pd.read_csv(file, header = None, sep = "\t")
    
    label = (re.search("binned_(.+?).bed", file).group(1))
    
    length = df[2] - df[1]
    mean = round(length.mean())
    std = round(length.std())
    median = round(length.median())
    
    list_length.append(length)
    list_labels.append(label)
    
    
    print(f"{file}\n"
          f"mean = {mean}\n"
          f"std = {std}\n"
          f"median = {median}\n")


plt.violinplot(list_length)

plt.ylabel("Bins_length")
plt.xlabel("Binnings")

plt.xticks(range(1, len(list_labels)+1),
           labels = list_labels)

plt.show()
plt.clf()


path = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/Dam_ID/c_elegans/sites/sites.bed"
df = pd.read_csv(path, header = None, sep = "\t")

length = df[2] - df[1]
mean = length.mean()
median = length.median()
std = length.std()




print(f"mean = {mean}\n"
      f"std = {std}\n"
      f"median = {median}\n")
