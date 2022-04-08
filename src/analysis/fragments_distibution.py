import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/N_2RG3_EKDL220001503-1a_H2C2YDSX3_L4_cov.bed",
                 sep = "\t",
                 header = None)

print(df[5].median())

bin_size = 5

plt.hist(df[5],
         bins = np.arange(0, 3000, bin_size))

plt.xticks(np.arange(0, 3000, 200))
plt.xlim(0, 3000)
plt.title(f"distribution of the size of the gatc sites in the yeast genome\n bin size = {bin_size}")

plt.show()
