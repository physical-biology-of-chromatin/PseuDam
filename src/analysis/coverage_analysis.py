import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

file = "/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/coverage_per_frag/N_3RG1D1_EKDL220001506-1a_H2KMHDSX3_L4_cov.bed"

df = pd.read_csv(file,
                 sep = "\t",
                 header = None)

sample_name = re.search("N_(.+?)_EKDL", file).group(1)

chrom_name = "ref|NC_001135|"

df_currated = df[df[0] == chrom_name] 

x = np.arange(1, df_currated[2].max(), 10)

y_median = df_currated[6].median() * (x / x)
y_mean = df_currated[6].mean() * (x / x)

plt.plot(x, y_mean, color = "red", label = "mean")
plt.plot(x, y_median, color = "green", label = "median")
plt.legend()

plt.bar(((df_currated[2] - df_currated[1]) / 2) + df_currated[1],
        df_currated[6],
        df_currated[5])

plt.ylim(0,1.2)

print(df_currated[6].median())


plt.title(f"{sample_name} chromosome III \n mean = {round(df_currated[6].mean(), 3)} | median = {round(df_currated[6].median(), 10)}")
plt.xlabel("% coverage")
plt.ylabel("position")
plt.show()