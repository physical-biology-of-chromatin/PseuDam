import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


df_noise = pd.read_csv("/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/Dam_alone.bed", 
                        header = None, 
                        sep = '\t')

df_signal = pd.read_csv("/datas/nathan/vscode_nextflow/nextflow-nathan/results/coverage/Dam_PC.bed", 
                          header = None, 
                          sep = '\t')

value_noise_list = list()
value_signal_list = list()
sites_list = list()
signal_norm_list = list()
bins_width_list = list()
j = 0
w = 0
i = 0
v = 0


signal_norm_ratio_list = list()
signal_norm_ratio_log_list = list()
signal_norm_sub_list = list()
signal_norm_reshape_list = list()
signal_noise_ratio_list = list()

for row_noise, row_signal in zip(df_noise.itertuples(), df_signal.itertuples()):
    
    if row_noise[7] == 0:
        j += 1
        
    if row_signal[7] == 0:
        w += 1
        
    if row_noise[7] == 0 and row_signal[7] == 0:
        v += 1
    
    if (row_noise[1] == "NT_033779.5" 
        and row_noise[2] > 14000000 
        and row_noise[3] < 16000000):
           
        value_noise = row_noise[7]
        value_noise_list.append(value_noise)

    if (row_signal[1] == "NT_033779.5" 
        and row_signal[2] > 14000000 
        and row_signal[3] < 16000000):
        
        site_signal = row_signal[2] + (row_signal[3] - row_signal[2]) 

        bin_width = row_signal[3] - row_signal[2]
        
        value_signal = row_signal[7]
        
        sites_list.append(site_signal)
        value_signal_list.append(value_signal)
        bins_width_list.append(bin_width)

    i +=1

print(f"there is {(j / i) * 100} de bins nuls dans le temoin ")
print(f"there is {(w / i) * 100} de bins nuls dans le PC ")
print(f"there is {(v / i) * 100} de bins nuls dans dans les 2 en meme temps")

mean_noise = df_noise[6].mean(skipna=True)
print(f"The mean covering of the noise sample is {mean_noise}")

mean_signal =  df_signal[6].mean(skipna=True)
print(f"The mean covering of the signal sample is {mean_noise}")

for value_signal, value_noise in zip(value_signal_list, value_noise_list):
    
    norm_value_signal = ((value_signal - min(value_signal_list))
                         /  (max(value_signal_list) - min(value_signal_list)))
    
    norm_value_noise = ((value_noise - min(value_noise_list))
                         /  (max(value_noise_list) - min(value_noise_list)))
    
    signal_norm_ratio = np.log2((value_signal + 1) / (value_noise + 1))
    signal_norm_ratio_log_list.append(signal_norm_ratio)
    
    signal_norm_ratio_log = (value_signal + 1) / (value_noise + 1)
    signal_norm_ratio_list.append(signal_norm_ratio_log)
    
    signal_norm_sub = value_signal - value_noise
    signal_norm_sub_list.append(signal_norm_sub)
    
    signal_norm_reshape = signal_norm_ratio * mean_signal
    signal_norm_reshape_list.append(signal_norm_reshape)
    
    signal_noise_ratio = ((value_signal + 1) - (value_noise + 1)) / (value_noise + 1)
    signal_noise_ratio_list.append(signal_noise_ratio)

plt.figure("ratio")
plt.bar(sites_list,
         height = signal_norm_ratio_list, 
         width = bins_width_list)
plt.ylabel("signal / noise")
plt.xlabel("positions")
plt.suptitle("barplot of the signal noise ratio")
path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
        +"analysis_test_data_ratio")
plt.savefig(path)


plt.figure("log2 ratio")
plt.bar(sites_list,
         height = signal_norm_ratio_log_list, 
         width = bins_width_list)
plt.ylabel("log2(signal / noise)")
plt.xlabel("positions")
plt.suptitle("barplot of the log2 of the signal noise ratio")
path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
        +"analysis_test_data_logratio")
plt.savefig(path)


plt.figure("signal noise ratio reshaped")
plt.bar(sites_list,
         height = signal_norm_reshape_list, 
         width = bins_width_list)
plt.ylabel("(signal / noise) * mean_signal")
plt.xlabel("positions")
plt.suptitle("barplot of the signal noise difference")
path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
        +"analysis_test_data_reshaped")
plt.savefig(path)


plt.figure("(signal - noise) / noise")
plt.bar(sites_list,
         height = signal_noise_ratio_list, 
         width = bins_width_list)
plt.ylabel("(signal - noise) / noise")
plt.xlabel("positions")
plt.suptitle("barplot of the signal noise ratio")
path = ("/datas/nathan/vscode_nextflow/nextflow-nathan/results/analysis/"
        +"analysis_test_data_signal-noise-ratio")
plt.savefig(path)