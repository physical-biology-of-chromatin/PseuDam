import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from time import time
import timeit
from matplotlib import cm

start = timeit.default_timer()

input0 = [[3,80000,100000]]
res = 6400
val = 'resolutions'

with h5py.File("/datas/nathan/vscode_nextflow/nextflow-nathan/data/HiC/GSE151553_A364_merged.mcool", "r") as h5f:
    
    # Get a h5py dataset object
    print("possible resolutions=",list(h5f[val]))
    
    chromid = h5f[val][str(res)]['bins']['chrom']
    
    bin_start = h5f[val][str(res)]['bins']['start']
    
    bin_end = h5f[val][str(res)]['bins']['end']
    
    bin_weigth = h5f[val][str(res)]['bins']['weight']
    
    chroms_name = h5f[val][str(res)]['chroms']['name']
    
    chroms_length = h5f[val][str(res)]['chroms']['length']
    
    chrom_offset =h5f[val][str(res)]['indexes']['chrom_offset']
    
    bin1_offset = h5f[val][str(res)]['indexes']['bin1_offset']
    
    bin1_t = h5f[val][str(res)]['pixels']['bin1_id']
    
    bin2_t = h5f[val][str(res)]['pixels']['bin2_id']
    
    count_t=h5f[val][str(res)]['pixels']['count']
    
    
    for inp in range(len(input0)):
        
        ch = input0[inp][0]
        start0 = input0[inp][1]
        end0   = input0[inp][2]
        start_bin = int(start0/res)
        end_bin = int(end0/res)
        length_matrix = end_bin-start_bin+1
        length_whole_chromosome = chroms_length[ch - 1]
        HiC = np.zeros((length_whole_chromosome, length_matrix))
        
        """
        for i in range(len(chroms_name)):
            
            if (chroms_name[i].decode('UTF-8')==ch):
                
                chID = i
                break
        """
        
        
        chID = input0[0][0] - 1
        start_bin = chrom_offset[chID]+start_bin
        end_bin = chrom_offset[chID]+end_bin
        print("start bin=",start_bin,",end bin=",end_bin,",total bins=",end_bin-start_bin)
        arr_start = bin1_offset[start_bin]
        arr_end = bin1_offset[end_bin]
        
        whole_start = bin1_offset[0]
        whole_stop = bin1_offset[int(length_whole_chromosome/res)]
        
        for i in range(arr_start,arr_end+1):
            

            bin2=bin2_t[i]
            
            count=count_t[i]
            
            if (bin2<=end_bin):
                
                for j in range(whole_start, whole_stop):
                    
                    bin1=bin1_t[j]
                
            
                    wt = bin_weigth[bin1]*bin_weigth[bin2]*count
                
                    HiC[bin1, bin2 - start_bin] = wt

        
        fig, ax = plt.subplots(1,len(input0),
                               sharex=False,
                               sharey=False)

        ax.imshow(np.log2(HiC),cmap='YlOrRd')
        ax.set_title('resolution='+str(res/1000)+'kb')
        ax.set_xlabel(f"{ch} : {str(start0/1000000)} - {str(end0/1000000)}Mb")
        ax.set_ylabel(f"{ch} : {str(start0/1000000)} - {str(end0/1000000)}Mb")
        ax.set_xticks(np.linspace(0, length_whole_chromosome, 11))
        
        stop = timeit.default_timer()

        print('Time: ', stop - start)  
        
        plt.show()
        
        a = 1