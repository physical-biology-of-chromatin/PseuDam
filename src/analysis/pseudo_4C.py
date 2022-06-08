
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from time import time
from matplotlib import cm




# [[chrom, start, stop]]
input0 = [['chrIII', 0, 12000000]]



res = 400
val = 'resolutions'

region = [80000, 100000, "chromIII"]
start_input = region[0]
end_input = region[1]

with h5py.File("/datas/nathan/Dam_ID_analysis/data/GSE151553_A364_merged.mcool", "r") as h5f:
    
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
        HiC = np.zeros((length_matrix,length_matrix))
        
        # Getting the desired region bins
        region_start = int(start_input/res)
        region_end = int(end_input/res)
        
        for i in range(len(chroms_name)):
            
            if (chroms_name[i].decode('UTF-8')==ch):
                chID = i
                break

        
        print("start bin=",start_bin,",end bin=",end_bin,",total bins=",end_bin-start_bin)
        
        arr_start = bin1_offset[start_bin]
        arr_end = bin1_offset[end_bin]
        
        
        # Calculating the offset for the deisred region
        
        region_start = chrom_offset[chID]+region_start
        region_end = chrom_offset[chID]+region_end
        
        arr_region_start = bin1_offset[region_start]
        arr_region_end = bin1_offset[region_end]
        
        for i in range(arr_start,arr_end+1):
            
            bin1=bin1_t[i]
            bin2=bin2_t[i]
            count=count_t[i]
            
            if i == arr_region_start:
                HiC_start = bin1 - start_bin
                
            if i == arr_region_end:
                HiC_end = bin1 - start_bin

            
            if (bin2<=end_bin):
                
                wt = bin_weigth[bin1]*bin_weigth[bin2]*count
                HiC[bin1-start_bin,bin2-start_bin] = wt
                HiC[bin2-start_bin,bin1-start_bin] = wt
        
             
        plt.figure("heatmap")
        HiC_map = np.log2(HiC[HiC_start:HiC_end])
        plt.imshow(HiC_map,cmap='YlOrRd')
        plt.title('resolution='+str(res/1000)+'kb')
        plt.xlabel(ch+':'+str(start0/1000)+'-'+str(end0/1000)+'Kb')
        plt.ylabel(ch+':'+str(start_input/1000)+'-'+str(end_input/1000)+'Kb')
        plt.xticks([])
        plt.yticks([])
        cbar = plt.colorbar()
        cbar.set_label('contact frequency (Log2)')
        
        plt.figure("graphe")
        
        
        HiC = HiC[HiC_start:HiC_end]
        
        np.savetxt("4C_400.txt.gz", HiC, delimiter = ",")
        
        
        plt.plot([np.nanmean(col) for col in HiC.T])
        
    plt.show()