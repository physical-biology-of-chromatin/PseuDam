import h5py
import numpy as np
import matplotlib.pyplot as plt


"""
Plots the desired region from a given hic h5 file
"""


# [[chrom name, start, stop]]
input0 = [["chrIII", 0, 320000]]

res = 400

val = "resolutions"

h5_path = "/datas/nathan/Dam_ID_analysis/data/GSE151553_A364_merged.mcool"

with h5py.File(h5_path, "r") as h5f:
    
    # Get a h5py dataset object
    
    print("possible resolutions=",list(h5f[val]))
    
    chromid = h5f[val][str(res)]["bins"]["chrom"]
    
    bin_start = h5f[val][str(res)]["bins"]["start"]
    
    bin_end = h5f[val][str(res)]["bins"]["end"]
    
    bin_weigth = h5f[val][str(res)]["bins"]["weight"]
    
    chroms_name = h5f[val][str(res)]["chroms"]["name"]
    
    chroms_length = h5f[val][str(res)]["chroms"]["length"]
    
    chrom_offset =h5f[val][str(res)]["indexes"]["chrom_offset"]
    
    bin1_offset = h5f[val][str(res)]["indexes"]["bin1_offset"]
    
    bin1_t = h5f[val][str(res)]["pixels"]["bin1_id"]
    
    bin2_t = h5f[val][str(res)]["pixels"]["bin2_id"]
    
    count_t=h5f[val][str(res)]["pixels"]["count"]
    
    for offset in chrom_offset:
        print(offset)
    
    
    for inp in range(len(input0)):
        
        ch = input0[inp][0]
        start0 = input0[inp][1]
        end0   = input0[inp][2]
        start_bin = int(start0/res)
        end_bin = int(end0/res)
        length_matrix = end_bin-start_bin+1
        HiC = np.zeros((length_matrix,length_matrix))
        
        # Finds the chrom from the name
        for i in range(len(chroms_name)):
            
            if (chroms_name[i].decode("UTF-8")==ch):
                chID = i
                break
            
        start_bin = chrom_offset[chID]+start_bin
        end_bin = chrom_offset[chID]+end_bin
        print("start bin=",start_bin,",end bin=",end_bin,",total bins=",end_bin-start_bin)
        arr_start = bin1_offset[start_bin]
        arr_end = bin1_offset[end_bin]
        
        for i in range(arr_start,arr_end+1):
            
            bin1=bin1_t[i]
            bin2=bin2_t[i]
            count=count_t[i]
            
            if (bin2<=end_bin):
                
                wt = bin_weigth[bin1]*bin_weigth[bin2]*count
                HiC[bin1-start_bin,bin2-start_bin] = wt
                HiC[bin2-start_bin,bin1-start_bin] = wt
       
                
        plt.subplot(1,len(input0),inp+1)
        
        plt.imshow(np.log2(HiC),cmap="YlOrRd")
        
        plt.title("resolution="+str(res/1000)+"kb")
        
        plt.xlabel(ch+":"+str(start0/1000000)+"-"+str(end0/1000000)+"Mb")
        plt.ylabel(ch+":"+str(start0/1000000)+"-"+str(end0/1000000)+"Mb")
        
        plt.xticks([])
        plt.yticks([])
        
        cbar = plt.colorbar()
        cbar.set_label("contact frequency (Log2)")
        
    plt.show()