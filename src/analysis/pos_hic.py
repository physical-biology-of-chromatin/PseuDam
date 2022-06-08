
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from time import time
from matplotlib import cm

input0 = ['chr1',650000,7270000]

res = 10000
val = 'resolutions'
normal = 'KR'

with h5py.File("../mouse_ES_micro-c/4DNFINNZDDXV.mcool", "r") as h5f:
    
    print('possible resolutions=',list(h5f[val]))
    
    chromid = h5f[val][str(res)]['bins']['chrom']
    
    bin_start = h5f[val][str(res)]['bins']['start']
    
    bin_end = h5f[val][str(res)]['bins']['end']
    
    bin_weigth = h5f[val][str(res)]['bins'][normal]
    
    bin_weight = pd.DataFrame({'KR':bin_weigth})
    
    chroms_name = h5f[val][str(res)]['chroms']['name']
    
    chroms_length = h5f[val][str(res)]['chroms']['length']
    
    chrom_offset =h5f[val][str(res)]['indexes']['chrom_offset']
    
    bin1_offset = h5f[val][str(res)]['indexes']['bin1_offset']
    
    bin1_t = h5f[val][str(res)]['pixels']['bin1_id']
    
    bin2_t = h5f[val][str(res)]['pixels']['bin2_id']
    
    count_t=h5f[val][str(res)]['pixels']['count']
    
    ch = input0[0]
    
    start0 = input0[1]
    
    end0   = input0[2]
    
    start_bin = int(start0/res)
    
    end_bin = int(end0/res)
    
    length_matrix = end_bin-start_bin+1
    
    
    for i in range(len(chroms_name)):
        if (chroms_name[i].decode('UTF-8')==ch):
            chID = i
            break
        
        
    start_bin = chrom_offset[chID]+start_bin
    end_bin = chrom_offset[chID]+end_bin
    
    arr_start = bin1_offset[start_bin]
    arr_end = bin1_offset[end_bin]+1
    
    df = pd.DataFrame({'bin1':bin1_t[arr_start:arr_end],'bin2':bin2_t[arr_start:arr_end],'raw':count_t[arr_start:arr_end]})
    
    sub_df2 = df[(df['bin1']>=start_bin)&(df['bin1']<=end_bin)]
    sub_df = sub_df2[(sub_df2['bin2']>=start_bin)&(sub_df2['bin2']<=end_bin)]
    sub_df.loc[:,'wt1'] = sub_df['bin1']
    sub_df = sub_df.set_index('wt1')
    sub_df = pd.merge(sub_df,bin_weight, left_index=True, right_index=True, how='left')
    sub_df.loc[:,'wt2'] = sub_df['bin2']
    sub_df = sub_df.set_index('wt2')
    sub_df = sub_df.rename(columns = {'KR': 'wt1'}, inplace = False)
    sub_df = pd.merge(sub_df,bin_weight, left_index=True, right_index=True, how='left')
    sub_df = sub_df.rename(columns = {'KR': 'wt2'}, inplace = False)
    sub_df = sub_df.sort_values(by=['bin1','bin2'])
    sub_df['norm'] = sub_df['wt1']*sub_df['wt2']*sub_df['raw']
    sub_df = sub_df.drop(['wt1','wt2'], 1)
    sub_df['gd'] = sub_df['bin2'] - sub_df['bin1']
    sub_df1 = sub_df.drop(['bin1','bin2','raw'], 1)
    exp_df = sub_df1.groupby('gd').sum()
    exp_df = exp_df.reset_index()
    exp_df['norm'] = exp_df['norm']/(length_matrix-exp_df['gd'])
    exp_df['gd'] = exp_df['gd']*(res/1000)
    
    plt.plot(exp_df['gd'][1:], exp_df['norm'][1:]/exp_df['norm'][1])
    
    exp_df = exp_df.rename(columns={"gd":"genomic_distance(kb)"})
    exp_df = exp_df.rename(columns={"norm":"contact_frequency"})
    exp_df.to_csv('contact_frequency_'+ch+'_'+str(start0/1000000)+'-'+str(end0/1000000)+'Mb_resolution'+str(res/1000)+'kb.txt', sep ='\t', encoding='utf-8',index=False)
    
    print(exp_df)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.title(ch+'_'+str(start0/1000000)+'-'+str(end0/1000000)+'Mb')
    plt.xlabel('genomic distance (kb)')
    plt.ylabel('contact probability')
    plt.show()