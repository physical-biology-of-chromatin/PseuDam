from packages.Dam_ID_utilities import chrom_rename

tissue_list = ["dpy7", "srf3"]

stage_list = ["L2", "L4"]

chrom_names = {
    "ENA|BX284601|BX284601.5" : "chrI",
    "ENA|BX284602|BX284602.5" : "chrII",
    "ENA|BX284603|BX284603.4" : "chrIII",
    "ENA|BX284604|BX284604.4" : "chrIV",
    "ENA|BX284605|BX284605.5" : "chrV",
    "ENA|BX284606|BX284606.5" : "chrX"
    }


for tissue in tissue_list:
    
    for stage in stage_list:
        
        path = f"/datas/nathan/Dam_ID_analysis/results/C_elegans/{tissue}/{stage}/lfc_{tissue}_{stage}.bedgraph"
        
        chrom_rename(path, chrom_names)
        
        path = f"/datas/nathan/Dam_ID_analysis/results/C_elegans/{tissue}/{stage}/qval_{tissue}_{stage}.bedgraph"
        
        chrom_rename(path, chrom_names)