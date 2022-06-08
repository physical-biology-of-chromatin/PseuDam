from packages.Dam_ID_utilities import chrom_rename
import re


####################################################################
#                                                                  #
#                         C elegans experiments                    #
#                                                                  #
####################################################################

def C_el():
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
        
def C_el_bin():

    epsilon = 10
    
    teta = 100

    chrom_names = {
        "chromosome_I" : "chrI",
        "chromosome_II" : "chrII",
        "chromosome_III" : "chrIII",
        "chromosome_IV" : "chrIV",
        "chromosome_V" : "chrV",
        "chromosome_X" : "chrX"
        }
        
    path = f"/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L4/binned_10_100/lfc_srf3_L4_binned_10_100.bedgraph"
    
    chrom_rename(path, chrom_names)
    
    path = f"/datas/nathan/Dam_ID_analysis/results/C_elegans/srf3/L4/binned_10_100/qval_srf3_L4_binned_10_100.bedgraph"
    
    chrom_rename(path, chrom_names)


def C_el_pipeline_nul():


    chrom_names = {
        "I" : "chrI",
        "II" : "chrII",
        "III" : "chrIII",
        "IV" : "chrIV",
        "V" : "chrV",
        "X" : "chrX"
        }

    path = "/datas/nathan/test_pipeline/srf3_L4/bedgraph_reference_pipeline.bedgraph"

    chrom_rename(path, chrom_names)



####################################################################
#                                                                  #
#                         DamC                                     #
#                                                                  #
####################################################################

def DamC():
    chrom_names = {
        "ref|NC_001133|" : "chrI",
        "ref|NC_001134|" : "chrII",
        "ref|NC_001135|" : "chrIII",
        "ref|NC_001136|" : "chrIV",
        "ref|NC_001137|" : "chrV",
        "ref|NC_001138|" : "chrVI",
        "ref|NC_001139|" : "chrVII",
        "ref|NC_001140|" : "chrVIII",
        "ref|NC_001141|" : "chrIX",
        "ref|NC_001142|" : "chrX",
        "ref|NC_001143|" : "chrXI",
        "ref|NC_001144|" : "chrXII",
        "ref|NC_001145|" : "chrXIII",
        "ref|NC_001146|" : "chrXIV",
        "ref|NC_001147|" : "chrXV",
        "ref|NC_001148|" : "chrXVI",
        "ref|NC_001224|" : "chrM",
        }


    chrom_names = {
        "chrom_I" : "chrI",
        "chrom_II" : "chrII",
        "chrom_III" : "chrIII",
        "chrom_IV" : "chrIV",
        "chrom_V" : "chrV",
        "chrom_VI" : "chrVI",
        "chrom_VII" : "chrVII",
        "chrom_VIII" : "chrVIII",
        "chrom_IX" : "chrIX",
        "chrom_X" : "chrX",
        "chrom_XI" : "chrXI",
        "chrom_XII" : "chrXII",
        "chrom_XIII" : "chrXIII",
        "chrom_XIV" : "chrXIV",
        "chrom_XV" : "chrXV",
        "chrom_XVI" : "chrXVI",
        "chrom_M" : "chrM",
        }


    chrom_names = {
        "ref|NC_001133|" : "chrI",
        "ref|NC_001134|" : "chrII",
        "ref|NC_001135|" : "chrIII",
        "ref|NC_001136|" : "chrIV",
        "ref|NC_001137|" : "chrV",
        "ref|NC_001138|" : "chrVI",
        "ref|NC_001139|" : "chrVII",
        "ref|NC_001140|" : "chrVIII",
        "ref|NC_001141|" : "chrIX",
        "ref|NC_001142|" : "chrX",
        "ref|NC_001143|" : "chrXI",
        "ref|NC_001144|" : "chrXII",
        "ref|NC_001145|" : "chrXIII",
        "ref|NC_001146|" : "chrXIV",
        "ref|NC_001147|" : "chrXV",
        "ref|NC_001148|" : "chrXVI",
        "ref|NC_001224|" : "chrM"
                                    
        }

            

            
    path = f"/datas/nathan/Dam_ID_analysis/results/total/counts_norm/3RG3D1_LFC.bedgraph"
        
    chrom_rename(path, chrom_names)
        
        
        
        
def sites_cel():

    chrom_names = {
        "ENA|BX284601|BX284601.5" : "chrI",
        "ENA|BX284602|BX284602.5" : "chrII",
        "ENA|BX284603|BX284603.4" : "chrIII",
        "ENA|BX284604|BX284604.4" : "chrIV",
        "ENA|BX284605|BX284605.5" : "chrV",
        "ENA|BX284606|BX284606.5" : "chrX"
        }

    
def renaming_bins_damc():
    
    names = {
   "chrom_I" : "ref|NC_001133|",
   "chrom_II" : "ref|NC_001134|", 
   "chrom_III" : "ref|NC_001135|",
   "chrom_IV" : "ref|NC_001136|", 
   "chrom_V" : "ref|NC_001137|",
   "chrom_VI" : "ref|NC_001138|", 
   "chrom_VII" : "ref|NC_001139|", 
   "chrom_VIII" : "ref|NC_001140|",
   "chrom_IX" : "ref|NC_001141|", 
   "chrom_X" : "ref|NC_001142|",
   "chrom_XI" : "ref|NC_001143|",
   "chrom_XII" : "ref|NC_001144|",
   "chrom_XIII" : "ref|NC_001145|",
   "chrom_XIV" : "ref|NC_001146|",
   "chrom_XV" : "ref|NC_001147|", 
   "chrom_XVI" : "ref|NC_001148|",
   "chrom_circular" : "ref|NC_001224|"
    }
    
    chrom_rename("/datas/nathan/vscode_nextflow/nextflow-nathan/data/sites/dam_total/sites_10_100.bed", names)
    
renaming_bins_damc()