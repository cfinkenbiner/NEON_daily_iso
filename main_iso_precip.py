#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 18:08:26 2021

@author: Catie Finkenbiner
@email: c.finkenbiner@gmail.com
"""

'''
1. Read in NEON precip. data and measured isotope datasets
2. Aggregate to biweekly scales
3. Calculate stats and estimate downscaled stats
4. Generate new time series + perform residual correction
'''

import precip_iso_downscaling_scripts as iso_ds
import pandas as pd
import numpy as np

### Create dataframe to write to and decide where to save it
df_d2H = pd.DataFrame(index=pd.to_datetime([])) # empty to begin
df_d18O = pd.DataFrame(index=pd.to_datetime([])) # empty to begin
out_dir = '/Users/catiefinkenbiner/Documents/NEON/NEON_daily_isotope_dataset/NEON_daily_iso/'


### location of NOEN 30 min Precip + Isotope Data
file_dir = '/Users/catiefinkenbiner/Documents/NEON/Downscaled_Product/Temporal_Downscaling_NEON_Pisotopes/DATA/'


### site list
site_list = ["DELA", "LENO", "TALL", "BARR", "BONA", "HEAL", "TOOL", "SRER","SYCA", "SJER", "ARIK", "CPER", "NIWO", "RMNP",
             "STER", "OSBS", "JERC", "PUUM", "KONZ", "UKFS", "SERC", "HARV", "UNDE", "BART", "NOGP", "WOOD", "BLUE", "OAES", "CUPE",
              "GUAN", "GUIL", "GRSM", "ORNL", "CLBJ", "PRIN", "MOAB", "ONAQ", "REDB", "BLAN", "MLBS", "SCBI", "WREF", "STEI", "YELL"]

for s in np.arange(len(site_list)):
    
    precip_30min_loc = file_dir+'PrecipData/'+str(site_list[s])+'PrecipData.xlsx'
    iso_p_loc = file_dir+'IsoData/'+str(site_list[s])+'IsoData.xlsx'
    
    precip_filter = 0.25 ## remove small P amts --> edit this if you desire
    
    
    ### runs downscaling method
    site_name = site_list[0]
    site_dailyP = iso_ds.read_p_data(site_list[s], precip_30min_loc, precip_filter, iso_p_loc, out_dir, residual_corr=True)
    
    if 'd2H' in site_dailyP.columns:
        
        df_d2H[site_list[s]] = site_dailyP['d2H']
        df_d18O[site_list[s]] = site_dailyP['d18O']
    else:
 
        df_d2H[site_list[s]] = np.nan
        df_d18O[site_list[s]] = np.nan


df_d2H.to_csv(out_dir+'daily_d2H.csv')
df_d18O.to_csv(out_dir+'daily_d18O.csv')
