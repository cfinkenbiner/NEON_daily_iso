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

import os.path
from os import path
import precip_iso_downscaling_scripts as iso_ds
import pandas as pd
import numpy as np

ROOTDIR = os.getcwd() # Home directory


### Create dataframe to write to and decide where to save it
df_d2H = pd.DataFrame(index=pd.to_datetime([])) 
df_d18O = pd.DataFrame(index=pd.to_datetime([])) 
df_meta = pd.DataFrame(index=pd.to_datetime([])) 

### Where would you like the daily P datasets saved to on your local computer?
out_dir = ROOTDIR


### Define location of NEON 30 min Precip + Isotope Data
file_dir = ROOTDIR + '/DATA/'


### NEON site list
site_list = ['ABBY','BARR','BART','BLAN','BONA','CLBJ','CPER','DCFS','DEJU','DELA',
             'DSNY','GRSM','GUAN','HARV','HEAL','JERC','JORN','KONA','KONZ','LAJA',
             'LENO','MLBS','MOAB','NIWO','NOGP','OAES','ONAQ','ORNL','OSBS','PUUM',
             'RMNP','SCBI','SERC','SJER','SOAP','SRER','STEI','STER','TALL','TEAK',
             'TOOL','TREE','UKFS','UNDE','WOOD','WREF','YELL']
            

### Iterate through site list        
for s in np.arange(len(site_list)):
    
    if path.exists(file_dir + 'PrecipData/' + str(site_list[s]) + 'PrecipData.xlsx'): 
        precip_30min_loc = file_dir + 'PrecipData/' + str(site_list[s]) + 'PrecipData.xlsx'
        iso_p_loc = file_dir + 'IsoData/' + str(site_list[s]) + 'IsoData.xlsx'
        
        precip_filter = 0.25 ## remove small P amts --> edit this if you desire
        
        ### run downscaling method (Finkenbiner et al. 2021)
        site_name = site_list[0]
        site_dailyP = iso_ds.read_p_data(site_list[s], precip_30min_loc, precip_filter, iso_p_loc, out_dir, residual_corr=True)
        
        
        ### save output to dataframe
        if 'd2H' in site_dailyP.columns:
            
            df_d2H[site_list[s]] = site_dailyP['d2H']
            df_d18O[site_list[s]] = site_dailyP['d18O']
            
        else:
     
            df_d2H[site_list[s]] = np.nan
            df_d18O[site_list[s]] = np.nan
            
        ### write metadata
        site_dailyP['meta'] = 1.0
        for i in np.arange(len(site_dailyP['Total P (unfiltered)'])):
            
            if (site_dailyP['Total P (unfiltered)'].iloc[i] >= precip_filter):
                site_dailyP['meta'].iloc[i] = 0.0
                
            elif ((0 < site_dailyP['Total P (unfiltered)'].iloc[i]) & 
                  (site_dailyP['Total P (unfiltered)'].iloc[i] < precip_filter)):
                site_dailyP['meta'].iloc[i] = 2.0
                
            else:
                site_dailyP['meta'].iloc[i] = 1.0
        
        df_meta[site_list[s]] = site_dailyP['meta']
        
    else: 
        print('Skipping ' + site_list[s] + ' because there is no precipitation or isotope data available in directory...')
        
        df_d2H[site_list[s]] = np.nan
        df_d18O[site_list[s]] = np.nan

df_d2H = df_d2H.rename_axis('Date')
df_d18O = df_d18O.rename_axis('Date')
df_meta = df_meta.rename_axis('Date')
df_meta = df_meta.fillna(1)

### export csv files
df_d2H.to_csv(out_dir + '/daily_p_d2H.csv')
df_d18O.to_csv(out_dir + '/daily_p_d18O.csv')
df_meta.to_csv(out_dir + '/daily_p_metadata.csv')
