#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 18:08:26 2021

@author: Catie Finkenbiner
@email: c.finkenbiner@gmail.com
"""

'''
THIS SCRIPT IS JUST FOR INVESTIGATING ISOTOPE VARIANCE!!!
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


### Where would you like the daily P datasets saved to on your local computer?
out_dir = ROOTDIR


### Define location of NEON 30 min Precip + Isotope Data
file_dir = '/Users/catiefinkenbiner/Documents/NEON/NEON_daily_isotope_dataset/NEON_daily_iso/DATA/'


### NEON site list
site_list = ['ABBY','BARR','BART','BLAN','BONA','CLBJ','CPER','DCFS','DEJU','DELA',
             'DSNY','GRSM','GUAN','HARV','HEAL','JERC','JORN','KONA','KONZ','LAJA',
             'LENO','MLBS','MOAB','NIWO','NOGP','OAES','ONAQ','ORNL','OSBS','PUUM',
             'RMNP','SCBI','SERC','SJER','SOAP','SRER','STEI','STER','TALL','TEAK',
             'TOOL','TREE','UKFS','UNDE','WOOD','WREF','YELL']

for ensemble in np.arange(1, 11):
    
    ### Iterate through site list        
    for s in np.arange(len(site_list)):
        if path.exists(file_dir + 'PrecipData/' + str(site_list[s]) + 'PrecipData.xlsx'): 
            precip_30min_loc = file_dir + 'PrecipData/' + str(site_list[s]) + 'PrecipData.xlsx'
            iso_p_loc = file_dir + 'IsoData/' + str(site_list[s]) + 'IsoData.xlsx'
            
            precip_filter = 0.25 ## [mm] remove small P amts/trace events
            
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
    
        else: 
            df_d2H[site_list[s]] = np.nan
            df_d18O[site_list[s]] = np.nan
    
    
    ### export csv files
    df_d2H.to_csv(out_dir + '/daily_d2H_'+str(ensemble)+'.csv')
    df_d18O.to_csv(out_dir + '/daily_d18O_'+str(ensemble)+'.csv')
