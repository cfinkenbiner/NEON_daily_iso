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

### location of data
file_dir = '/Users/catiefinkenbiner/Documents/NEON/Downscaled_Product/Temporal_Downscaling_NEON_Pisotopes/DATA/'

### site list
site_list = 'ONAQ'

### location of 30 min Precip data
precip_30min_loc = file_dir+'PrecipData/'+str(site_list)+'PrecipData.xlsx'
iso_p_loc = file_dir+'IsoData/'+str(site_list)+'IsoData.xlsx'

precip_filter = 0.05 ## remove small P amts

iso_ds.read_p_data(precip_30min_loc, precip_filter, iso_p_loc)
