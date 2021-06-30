#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 11:21:10 2021

@author: catiefinkenbiner
"""

### Average std per site

import pandas as pd
import numpy as np
import os.path
from os import path

ROOTDIR = os.getcwd() # Home directory
file_dir = ROOTDIR + '/isoP_stds/'

site_list = ['ABBY','BARR','BART','BLAN','BONA','CLBJ','CPER','DCFS','DEJU','DELA',
             'DSNY','GRSM','GUAN','HARV','HEAL','JERC','JORN','KONA','KONZ','LAJA',
             'LENO','MLBS','MOAB','NIWO','NOGP','OAES','ONAQ','ORNL','OSBS','PUUM',
             'RMNP','SCBI','SERC','SJER','SOAP','SRER','STEI','STER','TALL','TEAK',
             'TOOL','TREE','UKFS','UNDE','WOOD','WREF','YELL']

df_avg_std_d2H = pd.DataFrame(columns = site_list)
df_avg_std_d18O = pd.DataFrame(columns = site_list)

for s in site_list:
    df_temp1 = pd.DataFrame()
    df_temp2 = pd.DataFrame()
    
    for i in np.arange(1,11):
            
        d2H = pd.read_csv(file_dir + 'daily_d2H_' + str(i) + '.csv', index_col=0, parse_dates=True)
        d18O = pd.read_csv(file_dir + 'daily_d18O_' + str(i) + '.csv', index_col=0, parse_dates=True)
        
        df_temp1[str(i)] = d2H[s]
        df_temp2[str(i)] = d18O[s]
        
    df_avg_std_d2H[s] = df_temp1.std(axis=1) / np.sqrt(10)
    df_avg_std_d18O[s] = df_temp2.std(axis=1) / np.sqrt(10)
    
df_avg_std_d2H.to_csv(ROOTDIR + '/d2H_Pstds.csv')
df_avg_std_d18O.to_csv(ROOTDIR + '/d18O_Pstds.csv')
