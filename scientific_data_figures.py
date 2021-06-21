#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 15:36:57 2021

@author: catiefinkenbiner
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOTDIR = os.getcwd() # Home directory

#'/Users/catiefinkenbiner/Documents/NEON/NEON_daily_isotope_dataset/NEON_daily_iso/'
df_p_d2H = pd.read_csv(ROOTDIR + '/daily_d2H.csv', parse_dates=True, index_col=0)
df_p_d18O = pd.read_csv(ROOTDIR + '/daily_d18O.csv', parse_dates=True, index_col=0)


''' Figure 1. Time series of generated datasets '''
### location of NOEN 30 min Precip + Biweekly Isotope Data
df_biweekly_onaq = pd.read_excel(ROOTDIR + '/DATA/IsoData/ONAQIsoData.xlsx', parse_dates=True, index_col='collectDate')
df_biweekly_wref = pd.read_excel(ROOTDIR + '/DATA/IsoData/WREFIsoData.xlsx', parse_dates=True, index_col='collectDate')

### location of NOEN daily Precip Data (Lindsey's work)
df_daily_onaq = pd.read_csv(ROOTDIR + '/OUTPUT/ONAQ_daily_timeseries.csv', parse_dates=True, index_col=0)
df_daily_wref = pd.read_csv(ROOTDIR + '/OUTPUT/WREF_daily_timeseries.csv', parse_dates=True, index_col=0)


fig1, ax1 = plt.subplots(nrows=5, ncols=2, figsize=(12,15))

### Subplot [0,0] -- ONAQ
ax1[0,0].plot(df_p_d2H.index, df_p_d2H['ONAQ'], 'o', color='k', label='Daily')
ax1[0,0].plot(df_biweekly_onaq.index, df_biweekly_onaq['d2HWater'], '<', color='red', label='Biweekly')
ax00b = ax1[0,0].twinx()
ax00b.bar(df_daily_onaq.index, df_daily_onaq['Total P'], color='b', width=1, alpha=0.5)

ax1[0,0].set_title(r'ONAQ')
ax1[0,0].set_ylabel(r'$\delta^{2}H$ (‰)')

### Subplot [1,0]
ax1[1,0].plot(df_p_d18O.index, df_p_d18O['ONAQ'], 'o', color='k', label='Daily')
ax1[1,0].plot(df_biweekly_onaq.index, df_biweekly_onaq['d18OWater'], '<', color='red', label='Biweekly')
ax10b = ax1[1,0].twinx()
ax10b.bar(df_daily_onaq.index, df_daily_onaq['Total P'], color='b', width=1, alpha=0.5)

ax1[1,0].set_ylabel(r'$\delta^{18}O$ (‰)')


### Subplot [0,1] -- WREF
ax1[0,1].plot(df_p_d2H.index, df_p_d2H['WREF'], 'o', color='k', label='Daily')
ax1[0,1].plot(df_biweekly_wref.index, df_biweekly_wref['d2HWater'], '<', color='red', label='Biweekly')
ax01b = ax1[0,1].twinx()
ax01b.bar(df_daily_wref.index, df_daily_wref['Total P'], color='b', width=1, alpha=0.5)

ax1[0,1].set_title(r'WREF')
ax01b.set_ylabel('Precpitation (cm)', color='b')

### Subplot [1,1]
ax1[1,1].plot(df_p_d18O.index, df_p_d18O['WREF'], 'o', color='k', label='Daily')
ax1[1,1].plot(df_biweekly_wref.index, df_biweekly_wref['d18OWater'], '<', color='red', label='Biweekly')
ax11b = ax1[1,1].twinx()
ax11b.bar(df_daily_wref.index, df_daily_wref['Total P'], color='b', width=1, alpha=0.5)

ax11b.set_ylabel('Precpitation (cm)', color='b')

plt.tight_layout()
plt.show() ; plt.close()


''' Figure 2. Time series of generated datasets '''