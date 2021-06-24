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
import cartopy
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes

ROOTDIR = os.getcwd() # Home directory
#'/Users/catiefinkenbiner/Documents/NEON/NEON_daily_isotope_dataset/NEON_daily_iso/'

df_p_d2H = pd.read_csv(ROOTDIR + '/daily_d2H.csv', parse_dates=True, index_col=0)
df_p_d18O = pd.read_csv(ROOTDIR + '/daily_d18O.csv', parse_dates=True, index_col=0)


''' Figure 1. Time series of generated datasets '''
make_Fig1 = True
if make_Fig1 == True:
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
    ax1[0,0].legend()
    
    ### Subplot [1,0]
    ax1[1,0].plot(df_p_d18O.index, df_p_d18O['ONAQ'], 'o', color='k', label='Daily')
    ax1[1,0].plot(df_biweekly_onaq.index, df_biweekly_onaq['d18OWater'], '<', color='red', label='Biweekly')
    ax10b = ax1[1,0].twinx()
    ax10b.bar(df_daily_onaq.index, df_daily_onaq['Total P'], color='b', width=1, alpha=0.5)
    
    ax1[1,0].set_ylabel(r'$\delta^{18}O$ (‰)')
    ax1[1,0].legend()
    
    
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
    plt.savefig(ROOTDIR+'/Fig1.png', dpi=500)
    plt.show() ; plt.close()


''' Figure 2. Summary maps of different data products '''
make_Fig2 = True
if make_Fig2 == True:
    df_p_d2H = pd.read_csv(ROOTDIR + '/daily_d2H.csv', parse_dates=True, index_col=0)
    df_p_d18O = pd.read_csv(ROOTDIR + '/daily_d18O.csv', parse_dates=True, index_col=0)
    df_summary = pd.read_csv(ROOTDIR + '/Site_Summary_Table.csv') # for site lat/lon
    
    lat = []
    lon = []
    
    for i in df_p_d2H.columns:
        row = df_summary[df_summary['Site ID'].str.match(i)]
        
        if len(row) > 0:
            lat.append(float(row['Latitude'].values))
            lon.append(float(row['Longitude'].values))
        else:
            lat.append(-9999)
            lon.append(-9999)
    
    def create_map(lat, lon, values, label, current_subplot, ii):
        ax = current_subplot
        ax.set_extent([-170, -55, 10, 70], ccrs.Geodetic())
        gl = ax.gridlines(linestyle='--', draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right = False
        ax.coastlines()
        
        im = ax.scatter(lon, lat, 
                        alpha=0.85, 
                        s=100,
                        c = values,
                        cmap = plt.get_cmap("viridis"), 
                        transform = ccrs.PlateCarree())
        axgr.cbar_axes[ii].colorbar(im)
        axgr[ii].set_title(label)
        
    projection = ccrs.PlateCarree()
    axes_class = (GeoAxes, dict(map_projection = projection))   
    
    fig = plt.figure(1,figsize = (12,7))
    axgr = AxesGrid(fig, 111, axes_class = axes_class, 
                        nrows_ncols=(1,2), 
                        axes_pad=0.6,
                        cbar_location = 'right',
                        cbar_mode = 'each',
                        cbar_pad = 0.3,
                        cbar_size = '3%',
                        label_mode = '')
    
    create_map(np.array(lat), np.array(lon),
               df_p_d2H.mean(axis=0),
               'Average, $\delta^{2}$H (‰)', 
               axgr[0], 0)
    axgr[0].text(-165, 15, 'a)', fontsize=14)
    
    create_map(np.array(lat), np.array(lon),
               df_p_d18O.mean(axis=0),
               'Average, $\delta^{18}$O (‰)', 
               axgr[1], 1)
    axgr[1].text(-165, 15, 'b)', fontsize=14)
    
    plt.savefig(ROOTDIR+'/Fig2.png', dpi=500)
    plt.show()


''' Figure 3. Average std by site map '''
make_Fig3 = True
if make_Fig3 == True:
    df_avg_std_d2H = pd.read_csv(ROOTDIR + '/d2H_Pstds.csv', index_col=0)
    df_avg_std_d18O = pd.read_csv(ROOTDIR + '/d18O_Pstds.csv', index_col=0)
    df_summary = pd.read_csv(ROOTDIR + '/Site_Summary_Table.csv') # for site lat/lon
    
    lat = []
    lon = []
    
    for i in df_avg_std_d2H.columns:
        row = df_summary[df_summary['Site ID'].str.match(i)]
        
        if len(row) > 0:
            lat.append(float(row['Latitude'].values))
            lon.append(float(row['Longitude'].values))
        else:
            lat.append(-9999)
            lon.append(-9999)
    
    def create_map(lat, lon, values, label, current_subplot, ii):
        ax = current_subplot
        ax.set_extent([-170, -55, 10, 70], ccrs.Geodetic())
        gl = ax.gridlines(linestyle='--', draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right = False
        ax.coastlines()
        
        im = ax.scatter(lon, lat, 
                        alpha=0.85, 
                        s=100,
                        c = values,
                        cmap = plt.get_cmap("viridis"), 
                        transform = ccrs.PlateCarree())
        axgr.cbar_axes[ii].colorbar(im)
        axgr[ii].set_title(label)
        
    projection = ccrs.PlateCarree()
    axes_class = (GeoAxes, dict(map_projection = projection))   
    
    fig = plt.figure(1,figsize = (12,7))
    axgr = AxesGrid(fig, 111, axes_class = axes_class, 
                        nrows_ncols=(1,2), 
                        axes_pad=0.6,
                        cbar_location = 'right',
                        cbar_mode = 'each',
                        cbar_pad = 0.3,
                        cbar_size = '3%',
                        label_mode = '')
    
    create_map(np.array(lat), np.array(lon),
               df_avg_std_d2H.mean(axis=0),
               'Average Standard Deviation, $\delta^{2}$H (‰)', 
               axgr[0], 0)
    axgr[0].text(-165, 15, 'a)', fontsize=14)
    
    create_map(np.array(lat), np.array(lon),
               df_avg_std_d18O.mean(axis=0),
               'Average Standard Deviation, $\delta^{18}$O (‰)', 
               axgr[1], 1)
    axgr[1].text(-165, 15, 'b)', fontsize=14)
    
    plt.savefig(ROOTDIR+'/Fig3.png', dpi=500)
    plt.show()






