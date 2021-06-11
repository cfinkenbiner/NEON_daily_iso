#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 18:24:31 2021

@author: catiefinkenbiner
"""

import numpy as np
import pandas as pd
from datetime import datetime


''' Read time and calculate fractional days '''
def changeTimes_change_Pdata(siteDF):
    timeStrs = siteDF['endDateTime']
    newTimes = [] ; fracyr = []
  
    for t in timeStrs:
        cTime = datetime.strptime(str(t),"%Y-%m-%d %H:%M:%S")
        newTimes.append(cTime)
        fracyr.append(float(cTime.strftime('%j'))/365.25) # day of year
    siteDF['DateTime'] = newTimes
    siteDF['FracYear'] = fracyr
    return siteDF


def changeTimes_change_ISOdata(siteDF):
    timeStrs = siteDF['collectDate']
    newTimes = [] ; fracyr = []
  
    for t in timeStrs:
        cTime = datetime.strptime(str(t),"%Y-%m-%d %H:%M:%S")
        newTimes.append(cTime)
        fracyr.append(float(cTime.strftime('%j'))/365.25) # day of year
    siteDF['DateTime'] = newTimes
    siteDF['FracYear'] = fracyr
    return siteDF


''' Calculate an aggregated running mean of the time series '''
def getRunningMean(daylist,tsY,tsP,year_frac,w):
    tsYm = [] ; tsPm = [] ; tsXm = [] ; dl = []
    numList = np.arange(np.min(daylist),np.max(daylist)+w,w) # 'w' day agg. ts
    
    for i in np.arange(len(numList)-1):
        new_iso = [] ; new_p = []  ; fracyr = []     
        
        for z in np.arange(len(daylist)):
            
            if (daylist[z] >= numList[i]) & (daylist[z] < numList[i+1]):
                new_iso.append(tsY[z])
                new_p.append(tsP[z])
                fracyr.append(year_frac[z])
        
        new_iso = np.array(new_iso)
        new_p = np.array(new_p)
        
        if new_p.size > 0:
            if np.nanmean(new_p) > 0:
                tsPm.append(np.nanmean(new_p)) # calc n-day total
                ts = (np.nansum(new_iso*new_p)/np.nansum(new_p)) # P weighed total
                tsYm.append(ts)   
                tsXm.append(np.max(fracyr))
                dl.append(numList[i])
    
    return tsYm,tsPm,tsXm,dl


''' Read precipitation and isotope data '''
def read_p_data(precip_data_loc, precip_filter, iso_data_loc):

    # 30min Precipitation Data
    df_P30 = pd.read_excel(precip_data_loc,index=False)
    df_P30 = changeTimes_change_Pdata(df_P30)
    
    df_P30b = df_P30.set_index('DateTime')
    
    if 'secPrecipBulk' in df_P30:
        df_P30b.loc[df_P30b['secPrecipBulk'] < precip_filter,'secPrecipBulk'] = 0
        precip_daily = df_P30b['secPrecipBulk'].resample('D').sum() # sum to total daily P
    
    else:
        df_P30b.loc[df_P30b['priPrecipBulk'] < precip_filter,'priPrecipBulk'] = 0
        precip_daily = df_P30b['priPrecipBulk'].resample('D').sum() # sum to total daily P
    
    frac_year = df_P30b['FracYear'].resample('D').mean() # average daily year fraction

    daily_P = pd.DataFrame({'Total P': precip_daily,'FracYear': frac_year})
    daily_P['Total P'].replace(0, np.nan, inplace=True)
    
    print("Successfully read precipitation amount data...")

    # Biweekly Stable Water Isotope Data
    df_iso = pd.read_excel(iso_data_loc, index=False)
    df_iso = changeTimes_change_ISOdata(df_iso)

    # Create Biweekly Precipitation Amount Timeseries to Correspond to Recorded Isotope Values
    df_iso['setDate'] = pd.to_datetime(df_iso['setDate'])  
    df_iso['collectDate'] = pd.to_datetime(df_iso['collectDate']) 

    P14 = []
    for i in np.arange(len(df_iso['setDate'])):
        subset = ((df_P30['DateTime'] > df_iso['setDate'].iloc[i]) 
                & (df_P30['DateTime'] <= df_iso['collectDate'].iloc[i]))

        df_sub = df_P30.loc[subset]
        if 'secPrecipBulk' in df_P30:
            P14.append(np.nansum(df_sub['secPrecipBulk'].values))
        else:
            P14.append(np.nansum(df_sub['priPrecipBulk'].values))

    df_iso['Total P'] = P14
    del P14, i, subset
    
    print("Successfully read precipitation isotope data...")
    
if __name__ == '__main__':
    read_p_data(precip_data_loc, precip_filter, iso_data_loc)