#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 18:24:31 2021

@author: catiefinkenbiner
"""

import numpy as np
import pandas as pd
from datetime import datetime
import scipy as sp
from scipy import stats
from scipy import optimize


''' Read in date time and calculate fractional days for P and iso data'''
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


''' Calculate statistics '''
def calcSiteStats(site_stats, site, lamda, agglev, P, H, O):
    Pmu = np.mean(P)            # means
    Hmu = np.mean(H)
    Omu = np.mean(O)
    
    Psigma = np.std(P)           # standard deviations
    Hsigma = np.std(H)
    Osigma = np.std(O)
    
    PH_pearson = sp.stats.pearsonr(P,H)[0]     # Pearson correlation coeff.
    PO_pearson = sp.stats.pearsonr(P,O)[0]  
    HO_pearson = sp.stats.pearsonr(H,O)[0] 
    
    site_stats.append([site, agglev, lamda, Pmu, Hmu, Omu,
                       Psigma, Hsigma, Osigma,
                       PH_pearson, PO_pearson, HO_pearson])
    

''' Calculate sinisoidal function + parameters '''
def sine_func(t, amplitude, phase, offset):
    return amplitude * np.sin(2*np.pi*t - phase) + offset  # assume phase is decimal (year/year)

    
def sine_params(name, agglevel, frac_yr, P, H, O):  
    param_bounds=([0, -np.pi, -np.inf], [np.inf, np.pi, np.inf])               
    
    # d2H - find parameters for sine function based on fractional year 
    params1, params_covariance = optimize.curve_fit(sine_func, frac_yr, H, p0 = [np.std(H)*(2*np.sqrt(2)), np.pi/2,np.mean(H)], bounds = param_bounds)
    
    # d18O - find parameters for sine function based on fractional year 
    params2, params_covariance = optimize.curve_fit(sine_func, frac_yr, O, p0 = [np.std(O)*(2*np.sqrt(2)), np.pi/2,np.mean(O)], bounds = param_bounds)
      
    return params1, params2


''' Use prediction to create conditional copula generated values '''
def conditional_copula_ts(tsP, xday_stats, H_scale, O_scale):
    # Observed P statistics
    Pmu = np.mean(tsP)
    Psig = np.std(tsP)
    
    # Correlation Coefficients
    xday_rho1 = xday_stats[2,0] ; xday_rho2 = xday_stats[2,1] ; xday_rho3 = xday_stats[2,2] 
    
    # Calculate Sigma_bar
    sigma_bar = [[1 - xday_rho1*xday_rho1, xday_rho3 - xday_rho1*xday_rho2],
                 [xday_rho3 - xday_rho1*xday_rho2, 1 - xday_rho2*xday_rho2]]
    
    series = []
    for i in np.arange(len(tsP)):
        if tsP[i] > 0:
            a = (tsP[i] - Pmu) / Psig        
            H2mu_bar = xday_rho1 * a
            O18mu_bar = xday_rho2 * a
            mu_bar = np.array([H2mu_bar,O18mu_bar])

            X = sp.stats.multivariate_normal.rvs(mean= mu_bar,cov= sigma_bar) 
            X2 = sp.stats.norm.cdf(X[0], loc=0, scale=1) 
            X3 = sp.stats.norm.cdf(X[1], loc=0, scale=1)

            index = int(np.floor(X2 * len(H_scale)))
            newH = H_scale[index]

            index = int(np.floor(X3 * len(O_scale)))
            newO = O_scale[index]

            series.append([newH, newO]) 
        else:
            series.append([np.nan, np.nan])
    return series


''' Statistical Downscaling Method (Finkenbiner et al. 2021) '''
def downscaling_to_daily(sitename, daily_P, df_iso, out_dir_loc, residual_corr=True):
    print('Downscaling '+ sitename +' ...')
    
    '''
    Step 1 - Remove seasonal time series component (Section 2.b.1 in text)
    '''
    
    # Define Sinisoidal Functions describing seasonality
    df_iso = df_iso.sort_values('DateTime')
    df_iso = df_iso.dropna(subset=['Total P'])
    tsX = df_iso['FracYear'].values 
    tsP = df_iso['Total P'].values
    tsO = df_iso['d18OWater'].values
    tsH = df_iso['d2HWater'].values   

    dayslist = []
    for dt in np.arange(len(df_iso['DateTime'])):
        dayslist.append((df_iso['DateTime'].iloc[dt] - df_iso['DateTime'].iloc[0]).days)
    dayslist = np.array(dayslist)

    # lambda = precipitation frequency (see Eq. 4)
    p_events = df_iso[df_iso['Total P'].notna()]
    lamda = len(p_events['Total P'])/((daily_P.index.max() - daily_P.index.min()).days) 

    params1, params2 = sine_params(sitename, 14, tsX, tsP, tsH, tsO) # 14 = biweekly sample, sample frequency
    tsY_sine_wave = sine_func(tsX, params1[0], params1[1], params1[2]) # solve for amplitude, phase, offset
    adj_2H = np.array((tsH - tsY_sine_wave))    # remove seasonality from time series 

    tsY_sine_wave = sine_func(tsX, params2[0], params2[1], params2[2])     
    adj_18O = np.array((tsO - tsY_sine_wave)) 
    
    # Get biweekly site stats of stochastic component
    biweekly_stats = np.array([[np.mean(tsP), np.mean(adj_2H), np.mean(adj_18O)],
                     [np.std(tsP), np.std(adj_2H), np.std(adj_18O)],
                     [sp.stats.pearsonr(tsP,adj_2H)[0], sp.stats.pearsonr(tsP,adj_18O)[0], sp.stats.pearsonr(adj_2H,adj_18O)[0]]])
    
    '''
    Step 2 - Predict daily statistics from biweekly time series (Section 2.b.2 in text) Now we will need to aggregate the 
    stochastic biweekly time series - i.e. calculated weighed running means at biweekly (14-day), 28-day, 42-day, 
    56-day and 84-day intervals    
    '''

    # Define time series statistics 
    site_stats = [[sitename, 14, lamda,
                   biweekly_stats[0,0], biweekly_stats[0,1], biweekly_stats[0,2],
                   biweekly_stats[1,0], biweekly_stats[1,1], biweekly_stats[1,2],
                   biweekly_stats[2,0], biweekly_stats[2,1], biweekly_stats[2,2]]]
    site_stats_check = len(site_stats)

    for n in np.arange(28,85,14):
        xday_Hb, xday_Pb, xday_Xb, days = getRunningMean(np.array(dayslist), np.array(tsH), np.array(tsP), tsX, n)
        xday_Ob, xday_Pb, xday_Xb, days = getRunningMean(np.array(dayslist), np.array(tsO), np.array(tsP), tsX, n)
        
        if len(xday_Xb) > 4:

            xday_Xb = np.array(xday_Xb)
            xday_Pb = np.array(xday_Pb)
            xday_Hb = np.array(xday_Hb)
            xday_Ob = np.array(xday_Ob)

            params1a, params2a = sine_params(sitename, n, xday_Xb, xday_Pb, xday_Hb, xday_Ob)

            tsY_sine_wave = sine_func(xday_Xb, params1a[0],params1a[1],params1a[2])
            adj_2Hb = np.array((xday_Hb - tsY_sine_wave))                           

            tsY_sine_wave = sine_func(xday_Xb, params2a[0],params2a[1],params2a[2])     
            adj_18Ob = np.array((xday_Ob - tsY_sine_wave))
        
            if len(xday_Pb)>2 and len(adj_2Hb)>2 and len(adj_18Ob)>2:
                calcSiteStats(site_stats, sitename, lamda, n, xday_Pb, adj_2Hb, adj_18Ob)
            
    if len(site_stats) == site_stats_check:
        print("Site {} does not contain sufficient data.".format(sitename))
        
    else:        
        ### the stats are labeled with 'B' here because they are of the stochastic component - not the original time series
        Site_Stats = pd.DataFrame(site_stats, columns = ['site','agglev','lambda','PmuB','HmuB','OmuB','PsigB',
                                                         'HsigB','OsigB','PHpB','POpB','HOpB'])
    
        ### Now we apply Eq. 4 from Finkenbiner et al. 2021
        xaxis = np.array(Site_Stats['agglev'].values)
        yaxis1 = np.array(Site_Stats['HsigB'].values)
        yaxis2 = np.array(Site_Stats['OsigB'].values)

        def eq4(x, a, b):
            return b/(x*lamda)**a

        bounds = [[0.2,yaxis1[0]], [0.5,np.inf]]
        p1, p2 = optimize.curve_fit(eq4, xaxis, yaxis1, p0 = [0.3, yaxis1[0]], bounds=bounds)
        Hi = float(p1[1]) # estimated daily parameter

        bounds = [[0.2,yaxis2[0]], [0.5,np.inf]]
        p1, p2 = optimize.curve_fit(eq4, xaxis, yaxis2,  p0 = [0.3, yaxis2[0]], bounds=bounds)
        Oi = float(p1[1]) # estimated daily parameter
    
        
        '''
        Step 3 - Generate daily time series with estimated statistcs (Section 2.b.2 in text)
        Step 4 - Add in seasonal time series component (Section 2.b.2 in text)
        '''
    
        H_scale = np.sort(np.array(adj_2H) * Hi / Site_Stats['HsigB'].iloc[0])
        O_scale = np.sort(np.array(adj_18O) * Oi / Site_Stats['OsigB'].iloc[0])

        copula_stats = np.matrix([[0, 0, 0], [np.std(tsP), Hi, Oi],
                                    [Site_Stats['PHpB'].iloc[0], 
                                     Site_Stats['POpB'].iloc[0],
                                     Site_Stats['HOpB'].iloc[0]]])    

        new_ts = conditional_copula_ts(daily_P['Total P'], copula_stats, H_scale, O_scale)
        y = np.array([np.array(xi) for xi in new_ts])   

        # Add back in n-day sine function here:
        tsH_daily = y[:,0] + sine_func(daily_P['FracYear'], params1[0], params1[1], params1[2])    
        tsO_daily = y[:,1] + sine_func(daily_P['FracYear'], params2[0], params2[1], params2[2])
    
        daily_P['d2H'] = tsH_daily
        daily_P['d18O'] = tsO_daily 

        
        ''' Residual Correction '''
        if residual_corr == True:
    
            # goes through each biweekly period in the isotope data
            for n in np.arange(len(df_iso['setDate'])):
                start_date = df_iso['setDate'].iloc[n]
                end_date = df_iso['collectDate'].iloc[n]
            
                # find observed H and O
                obs_d2Hw = df_iso['d2HWater'].iloc[n]
                obs_d18Ow = df_iso['d18OWater'].iloc[n]
            
                # find corresponding data in synthetic time series
                mask = (daily_P.index >= start_date) & (daily_P.index <= end_date)
                df_n = daily_P.loc[mask]
                
                syn_d2Hw = np.nansum(df_n['d2H'].values * df_n['Total P'].values) / np.nansum(df_n['Total P'].values)
                syn_d18Ow = np.nansum(df_n['d18O'].values * df_n['Total P'].values) / np.nansum(df_n['Total P'].values)
            
                # find difference between observed and synthetic
                diff_d2Hw = obs_d2Hw - syn_d2Hw
                diff_d18Ow = obs_d18Ow - syn_d18Ow
                
                daily_P['d2H'].loc[mask] = daily_P['d2H'].loc[mask] + diff_d2Hw
                daily_P['d18O'].loc[mask] = daily_P['d18O'].loc[mask] + diff_d18Ow        
                
            print('Done.')
        
        else:
            print('Done.')
    
    return daily_P
    

''' Read .xlsx precipitation and isotope data '''
def read_p_data(sitename, precip_data_loc, precip_filter, iso_data_loc, out_dir_loc, residual_corr):

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
    
    print("Successfully read "+ sitename+ " precipitation amount data...")

    # Biweekly Stable Water Isotope Data
    df_iso = pd.read_excel(iso_data_loc, index=False)
    df_iso = changeTimes_change_ISOdata(df_iso)

    ## Create Biweekly Precipitation Amount Timeseries to Correspond to Recorded Isotope Values
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
    
    print("Successfully read "+ sitename+ " precipitation isotope data...")
    
    
    ''' Statistical Downscaling Method (Finkenbiner et al. 2021) '''
    return downscaling_to_daily(sitename, daily_P, df_iso, out_dir_loc)
    

if __name__ == '__main__':
    read_p_data(sitename, precip_data_loc, precip_filter, iso_data_loc, out_dir_loc, residual_corr)