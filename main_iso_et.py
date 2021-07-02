# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 14:55:36 2021

@author: libon
"""


import matplotlib.pyplot as plt
from functools import reduce
import datetime
import numpy as np
import pandas as pd
import copy
import os


'''This script read in the keelinged output files and 
    generate the flags, filled values is -9999.
    Users should filter these values on the istope dataset and the flag dataset before
    qc the dataset
    
    xxxx_qc : xxxx is the four letter abbreviation of NEON sites, 
                   0: the data point has the recommended quality,
                   1: no calibration isotope
                   2: this point has number of keeling points >= 5
                   3: the keeling r2 is smaller than 0.9
                   4: this point is beyond the IQR (IQR is based on rsq >= 0.9, # keeling >= 5)
'''




def getTimeSeriesPlot(site,name,df):
    cols = df.columns[1:]
    fig, ax = plt.subplots(3,figsize= (70, 60))
    for i in np.arange(len(cols)):
        x = np.arange(len(df))
        ax[i].set_xlim([0,len(df)])
        ax[i].set(xticks=x, xticklabels = df['date'])
        ax[i].scatter(x, df[cols[i]], c = 'r',
                      s = 300, alpha = 0.8)
        for n, label in enumerate(ax[i].xaxis.get_ticklabels()):
            if n % 200 != 0:
                label.set_visible(False)
        ax[i].tick_params("both", labelsize = 70)
        ax[i].set_ylabel(cols[i],fontsize = 80)
    fig.suptitle(name + " (" + site + ')', x = 0.5, y = 0.9, fontsize= 80)

       
def getCleanedIsotope(rawIsoPath, whichIsotope,
                      outputIsoPath = None,
                      isTimeSeries = False): 
                                                         
    '''     rawIsoPath  : the absolute path of istope .csv files
          whichIsotope  : the isotope you need please either specify H2, O18 or C13 
         outputIsoPath  : where to save the processed dataset
          isTimeSeries  : if the user need the time series plot for the isotope 
          
    '''
    ##read in the raw isotope dataset 
    assert whichIsotope in ['C13', 'O18', 'H2'], "Please name whichIsotope to be one of [C13, O18, H2]"
    rawIsotope = pd.read_csv(rawIsoPath, index_col=0) 
    # isotopeSite = list(rawIsotope['site'].unique())
    listOfDataDF = []
    listOfFlagDf = []
    isotopeSites = ["BONA","CLBJ","CPER","GUAN","HARV","KONZ",
                    "NIWO","ONAQ","ORNL","OSBS","PUUM","SCBI",
                    "SJER","SRER","TALL","TOOL","UNDE","WOOD",
                    "WREF","YELL","ABBY","BARR","BART","BLAN",
                    "DELA","DSNY","GRSM","HEAL","JERC","JORN",
                    "KONA","LAJA","LENO","MLBS","MOAB","NOGP",
                    "OAES","RMNP","SERC","SOAP","STEI","STER",
                    "TEAK","TREE","UKFS","DCFS","DEJU"]
    isotopeSites.sort()
    print(isotopeSites)
    print(rawIsotope.duplicated().shape[0], rawIsotope.shape[0])
    for site in isotopeSites:
        
        iso = rawIsotope[rawIsotope['site'] == site]
        if iso.empty:
            raw = pd.DataFrame(pd.date_range('2017-01-01', '2020-01-01'), columns=['date'])
            raw['date'] = raw['date'].dt.date
            raw[site] = raw.shape[0]*[np.nan]
            listOfDataDF.append(raw[['date', site]])
            listOfFlagDf.append(raw[['date', site ]])
                                   
                                    
            continue
          
           
        ########raw isotope dataset ########
        isotopeFrame = copy.deepcopy(iso) 
        isotopeFrame['date'] = pd.to_datetime(isotopeFrame['date'], format = '%Y-%m-%d %H:%M:%S')
        isotopeFrame['date'] = isotopeFrame['date'].dt.date
        #######################################
        isoRawData = copy.deepcopy(isotopeFrame)
        isoRawData.rename(columns = {'flux':site}, inplace=True)        
        listOfDataDF.append(isoRawData[['date', site]])
        
        #######Generating isotope flags######
        isoFlagData = copy.deepcopy(isotopeFrame)
        ##number of points to generate keeling plot nptsReg >= 5 is 0  1 otherwise
        isoFlagData[site] = np.where(isoFlagData['nptsReg'] >= 5, 0, 2) 
        
        ##if number of points is good then assign rsq < 0.9 to be 2 
        isoFlagData.loc[(isoFlagData[site] == 0) & (isoFlagData['rsq'] < 0.9),
                        site] = 3
        

        isoFlagData.loc[(isoFlagData[site] != 0),'flux'] = np.nan   
                     
        p75, p25 = np.nanpercentile(isoFlagData['flux'], [75, 25]) 
        IQR = p75 - p25
        
        isoFlagData.loc[(isoFlagData[site] == 0) , site] = isoFlagData.loc[(isoFlagData[site] == 0) , 'flux'].apply(lambda x: 0 if (x >= p25 - 1.5*IQR and x <= p75 + 1.5*IQR) else 4)
        
        if isTimeSeries:
            rawSeries = isoRawData[['date', site]].copy()
            rawSeries.rename(columns = {site:'raw'}, inplace = True)

            afterqcFlag = isoFlagData[['date', 'flux']].copy()
            afterAllFlag = isoFlagData.loc[isoFlagData[site] == 0, ['date','flux']].copy()
            afterAllFlag.rename(columns = {'flux':'after flag'}, inplace = True)
            
            merged_ = reduce(lambda x, y: pd.merge(x, y, on = 'date', how = 'outer'), 
                              [rawSeries, afterqcFlag, afterAllFlag])  
            getTimeSeriesPlot(site,whichIsotope,merged_)

    ##Overall time series starts and ends        
        listOfFlagDf.append(isoFlagData[['date', site]])
                                             

      ################flag datafrane  
    flags = reduce(lambda x, y: pd.merge(x, y, on = 'date', how = 'outer'), 
                              listOfFlagDf)  
    Data =  reduce(lambda x, y: pd.merge(x, y, on = 'date', how = 'outer'), 
                                listOfDataDF)
    assert min(flags['date']) == min(Data['date']), "start time of flags and data should be the same"
    assert max(flags['date']) == max(Data['date']), "end time of flags and data should be the same"
    
    startDate = min(flags['date'])
    endDate = max(flags['date'])
    
    dts = pd.DataFrame({'date':pd.date_range(startDate, endDate)})
    dts['date'] = dts['date'].dt.date
    
    flags = pd.merge(flags, dts, on = 'date', how = 'outer')
    flags = flags.fillna(1)
    print(flags.drop_duplicates().shape[0], flags.shape[0])
    flags.drop_duplicates(subset=None, keep='first', 
                          inplace=True, ignore_index=True)
    flags.sort_values(by=['date'], ignore_index = True, inplace=True)
    
    
    ################isotope dataframe
    Data = Data.merge(dts, on = 'date', how = 'outer')
    Data = Data.fillna(-9999)
    Data.sort_values(by=['date'], ignore_index = True, inplace=True)
    print(Data.drop_duplicates().shape[0], Data.shape[0])
    Data.drop_duplicates(subset=None, keep='first', 
                          inplace=True, ignore_index=True)
    Data.sort_values(by=['date'], ignore_index = True, inplace=True)
    
    if outputIsoPath is not None:
        Data.to_csv(outputIsoPath + whichIsotope + '_data.csv', index = False)
        flags.to_csv(outputIsoPath + whichIsotope + '_qc.csv', index = False)
    else:
        return Data, flags
        # return None
       

if __name__ == '__main__':
  C13_iso = getCleanedIsotope('D:/NEON_daily_iso/keeling isotopes/et_C13_iso.csv'
                          ,'C13',
                          'D:/NEON_daily_iso/Output/')
                         
  
  H2_iso = getCleanedIsotope('D:/NEON_daily_iso/keeling isotopes/et_H2_iso.csv'
                           ,'H2', 
                          'D:/NEON_daily_iso/Output/')
  
  O18_iso = getCleanedIsotope('D:/NEON_daily_iso/keeling isotopes/et_O18_iso.csv'
                           ,'O18', 
                          'D:/NEON_daily_iso/Output/')
                         
    
    



    
    
    
    