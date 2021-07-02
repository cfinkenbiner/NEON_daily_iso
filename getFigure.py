# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 08:05:01 2021

@author: libon
"""

import numpy as np
import pandas as pd
import os 


#######This script used to generate statistic as well as validate the quality flags 
dataPath = 'D:/NEON_daily_iso/keeling isotopes/'

carbon = pd.read_csv(dataPath + 'et_C13_iso.csv', index_col=0)
hydrogen = pd.read_csv(dataPath + 'et_H2_iso.csv', index_col=0)
oxygen = pd.read_csv(dataPath + 'et_O18_iso.csv', index_col=0)

newdfPath = 'D:/NEON_daily_iso/Output/'

newcarb = pd.read_csv(newdfPath + 'C13_data.csv') 
newcarbFlag = pd.read_csv(newdfPath + 'C13_qc.csv') 

newO = pd.read_csv(newdfPath + 'O18_data.csv') 
newOFlag = pd.read_csv(newdfPath + 'O18_qc.csv') 

newH = pd.read_csv(newdfPath + 'H2_data.csv') 
newHFlag = pd.read_csv(newdfPath + 'H2_qc.csv')



neonList = ["BONA","CLBJ","CPER","GUAN","HARV","KONZ",
            "NIWO","ONAQ","ORNL","OSBS","PUUM","SCBI",
            "SJER","SRER","TALL","TOOL","UNDE","WOOD",
            "WREF","YELL","ABBY","BARR","BART","BLAN",
            "DELA","DSNY","GRSM","HEAL","JERC","JORN",
            "KONA","LAJA","LENO","MLBS","MOAB","NOGP",
            "OAES","RMNP","SERC","SOAP","STEI","STER",
            "TEAK","TREE","UKFS","DCFS","DEJU"]
neonList.sort()




def getProcessed(df):
    dff = df.copy()
    dff = dff[(dff['rsq'] >= 0.9)&(dff['nptsReg'] >=5)]
    p75, p25 = np.nanpercentile(dff['flux'], [75, 25]) 
    IQR = p75 - p25
    
    return dff[(dff['flux'] >= p25 - 1.5*IQR) & (dff['flux'] <= p75 + 1.5*IQR)]

def getOriginIso(df,dfFlag,siteName):
    cData = df[['date',siteName]]
    
    cFlag = dfFlag[['date',siteName]]
    cMerge = pd.merge(cData, cFlag, on = 'date')
    cMerge.rename(columns = {siteName + '_x': siteName,
                             siteName + '_y': siteName + '_qc'} ,inplace = True)
    cMerge = cMerge[cMerge[siteName + '_qc'] == 0]
    cMerge.sort_values(by = ['date'],
                       ignore_index = True,inplace = True, )
    return cMerge



colNames= ['site', 
           'O18_ET_mean','O18_ET_std','O18_slope_mean',
           'O18_slope_std','O18_goodpoints', 'O18_frac_good_overall',
           'H2_ET_mean','H2_ET_std','H2_slope_mean',
           'H2_slope_std','H2_goodpoints','H2_frac_good_overall',
           'C13_mean', 'C13_std', 'C13_slope_mean',
           'C13_slope_std', 'C13_goodpoints', 'C13_frac_good_overall']

dfSite = pd.DataFrame(columns= colNames)

for i in np.arange(len(neonList)):
    print("Working on site {}".format(neonList[i]))

    isoO18 = oxygen[oxygen['site'] == neonList[i]]
    isoH2 = hydrogen[hydrogen['site'] == neonList[i]]
    isoC13 = carbon[carbon['site'] == neonList[i]]
    inSiteList = [neonList[i]]
    
    fH2 = getOriginIso(newH,newHFlag, neonList[i])                                                  
    fC13 = getOriginIso(newcarb,newcarbFlag, neonList[i])    
    fO18 = getOriginIso(newO,newOFlag, neonList[i])   
    
    if isoO18.empty:
        print("The original isoO18 is empty")
        print("is the flagged data empty ? : {}".format(fO18.empty))
        inSiteList.extend([np.nan]*6)
    else:
        isoO18_p = getProcessed(isoO18)
        assert isoO18_p.shape[0] == fO18.shape[0]
        
        isoO18_p.sort_values(by=['date'], 
                             ignore_index = True, inplace=True)
        O18TimeMean = isoO18_p['flux'].mean() 
        O18TimeStd = isoO18_p['flux'].std() 
        O18SlopeMean = isoO18_p['slopeStd'].mean() 
        O18SlopeStd = isoO18_p['slopeStd'].std() 
        
        assert O18TimeMean == fO18[neonList[i]].mean()
        assert O18TimeStd == fO18[neonList[i]].std()    
                                                       
        inSiteList.extend([O18TimeMean, O18TimeStd, O18SlopeMean,
                           O18SlopeStd, isoO18_p.shape[0], isoO18_p.shape[0]/newO.shape[0]])

    if isoH2.empty:
        inSiteList.extend([np.nan]*6)
        print("The original isoH2 is empty")
        print("is the flagged data empty ? : {}".format(fH2.empty))
    else:
        isoH2_p = getProcessed(isoH2)
        isoH2_p.sort_values(by=['date'], 
                             ignore_index = True, inplace=True)
        assert isoH2_p.shape[0] == fH2.shape[0]
                                                                     
        
        # print(isoH2_p.shape[0])
        H2TimeMean = isoH2_p['flux'].mean() 
        H2TimeStd = isoH2_p['flux'].std() 
        H2SlopeMean = isoH2_p['slopeStd'].mean() 
        H2SlopeStd = isoH2_p['slopeStd'].std() 
        
        assert H2TimeMean == fH2[neonList[i]].mean()
        assert H2TimeStd == fH2[neonList[i]].std()  
        
        inSiteList.extend([H2TimeMean, H2TimeStd, H2SlopeMean,
                           H2SlopeStd, isoH2_p.shape[0], isoH2_p.shape[0]/newH.shape[0]])
   
    if isoC13.empty:
        inSiteList.extend([np.nan]*6)

    else:
        isoC13_p = getProcessed(isoC13)
        isoC13_p.sort_values(by=['date'], 
                             ignore_index = True, inplace=True)
        
        assert isoC13_p.shape[0] == fC13.shape[0]
      
        
        C13TimeMean = isoC13_p['flux'].mean() 
        C13TimeStd = isoC13_p['flux'].std() 
        
        assert C13TimeMean == fC13[neonList[i]].mean()
        assert C13TimeStd == fC13[neonList[i]].std()
        C13SlopeMean = isoC13_p['slopeStd'].mean()
        C13SlopeStd = isoC13_p['slopeStd'].std() 
        inSiteList.extend([C13TimeMean, C13TimeStd, C13SlopeMean,
                           C13SlopeStd, isoC13_p.shape[0], isoC13_p.shape[0]/newcarb.shape[0]])
    dfSite.loc[i,:] = inSiteList
    print('---------------')


dfSite.to_csv('D:/NEON_daily_iso/Figure Stats/ET_iso_stats.csv',
              index = False)













