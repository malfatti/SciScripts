# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 08:59:41 2016

@author: cerebro
"""

import h5py

#%% Load SoundMeasurement files
#FileName = '20160418180701-SoundMeasurement/20160418180701-SoundMeasurement.hdf5'

def SoundMeasurement(FileName):
    DataInfo = {}; SoundIntensity = {}
    with h5py.File(FileName) as h5:
        SoundRec = [0]*len(DataInfo['NoiseFrequency'])
        
        for Freq in range(len(SoundRec)):
            SoundRec[Freq] = [0]*len(DataInfo['SoundAmpF'])
            Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' + \
                  str(DataInfo['NoiseFrequency'][Freq][1])
            
            for AmpF in range(len(SoundRec[Freq])):
                aKey = str(DataInfo['SoundAmpF'][AmpF])
                SoundRec[Freq][AmpF] = list(h5['SoundRec'][Key][aKey])
        
        for Key,Val in h5['SoundRec'].attrs.items():
            DataInfo[Key] = Val
        
        for Freq in h5['SoundIntensity'].keys():
            SoundIntensity[Key] = {}
            
            for AmpF in h5['SoundIntensity'][Key].keys():
                SoundIntensity[Key][AmpF] = h5['SoundIntensity'][Key][AmpF][0]
    
    return(DataInfo, SoundRec, SoundIntensity)

