#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPIAS analysis
"""

## Set experiment details
Animal = 'GPIAZon'
Exp = 'GPIAZon_NaCl'
RecFolder = '2017-04-11_14-18-12_GPIAZon_NaCln01'
ExpFile = '20170411141734-GPIAZon_NaCln01-GPIAS.hdf5'

PiezoCh = [8]
TTLCh = 3
GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [100, 300]     # frequency for filter
FilterOrder = 3       # butter order
AnalogTTLs = True
Override = {}

import DataAnalysis, Hdf5F
import numpy as np

DataFolder = Animal + '/' + Exp + '/' + RecFolder
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'
AnalysisKey = Exp + '/' + RecFolder.split('/')[-1]

Data = Hdf5F.LoadOEKwik(DataFolder, AnalogTTLs, 'Bits')[0]
DataInfo = Hdf5F.LoadDict('/DataInfo', Animal + '/' + Exp + '/' + ExpFile)
Proc = Hdf5F.GetProc(Data, 'OE')
Rate = Data[Proc]['info']['0']['sample_rate']

for Rec in Data[Proc]['data'].keys():
    BitVolts = 10000/(2**16)
    Data[Proc]['data'][Rec] = Data[Proc]['data'][Rec] * BitVolts

#DataInfo['PiezoCh'] = PiezoCh
#DataInfo['TTLCh'] = TTLCh
for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
    DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, Animal + '/' + Exp + '/' + ExpFile)

# Fix stupid breaks in recs
def CheckGPIASRecs(Data, SizeLimits):
    ToCheck = [Rec for Rec in Data.keys() 
                   if len(Data[Rec])<min(SizeLimits) or len(Data[Rec])>max(SizeLimits)]
    
    if ToCheck:
        Params = {'backend': 'Qt5Agg'}
        from matplotlib import rcParams; rcParams.update(Params)
        import matplotlib.pyplot as plt
        
        for Rec in ToCheck:
            print('Showing rec ' + Rec)
            plt.plot(Data[Rec])
            plt.show()
        
        return(ToCheck)
    else:
        print('All recs within expected size.')
        return(None)

CheckGPIASRecs(Data[Proc]['data'], [30000, 100000])

# for np.delete(DataInfo['FreqOrder'], 60, 0)
#DataInfo['FreqOrder'] = np.delete(DataInfo['FreqOrder'], [59,60], 0)
#for Key in Data[Proc].keys(): del(Data[Proc][Key]['59'])

DataAnalysis.GPIASAnalysis(Data[Proc]['data'], DataInfo, Rate, AnalysisFile, AnalysisKey, 
                           GPIASTimeBeforeTTL, GPIASTimeAfterTTL, FilterFreq, 
                           FilterOrder, AnalogTTLs, Override)