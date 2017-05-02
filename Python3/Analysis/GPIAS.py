#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPIAS analysis
"""
#%% Import

import DataAnalysis, Hdf5F
import os
from glob import glob

#%% Batch
Animal = 'GPIAZon'
Exp = 'GPIAZon_SSal'
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'

GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [100, 300]     # frequency for filter
FilterOrder = 3       # butter order
PiezoCh = [8]
TTLCh = 1

Paths = glob(Animal + '/' + Exp + '/2017-*'); Paths.sort()
Files = glob(Animal + '/' + Exp + '/201704*'); Files.sort()

for _ in [5, 4, 0] : del(Paths[_], Files[_])


for Ind, DataPath in enumerate(Paths):
    RecFolder = DataPath.split('/')[-1]
    AnalysisKey = Exp + '/' + RecFolder.split('/')[-1]
    FigName = '/'.join([Animal, Exp, 'Figs', RecFolder+'-GPIAS'])
    FigPath = '/'.join(FigName.split('/')[:-1])
    os.makedirs(FigPath, exist_ok=True)
    
    Data = Hdf5F.LoadOEKwik(DataPath, AnalogTTLs=True, Unit='Bits')[0]
    DataInfo = Hdf5F.LoadDict('/DataInfo', Files[Ind])
    Proc = Hdf5F.GetProc(Data, 'OE')
    Rate = Data[Proc]['info']['0']['sample_rate']
    
    # Test recs
#    SOAB = DataAnalysis.GPIAS.CheckGPIASRecs(Data[Proc]['data'], [60000, 100000])
#    if SOAB: 
#        print(); print(DataPath)
#        print(); print(SOAB); print()
#    else: print(); print(DataPath, 'clear.'); print()
        
    for Rec in Data[Proc]['data'].keys():
        BitVolts = 10000/(2**16)
        Data[Proc]['data'][Rec] = Data[Proc]['data'][Rec] * BitVolts
    
    DataInfo['PiezoCh'] = PiezoCh
    DataInfo['TTLCh'] = TTLCh
    
    for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
        DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, Files[Ind])
    
    DataInfo['FreqOrder'][-3:][:,1] = -2
    
    GPIAS, XValues = DataAnalysis.GPIAS.Analysis(
                         Data[Proc]['data'], DataInfo, Rate, AnalysisFile, 
                         AnalysisKey, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, 
                         FilterFreq, FilterOrder, Return=True)
    
    
    DataAnalysis.Plot.GPIAS(GPIAS, XValues, DataInfo['SoundLoudPulseDur'], 
                            FigName, Save=True, Visible=True)


#%% Individual

Animal = 'GPIAZon'
Exp = 'GPIAZon_SSal'
RecFolder = '2017-04-15_14-28-03_GPIAZon_SSaln04'
ExpFile = '20170415142722-GPIAZon_SSaln04-GPIAS.hdf5'

GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [100, 300]     # frequency for filter
FilterOrder = 3       # butter order
PiezoCh = [8]
TTLCh = 1

DataPath = Animal + '/' + Exp + '/' + RecFolder
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'
AnalysisKey = Exp + '/' + RecFolder.split('/')[-1]
FigName = '/'.join([Animal, Exp, 'Figs', RecFolder+'-GPIAS'])
FigPath = '/'.join(FigName.split('/')[:-1])
os.makedirs(FigPath, exist_ok=True)

Data = Hdf5F.LoadOEKwik(DataPath, AnalogTTLs=True, Unit='Bits')[0]
DataInfo = Hdf5F.LoadDict('/DataInfo', Animal + '/' + Exp + '/' + ExpFile)
Proc = Hdf5F.GetProc(Data, 'OE')
Rate = Data[Proc]['info']['0']['sample_rate']

for Rec in Data[Proc]['data'].keys():
    BitVolts = 10000/(2**16)
    Data[Proc]['data'][Rec] = Data[Proc]['data'][Rec] * BitVolts

DataInfo['PiezoCh'] = PiezoCh
DataInfo['TTLCh'] = TTLCh

for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
    DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, Animal + '/' + Exp + '/' + ExpFile)

DataInfo['FreqOrder'][-3:][:,1] = -2

## Fix stupid breaks in recs
#import numpy as np
#SOAB = DataAnalysis.GPIAS.CheckGPIASRecs(Data[Proc]['data'], [65000, 100000]); SOAB.sort()
#Params = {'backend': 'TkAgg'}
#from matplotlib import rcParams; rcParams.update(Params)
#from matplotlib import pyplot as plt
#R = {}
#for r in ['1', '3', '4', '5'] + \
#         [str(_) for _ in range(10, 23)] + \
#         [str(_) for _ in range(32, 39)] + \
#         [str(_) for _ in range(47, 58)] + \
#         [str(_) for _ in range(68, 73)] + \
#         [str(_) for _ in range(78, 82)] + \
#         [str(_) for _ in range(85, 90)]: 
#    print(r)
#    R[r] = DataAnalysis.PSD(Data[Proc]['data'][r][:,2], Rate)
#    plt.figure(); plt.semilogy(R[r][0], R[r][1])
#    plt.figure(); plt.plot(Data[Proc]['data'][r][:,2])
#    plt.show()

# 20170411141734-GPIAZon_NaCln01-GPIAS.hdf5
#DataInfo['FreqOrder'] = np.delete(DataInfo['FreqOrder'], [60], 0)
#for Key in Data[Proc].keys(): del(Data[Proc][Key]['59'])

# 20170411160111-GPIAZon_NaCln03-GPIAS.hdf5
#for Key in Data[Proc].keys(): del(Data[Proc][Key]['110'])

# 20170411165239-GPIAZon_NaCln04-GPIAS.hdf5
#for Key in Data[Proc].keys(): del(Data[Proc][Key]['15'])

# 20170413124748-GPIAZon_SSaln01-GPIAS.hdf5
#SOAB = ['51', '46', '36', '54', '38', '58', '23', '63', '70', '48', '34']
#ToDel = [_ for _ in range(35,40)] + \
#        [_ for _ in range(42,45)] + \
#        [_ for _ in range(47,55)]
#DataInfo['FreqOrder'] = np.delete(DataInfo['FreqOrder'], ToDel, 0)
#for Key in Data[Proc].keys(): 
#    for r in SOAB: del(Data[Proc][Key][r])

# 20170415142722-GPIAZon_SSaln04-GPIAS.hdf5
#for Key in Data[Proc].keys(): 
#    for r in SOAB: del(Data[Proc][Key][r])

## Run Analysis
GPIAS, XValues = DataAnalysis.GPIAS.Analysis(
                     Data[Proc]['data'], DataInfo, Rate, AnalysisFile, 
                     AnalysisKey, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, 
                     FilterFreq, FilterOrder, Return=True)


DataAnalysis.Plot.GPIAS(GPIAS, XValues, DataInfo['SoundLoudPulseDur'], 
                        FigName, Save=True, Visible=True)