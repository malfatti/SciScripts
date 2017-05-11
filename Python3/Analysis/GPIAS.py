#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPIAS analysis
"""
#%% Import
import DataAnalysis, GPIAZon, Hdf5F
import numpy as np
from glob import glob
from os import makedirs

#%% Group
AnalysisFile = 'GPIAZon/GPIAZon-Analysis.hdf5'
Groups = ['GPIAZon_NaCl', 'GPIAZon_SSal']
ExpList = ['NaCl', 'SSal']#, 'Atr']
Exps = {'GPIAZon_NaCl': ['NaCl', 'SSal', 'Atr'],
        'GPIAZon_SSal': ['SSal', 'Atr', 'NaCl']}

Save = False; Invalid = False
DiffThr = 0.6; InvalidThr = 0.1

Index = GPIAZon.GetIndex(Groups, Exps, AnalysisFile)
Diff = GPIAZon.GetDiff(Index, Groups, ExpList, DiffThr, Invalid, InvalidThr)

MeansV, MeansFull = GPIAZon.GetMeansV(Index, Groups, ExpList)
NaCl, SSal = GPIAZon.GetIndFreqMean(Index, Diff, Groups)

Pairs, PairsFull = GPIAZon.GetPairs(MeansV, MeansFull, Groups, ExpList)
Pairs, PairsFull = GPIAZon.ClearPairs(Pairs, PairsFull)

YMax = 0.4
GPIAZon.Plot.Index_Freq_Exp_BP(MeansFull, PairsFull, ExpList, Save)
GPIAZon.Plot.Index_Exp_BP(NaCl, SSal, ExpList, YMax, Invalid, Save)


#%% Batch
Animal = 'GPIAZon'
Exp = 'GPIAZon_NaCl'
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'

GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [100, 300]     # frequency for filter
FilterOrder = 3       # butter order

Paths = glob(Animal + '/' + Exp + '/2017-*'); Paths.sort()
Files = glob(Animal + '/' + Exp + '/20170*'); Files.sort()

# NaCl
del(Paths[3], Paths[2], Paths[0]); del(Files[3], Files[2], Files[0])
# SSal
#del(Paths[-2:], Paths[5], Paths[0]); del(Files[-2:], Files[5], Files[0])

for Ind, DataPath in enumerate(Paths):
    RecFolder = DataPath.split('/')[-1]
    AnalysisKey = Exp + '/' + RecFolder.split('/')[-1]
    FigName = '/'.join([Animal, Exp, 'Figs', RecFolder+'-GPIAS'])
    FigPath = '/'.join(FigName.split('/')[:-1])
    makedirs(FigPath, exist_ok=True)
    
    Data = Hdf5F.LoadOEKwik(DataPath, AnalogTTLs=True, Unit='Bits')[0]
    DataInfo = Hdf5F.LoadDict('/DataInfo', Files[Ind])
    Proc = Hdf5F.GetProc(Data, 'OE')
    Rate = Data[Proc]['info']['0']['sample_rate']
    
    # Test TTLCh
#    plt.plot(Data[Proc]['data']['0'][:,TTLCh-1]); plt.show()
#    print(DataInfo['TTLCh'])
    
    # Test recs
#    SOAB = DataAnalysis.GPIAS.CheckGPIASRecs(Data[Proc]['data'], [65000, 100000])
#    if SOAB: 
#        print(); print(DataPath)
#        print(); print(SOAB); print()
#    else: print(); print(DataPath, 'clear.'); print()
    
    for Rec in Data[Proc]['data'].keys():
        BitVolts = 10000/(2**16)
        Data[Proc]['data'][Rec] = Data[Proc]['data'][Rec] * BitVolts
    
    DataInfo['PiezoCh'] = [int(DataInfo['PiezoCh'])]
    DataInfo['TTLCh'] = int(DataInfo['TTLCh'])
    
    for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
        DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, Files[Ind])
    
    DataInfo['FreqOrder'][-3:][:,1] = -2
    
    GPIAS, XValues = DataAnalysis.GPIAS.Analysis(
                         Data[Proc]['data'], DataInfo, Rate, AnalysisFile, 
                         AnalysisKey, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, 
                         FilterFreq, FilterOrder, Return=True)
    
    
    DataAnalysis.Plot.GPIAS(GPIAS, XValues, DataInfo['SoundLoudPulseDur'], 
                            FigName, Save=True, Visible=False)
    
    del(GPIAS, XValues)


#%% Batch broken
Animal = 'GPIAZon'
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'

GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [100, 300]     # frequency for filter
FilterOrder = 3       # butter order

Paths = [
    'GPIAZon/GPIAZon_NaCl/2017-04-11_14-18-12_GPIAZon_NaCln01',
    'GPIAZon/GPIAZon_NaCl/2017-04-11_16-01-38_GPIAZon_NaCln03',
    'GPIAZon/GPIAZon_NaCl/2017-04-11_16-53-06_GPIAZon_NaCln04',
    'GPIAZon/GPIAZon_SSal/2017-04-13_12-48-26_GPIAZon_SSaln01',
    'GPIAZon/GPIAZon_SSal/2017-04-15_14-28-03_GPIAZon_SSaln04',
    'GPIAZon/GPIAZon_SSal/2017-04-21_14-47-12_GPIAZon_SSaln04',
    'GPIAZon/GPIAZon_SSal/2017-04-24_13-54-34_GPIAZon_SSaln05'
    ]

Files = [
    'GPIAZon/GPIAZon_NaCl/20170411141734-GPIAZon_NaCln01-GPIAS.hdf5',
    'GPIAZon/GPIAZon_NaCl/20170411160111-GPIAZon_NaCln03-GPIAS.hdf5',
    'GPIAZon/GPIAZon_NaCl/20170411165239-GPIAZon_NaCln04-GPIAS.hdf5',
    'GPIAZon/GPIAZon_SSal/20170413124748-GPIAZon_SSaln01-GPIAS.hdf5',
    'GPIAZon/GPIAZon_SSal/20170415142722-GPIAZon_SSaln04-GPIAS.hdf5',
    'GPIAZon/GPIAZon_SSal/20170421144530-GPIAZon_SSaln04-GPIAS.hdf5',
    'GPIAZon/GPIAZon_SSal/20170424135406-GPIAZon_SSaln05-GPIAS.hdf5'
    ]

ToDelete = [
    [['59'], [60]],
    [['110'], []],
    [['15'], []],
    [['51', '46', '36', '54', '38', '58', '23', '63', '70', '48', '34'], [35, 36, 37, 38, 39, 42, 43, 44, 47, 48, 49, 50, 51, 52, 53, 54]],
    [['13', '14', '16', '2', '20', '33', '37', '49', '50', '54', '55', '56', '70', '71', '80', '87'], []],
    [['49'], []],
    [['53'], []]
    ]

for Ind, DataPath in enumerate(Paths):
    Exp = DataPath.split('/')[1]
    RecFolder = DataPath.split('/')[-1]
    AnalysisKey = Exp + '/' + RecFolder.split('/')[-1]
    FigName = '/'.join([Animal, Exp, 'Figs', RecFolder+'-GPIAS'])
    FigPath = '/'.join(FigName.split('/')[:-1])
    makedirs(FigPath, exist_ok=True)
    
    Data = Hdf5F.LoadOEKwik(DataPath, AnalogTTLs=True, Unit='Bits')[0]
    DataInfo = Hdf5F.LoadDict('/DataInfo', Files[Ind])
    Proc = Hdf5F.GetProc(Data, 'OE')
    Rate = Data[Proc]['info']['0']['sample_rate']
    
    for Rec in Data[Proc]['data'].keys():
        BitVolts = 10000/(2**16)
        Data[Proc]['data'][Rec] = Data[Proc]['data'][Rec] * BitVolts
    
    DataInfo['PiezoCh'] = [int(DataInfo['PiezoCh'])]
    DataInfo['TTLCh'] = int(DataInfo['TTLCh'])
    
    for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
        DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, Files[Ind])
    
    DataInfo['FreqOrder'][-3:][:,1] = -2
    
    DataInfo['FreqOrder'] = np.delete(DataInfo['FreqOrder'], ToDelete[Ind][1], 0)
    for Key in Data[Proc].keys(): 
        for r in ToDelete[Ind][0]: del(Data[Proc][Key][r])
    
    GPIAS, XValues = DataAnalysis.GPIAS.Analysis(
                         Data[Proc]['data'], DataInfo, Rate, AnalysisFile, 
                         AnalysisKey, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, 
                         FilterFreq, FilterOrder, Return=True)
    
    
    DataAnalysis.Plot.GPIAS(GPIAS, XValues, DataInfo['SoundLoudPulseDur'], 
                            FigName, Save=True, Visible=True)
    
    del(GPIAS, XValues)


#%% Individual

Animal = 'GPIAZon'
Exp = 'GPIAZon_NaCl'
RecFolder = '2017-04-11_16-53-06_GPIAZon_NaCln04'
ExpFile = '20170411165239-GPIAZon_NaCln04-GPIAS.hdf5'

GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [100, 300]     # frequency for filter
FilterOrder = 3       # butter order
PiezoCh = [8]
TTLCh = 6

DataPath = Animal + '/' + Exp + '/' + RecFolder
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'
AnalysisKey = Exp + '/' + RecFolder.split('/')[-1]
FigName = '/'.join([Animal, Exp, 'Figs', RecFolder+'-GPIAS'])
FigPath = '/'.join(FigName.split('/')[:-1])
makedirs(FigPath, exist_ok=True)

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
#SOAB = ['13', '14', '16', '2', '20', '33', '37', '49', '50', '54', '55', '56', '70', '71', '80', '87']
#for Key in Data[Proc].keys(): 
#    for r in SOAB: del(Data[Proc][Key][r])

## Run Analysis
GPIAS, XValues = DataAnalysis.GPIAS.Analysis(
                     Data[Proc]['data'], DataInfo, Rate, AnalysisFile, 
                     AnalysisKey, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, 
                     FilterFreq, FilterOrder, Return=True)


DataAnalysis.Plot.GPIAS(GPIAS, XValues, DataInfo['SoundLoudPulseDur'], 
                        FigName, Save=True, Visible=True)