#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPIAS analysis
"""

## Set experiment details
Animal = 'TestGPIAS_01'
Exp = 'TestGPIASn01-20170405-GPIAS'
RecFolder = 'KwikFiles/2017-04-05_15-44-55'
ExpFile = '20170405152613-Test01-GPIAS.hdf5'

PiezoCh = [8]
TTLCh = 3
GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [70, 300]     # frequency for filter
FilterOrder = 3       # butter order
AnalogTTLs = True
Override = {}

import DataAnalysis, Hdf5F

DataFolder = Animal + '/' + Exp + '/' + RecFolder
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'
AnalysisKey = Exp + '/' + RecFolder.split('/')[-1]

Data = Hdf5F.LoadOEKwik(DataFolder, AnalogTTLs, 'uV')[0]
DataInfo = Hdf5F.LoadDict('/DataInfo', Animal + '/' + Exp + '/' + ExpFile)
Proc = Hdf5F.GetProc(Data, 'OE')
Rate = Data[Proc]['info']['0']['sample_rate']

DataInfo['PiezoCh'] = PiezoCh
DataInfo['TTLCh'] = TTLCh
for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
    DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, Animal + '/' + Exp + '/' + ExpFile)


DataAnalysis.GPIASAnalysis(Data[Proc]['data'], DataInfo, Rate, AnalysisFile, AnalysisKey, 
                           GPIASTimeBeforeTTL, GPIASTimeAfterTTL, FilterFreq, 
                           FilterOrder, AnalogTTLs, Override)