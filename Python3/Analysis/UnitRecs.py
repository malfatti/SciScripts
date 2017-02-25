#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 22:25:17 2017

@author: malfatti
"""
#%% Setting
from glob import glob

import numpy as np
import os

import DataAnalysis
import IntanBin

#%% Old Arch recs w/o TTLs
Here = os.getcwd(); Path = 'SepCh/'
ClusterPath = Here + '/' + Path

RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
A16 = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

ChannelMap = DataAnalysis.RemapChannels(A16['Tip'], A16['Head'], RHAHeadstage)

Rate = np.array([25000]); AnalogTTLs = True; Override={}

## Arch3_01
Files = [['Intan/Arch45S_153244.int', [1652313]]]
### Arch3_02
#Files = [['IntanSound/Arch3n2-43S_184905.int', [154934, 4121425]],
#         ['IntanSound/Arch3n2-35S_175902.int', [1702747]],
##         'IntanSound/Arch3n2-40S_184417.int',
#         ['IntanSound/Arch3n2-45S_191253.int', [1283124]]]
## Arch3_03
#Files = [['IntanSound/Arch3n3-40S_123259.int', [226763]],
#         ['IntanSound/Arch3n3-43S_125706.int', [604540]],
#         ['IntanSound/Arch3n3-45S_130436.int', [327500]]]

AnalysisFile = './CaMKIIaArch3_02.hdf5'

for File in Files:
#    Data = {}; Rec = File[0].split('-')[-1].split('_')[0] # Arch3_0{2,3}
    Data = {}; Rec = File[0].split('_')[0][10:] # Arch3_01
    Data[Rec] = IntanBin.IntLoad(File[0])[0]
    Override['RecS'] = Rec
    
    Clusters = DataAnalysis.ClusterizeSpks(Data, Rate, ClusterPath, AnalysisFile, 
                              ChannelMap, Override, Return=True)
    
#    Clusters = Hdf5F.LoadClusters(AnalysisFile)
    DataAnalysis.UnitsSpks(Clusters, AnalysisFile, Override)
    
    TimeBeforeTTL = 0; TimeAfterTTL = 300
    Override['TTLs'] = DataAnalysis.GenerateFakeTTLsRising(File[1], int(0.003*Rate), 200, int(0.090*Rate))
    
    DataAnalysis.UnitsPSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
              AnalogTTLs, Override)
    
    DataAnalysis.Plot.UnitsSpksPSTH(AnalysisFile, 'svg', Override)


#%% Old arch SoundLight

Here = os.getcwd(); Path = 'SepCh/'; ClusterPath = Here + '/' + Path

RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
A16 = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

ChannelMap = DataAnalysis.RemapChannels(A16['Tip'], A16['Head'], RHAHeadstage)

Rate = np.array([25000]); AnalogTTLs = True; Override={}

Files = glob('IntanSound/*LS*') + glob('IntanSound/*SL*')

AnalysisFile = './CaMKIIaArch3_02.hdf5'

for File in Files:
#    Data = {}; Rec = File.split('-')[-1].split('_')[0]
    Data = {}; Rec = File.split('_')[0][10:] # Arch3_01
    Data[Rec] = IntanBin.IntLoad(File)[0]
    Override['RecS'] = Rec
    
    Clusters = DataAnalysis.ClusterizeSpks(Data, Rate, ClusterPath, AnalysisFile, 
                              ChannelMap, Override, Return=True)
    
#    Clusters = Hdf5F.LoadClusters(AnalysisFile)
    DataAnalysis.UnitsSpks(Clusters, AnalysisFile, Override)
    
    
    Starts = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Data[Rec][:, 16])
    Override['TTLs'] = DataAnalysis.GenerateFakeTTLsRising(Starts, int(0.003*Rate), 200, int(0.090*Rate))
    TimeBeforeTTL = 0; TimeAfterTTL = 300
    
    DataAnalysis.UnitsPSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
              AnalogTTLs, Override)

    DataAnalysis.Plot.UnitsSpksPSTH(AnalysisFile, 'svg', Override)


#%% Test Different BinSizes
