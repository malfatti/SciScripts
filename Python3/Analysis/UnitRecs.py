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
import Hdf5F
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

AnalysisFile = './CaMKIIaArch3_01.hdf5'
TimeBeforeTTL = 0; TimeAfterTTL = 300

for File in Files:
#    Data = {}; Rec = File[0].split('-')[-1].split('_')[0] # Arch3_0{2,3}
    Data = {}; Rec = File[0].split('_')[0][10:] # Arch3_01
    Data[Rec] = IntanBin.IntLoad(File[0])[0]
    Override['RecS'] = Rec
    
    Clusters = DataAnalysis.ClusterizeSpks(Data, Rate, ClusterPath, AnalysisFile, 
                              ChannelMap, Override, Return=True)
    
#    Clusters = Hdf5F.LoadClusters(AnalysisFile)
    DataAnalysis.UnitsSpks(Clusters, AnalysisFile, Override)
#    
    Override['TTLs'] = DataAnalysis.GenerateFakeTTLsRising(File[1], int(0.003*Rate), 200, int(0.090*Rate))
#    
    DataAnalysis.UnitsPSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
                           AnalogTTLs, Override)
    
    Override['XValues'] = np.arange(TimeBeforeTTL, TimeAfterTTL, 1)
    DataAnalysis.Plot.UnitsSpksPSTH(AnalysisFile, 'svg', Override)
    
    XValuesList = [np.arange(TimeBeforeTTL, TimeAfterTTL, _) for _ in [0.3, 1, 3, 5]]
    DataAnalysis.Plot.UnitsPSTHTestBin(XValuesList, AnalysisFile, 'svg', Override)


#%% Old arch SoundLight

Here = os.getcwd(); Path = 'SepCh/'; ClusterPath = Here + '/' + Path

RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
A16 = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

ChannelMap = DataAnalysis.RemapChannels(A16['Tip'], A16['Head'], RHAHeadstage)

Rate = np.array([25000]); AnalogTTLs = True; Override={}

Files = glob('Intan/*LS*') + glob('Intan/*SL*')

AnalysisFile = './CaMKIIaArch3_01.hdf5'
TimeBeforeTTL = 0; TimeAfterTTL = 300

for File in Files:
    Data = {}; Rec = File.split('-')[-1].split('_')[0]
#    Data = {}; Rec = File.split('_')[0][10:] # Arch3_01
    Data[Rec] = IntanBin.IntLoad(File)[0]
    Override['RecS'] = Rec
    
    Clusters = DataAnalysis.ClusterizeSpks(Data, Rate, ClusterPath, AnalysisFile, 
                              ChannelMap, Override, Return=True)
    
#    Clusters = Hdf5F.LoadClusters(AnalysisFile)
    DataAnalysis.UnitsSpks(Clusters, AnalysisFile, Override)
#    
    Starts = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Data[Rec][:, 16])
    Override['TTLs'] = DataAnalysis.GenerateFakeTTLsRising(Starts, int(0.003*Rate), 200, int(0.090*Rate))
    
    DataAnalysis.UnitsPSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
              AnalogTTLs, Override)
    
    Override['XValues'] = np.arange(TimeBeforeTTL, TimeAfterTTL, 1)
    DataAnalysis.Plot.UnitsSpksPSTH(AnalysisFile, 'svg', Override)
    
    XValuesList = [np.arange(TimeBeforeTTL, TimeAfterTTL, _) for _ in [0.3, 1, 3, 5]]
    DataAnalysis.Plot.UnitsPSTHTestBin(XValuesList, AnalysisFile, 'svg', Override)


#%% CaMKIIaArch3_04
Here = os.getcwd(); Path = 'SepCh/'; ClusterPath = Here + '/' + Path

RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
A16 = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

ChannelMap = DataAnalysis.RemapChannels(A16['Tip'], A16['Head'], RHAHeadstage)

Rate = np.array([25000]); AnalogTTLs = True; Override={}

Files = glob('Intan/*.int'); Files.sort()

AnalysisFile = './CaMKIIaArch3_04.hdf5'
TimeBeforeTTL = 0; TimeAfterTTL = 300

for File in Files:
    Data = {}; Rec = File.split('-')[-1].split('_')[0]
    Data[Rec] = IntanBin.IntLoad(File)[0]
    Override['RecS'] = Rec
    
    Clusters = DataAnalysis.ClusterizeSpks(Data, Rate, ClusterPath, AnalysisFile, 
                              ChannelMap, Override, Return=True)
    
#    Clusters = Hdf5F.LoadClusters(AnalysisFile)
    DataAnalysis.UnitsSpks(Clusters, AnalysisFile, Override)
    
    Threshold = 2.5
    Override['TTLs'] = [Ind for Ind, El in enumerate(Data[Rec][:, 16]) 
                            if El > Threshold 
                            if Data[Rec][:, 16][Ind-1] < Threshold]
    
    DataAnalysis.UnitsPSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
              AnalogTTLs, Override)
    
    Override['XValues'] = np.arange(TimeBeforeTTL, TimeAfterTTL, 1)
    DataAnalysis.Plot.UnitsSpksPSTH(AnalysisFile, 'svg', Override)
    
    XValuesList = [np.arange(TimeBeforeTTL, TimeAfterTTL, _) for _ in [0.3, 1, 3, 5]]
    DataAnalysis.Plot.UnitsPSTHTestBin(XValuesList, AnalysisFile, 'svg', Override)


#%% Units per Rec per Ch

import pandas as pd

Title = 'CaMKIIaeArch3_04'
Index = ['Ch', 'UnitsNo', 'Responding']

Table = {'CaMKIIaeArch3_04': 
            [[1, 1, 1],
             [2, 0, 0],
             [3, 1, 1],
             [4, 1, 0],
             [5, 0, 0],
             [6, 2, 0],
             [7, 1, 0],
             [8, 1, 1],
             [9, 1, 0],
             [10, 1, 0],
             [11, 2, 1],
             [12, 2, 0],
             [13, 3, 0],
             [14, 1, 0],
             [15, 2, 0],
             [16, 0, 0]],
        
        'CaMKIIaeArch3_03':
            
