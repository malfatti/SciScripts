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
    
    XValuesList = [np.arange(TimeBeforeTTL, TimeAfterTTL, _) for _ in [0.3, 1, 2, 3]]
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
#    Data = {}; Rec = File.split('-')[-1].split('_')[0]
    Data = {}; Rec = File.split('_')[0][10:] # Arch3_01
#    Data[Rec] = IntanBin.IntLoad(File)[0]
    Override['RecS'] = Rec
#    
#    Clusters = DataAnalysis.ClusterizeSpks(Data, Rate, ClusterPath, AnalysisFile, 
#                              ChannelMap, Override, Return=True)
#    
##    Clusters = Hdf5F.LoadClusters(AnalysisFile)
#    DataAnalysis.UnitsSpks(Clusters, AnalysisFile, Override)
##    
    Starts = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Data[Rec][:, 16])
    Override['TTLs'] = DataAnalysis.GenerateFakeTTLsRising(Starts, int(0.003*Rate), 200, int(0.090*Rate))
#    
#    DataAnalysis.UnitsPSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
#              AnalogTTLs, Override)
#    
#    Override['XValues'] = np.arange(TimeBeforeTTL, TimeAfterTTL, 1)
#    DataAnalysis.Plot.UnitsSpksPSTH(AnalysisFile, 'svg', Override)
    
    XValuesList = [np.arange(TimeBeforeTTL, TimeAfterTTL, _) for _ in [0.3, 1, 2, 3]]
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
    if Rec[-1] == 'L': continue
    Data[Rec] = IntanBin.IntLoad(File)[0]
    Override['RecS'] = Rec
    
#    Clusters = DataAnalysis.ClusterizeSpks(Data, Rate, ClusterPath, AnalysisFile, 
#                              ChannelMap, Override, Return=True)
    
    Clusters = Hdf5F.LoadClusters(AnalysisFile)
    DataAnalysis.UnitsSpks(Clusters, AnalysisFile, Override)
    
    Threshold = 1.5
    Override['TTLs'] = [Ind for Ind, El in enumerate(Data[Rec][:, 16]) 
                            if El > Threshold 
                            if Data[Rec][:, 16][Ind-1] < Threshold]
    
    DataAnalysis.UnitsPSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
              AnalogTTLs, Override)
    
    Override['XValues'] = np.arange(TimeBeforeTTL, TimeAfterTTL, 1)
    DataAnalysis.Plot.UnitsSpksPSTH(AnalysisFile, 'svg', Override)
    
    XValuesList = [np.arange(TimeBeforeTTL, TimeAfterTTL, _) for _ in [0.3, 1, 2, 3]]
    DataAnalysis.Plot.UnitsPSTHTestBin(XValuesList, AnalysisFile, 'svg', Override)


#%% Units per Rec per Ch

import numpy as np
import pandas as pd
from copy import deepcopy
from datetime import datetime

Index = ['Ch', 'UnitsNo', 'Responding']

Table = {'CaMKIIaeArch3_04': {
             '2000': {
                 # Ch11 maybe
                 'Sound-Light': [
                     [1, 1, 1],
                     [3, 1, 1],
                     [4, 1, 0],
                     [6, 2, 0],
                     [7, 1, 0],
                     [8, 1, 1],
                     [9, 1, 0],
                     [10, 1, 0],
                     [11, 2, 1], # Nice
                     [12, 2, 0],
                     [13, 3, 0],
                     [14, 1, 0],
                     [15, 2, 0]],
                'Sound':[
                     [1, 1, 1],
                     [4, 1, 1],
                     [6, 1, 1],
                     [7, 1, 1],
                     [9, 1, 1],
                     [10, 1, 1],
                     [11, 1, 1], # Best 
                     [12, 1, 0],
                     [13, 1, 0],
                     [14, 1, 1], # Ugly
                     [15, 1, 0]]
                     },
             '3500': {
                 'Sound-Light': [
                     [1, 2, 0],
                     [3, 2, 0],
                     [4, 2, 0],
                     [5, 2, 0],
                     [6, 2, 0],
                     [7, 1, 1],
                     [8, 2, 0],
                     [10, 2, 0],
                     [12, 3, 0],
                     [13, 1, 0],
                     [15, 2, 0],
                     [16, 2, 0]],
                 'Sound': [
                     [1, 1, 1], # Each 50ms?
                     [3, 1, 1],
                     [4, 1, 1], # Nice
                     [5, 1, 1], # Best ever, reproduces masking result :D
                     [6, 1, 1], # Nice +-
                     [8, 2, 1], # Resembles Ch5
                     [10, 1, 1],
                     [12, 2, 1],
                     [13, 1, 1],
                     [15, 1, 1]] # Sine wave-like response?
                     },
             '4300': {
                 'Sound-Light': [
                     [1, 2, 0],
                     [3, 1, 0],
                     [4, 1, 0],
                     [5, 1, 1],
                     [6, 1, 1],
                     [7, 3, 1],
                     [9, 1, 1],
                     [10, 1, 1],
                     [11, 1, 0],
                     [14, 1, 0]],
                 'Sound': [
                     [1, 2, 0],
                     [3, 2, 1],
                     [4, 1, 1],
                     [5, 1, 1], # Nice
                     [6, 1, 1], # Ugly
                     [7, 2, 0],
                     [8, 1, 1],
                     [9, 2, 0],
                     [10, 1, 1],
                     [11, 1, 1],
                     [14, 1, 0],
                     [15, 1, 1]] # Ugly
                     },
             '4500': {
                 # Ch06,08 maybe
                 'Sound-Light': [
                     [1, 1, 1],
                     [3, 1, 1],
                     [4, 1, 1],
                     [5, 2, 1],
                     [6, 1, 0],
                     [7, 1, 1],
                     [8, 4, 1],
                     [9, 1, 0],
                     [11, 2, 1],
                     [12, 1, 1],
                     [14, 1, 0],
                     [15, 1, 0]],
                 'Sound': [
                     [1, 1, 1],
                     [3, 1, 1],
                     [4, 1, 1],
                     [5, 2, 0],
                     [6, 1, 1],
                     [7, 1, 0],
                     [8, 2, 2], # Sine wave-like response?
                     [9, 1, 1],
                     [10, 1, 0],
                     [11, 2, 1],
                     [12, 1, 1],
                     [14, 1, 1],
                     [15, 1, 0]]
                     }
                 },
        'CaMKIIaeArch3_01': {
             '3670': {
                 'Sound-Light': [
                     [6, 1, 1],
                     [12, 1, 1],
                     [13, 2, 0],
                     [15, 2, 2],
                     ]
             },
             '4000': {
                 'Sound-Light': [
                     [4, 1, 1],
                     [5, 3, 2],
                     [6, 1, 1],
                     [7, 4, 3],
                     [8, 1, 1],
                     [9, 3, 1],
                     [10, 2, 1],
                     [11, 6, 3], # YeahRight
                     [12, 2, 2],
                     [13, 2, 1],
                     [14, 2, 1],
                     [15, 2, 2]] # Weak
             },
             '4300': {
                 'Sound-Light': [
                     [3, 1, 1],
                     [4, 1, 1],
                     [5, 4, 2],
                     [6, 1, 1],
                     [7, 3, 2],
                     [8, 3, 1],
                     [9, 2, 1],
                     [10, 4, 3],
                     [11, 2, 2], # Nice
                     [12, 1, 1], # Suspiciously like sine wave...
                     [13, 3, 2],
                     [14, 1, 0],
                     [15, 1, 1]] 
             }
        },
        'CaMKIIaeArch3_02': {
             '3500': {
                 # Ch01,04,07-,13+
                 'Sound-Light': [
                     [1, 1, 0],
                     [3, 1, 1],
                     [4, 1, 0],
                     [5, 1, 0],
                     [6, 2, 0],
                     [7, 1, 1],
                     [8, 1, 1],
                     [9, 1, 1],
                     [10, 1, 1], # Better than SoundOnly
                     [11, 1, 0],
                     [12, 1, 0],
                     [13, 1, 0],
                     [15, 1, 1]],
                'Sound': [
                     [1, 1, 1],
                     [3, 1, 0],
                     [4, 2, 1],
                     [5, 1, 0],
                     [6, 1, 0],
                     [7, 1, 1],
                     [8, 1, 1],
                     [9, 1, 0],
                     [10, 1, 1],
                     [11, 1, 0],
                     [12, 1, 0],
                     [13, 1, 1],
                     [14, 1, 0],
                     [15, 1, 1]]
             },
             '4300': {
                 #Ch13
                 'Sound-Light': [
                     [1, 1, 0],
                     [3, 1, 0],
                     [4, 1, 0],
                     [5, 2, 0],
                     [6, 1, 0],
                     [7, 1, 1],
                     [8, 1, 1],
                     [9, 1, 0],
                     [10, 1, 0],
                     [11, 1, 0],
                     [12, 1, 0],
                     [13, 2, 0],
                     [14, 1, 1],
                     [15, 2, 0]],
                'Sound': [
                     [1, 1, 0],
                     [3, 1, 0],
                     [4, 2, 0],
                     [5, 1, 0],
                     [6, 2, 0],
                     [7, 1, 1],
                     [8, 1, 1],
                     [9, 1, 1],
                     [10, 1, 1],
                     [11, 1, 0],
                     [12, 1, 0],
                     [13, 3, 2],
                     [14, 1, 1],
                     [15, 1, 0]]
            },
             '4500': {
                 #Ch12-
                 'Sound-Light': [
                     [1, 2, 0],
                     [3, 1, 0],
                     [5, 1, 0],
                     [6, 3, 1],
                     [7, 1, 1],
                     [8, 1, 1],
                     [9, 1, 1],
                     [10, 1, 0],
                     [11, 1, 0],
                     [12, 2, 2],
                     [13, 2, 2],
                     [14, 1, 1],
                     [15, 3, 3]],
                'Sound': [
                     [1, 2, 0],
                     [3, 1, 1],
                     [5, 1, 0],
                     [6, 3, 2],
                     [7, 1, 1],
                     [8, 1, 1],
                     [9, 1, 1],
                     [10, 1, 1],
                     [11, 1, 1],
                     [12, 3, 2],
                     [13, 3, 0],
                     [14, 1, 1],
                     [15, 3, 2]]
            },
    }
}
        
#        'CaMKIIaeArch3_03':
#Now = datetime.now().strftime("%Y%m%d%H%M%S")
#TxtFile = open('UnitRecs_'+Now+'.txt', 'w')
TxtFile = open('/home/malfatti/Nebula/Documents/PhD/Notes/Experiments/UnitRecs', 'w')

Animals = list(Table.keys()); Animals.sort()
for Animal in Animals:
    Coords = list(Table[Animal].keys()); Coords.sort()
    TxtFile.write('""" ' + Animal + ' """')
    
    for Coord in Coords:
        Stims = list(Table[Animal][Coord].keys()); Stims.sort()
        TxtFile.write(r'''
''')
        TxtFile.write('## ' + Coord + '\n')
        
        for Stim in Stims:
            if not Table[Animal][Coord][Stim]: continue
            Data = pd.DataFrame(data=Table[Animal][Coord][Stim], columns=Index)
            
            TxtFile.write(Stim + '\n')
            TxtFile.write(Data.to_string())
            TxtFile.write(r'''
''')

TxtFile.close()

NPTable = deepcopy(Table)
Animals = list(Table.keys()); Animals.sort()
for Animal in Animals:
    Coords = list(Table[Animal].keys()); Coords.sort()
    
    for Coord in Coords:
        Stims = list(Table[Animal][Coord].keys()); Stims.sort()
        
        for Stim in Stims:
            if not Table[Animal][Coord][Stim]: continue
            
            NPTable[Animal][Coord][Stim] = np.array(Table[Animal][Coord][Stim])

## Total Units
Total = {}
Animals = list(Table.keys()); Animals.sort()
for Animal in Animals:
    Coords = list(Table[Animal].keys()); Coords.sort()
    Total[Animal] = {}
    
    for Coord in Coords:
        Stims = list(Table[Animal][Coord].keys()); Stims.sort()
        Total[Animal][Coord] = {}
        
        for Stim in Stims:
            if not Table[Animal][Coord][Stim]: continue
            
            Total[Animal][Coord][Stim] = [len(NPTable[Animal][Coord][Stim][:,0]),
                                          sum(NPTable[Animal][Coord][Stim][:,1]),
                                          sum(NPTable[Animal][Coord][Stim][:,2])]

## No of Ch with units
Animals = list(Table.keys()); Animals.sort()
for Animal in Animals:
    Coords = list(Table[Animal].keys()); Coords.sort()
    Total[Animal] = {}
    
    for Coord in Coords:
        Stims = list(Table[Animal][Coord].keys()); Stims.sort()
        if len(Stims) > 2: continue
        Total[Animal][Coord] = {}
        
        for Stim in Stims:
            if not Table[Animal][Coord][Stim]: continue
            
            