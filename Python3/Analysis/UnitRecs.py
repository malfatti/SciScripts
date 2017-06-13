#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 22:25:17 2017

@author: malfatti
"""
#%% Setting
#import DataAnalysis, Hdf5F, IntanBin, Klusta, os, TarBinPy
import os
import numpy as np
#import pandas as pd

from DataAnalysis import DataAnalysis, Plot, Units
from glob import glob
from IO import Bin, Hdf5, Intan, Klusta
from IO.Txt import DictRead
#from klusta.kwik import KwikModel

#Params = {'backend': 'TkAgg'}
#from matplotlib import rcParams; rcParams.update(Params)
#from matplotlib.gridspec import GridSpec
#from matplotlib import pyplot as plt


#%% Klusta
Board = 'OE'; Rec='0'
TimeBeforeTTL = 0; TimeAfterTTL = 300; BinSize = 3

CustomAdaptor = [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4]
A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
Map = DataAnalysis.RemapChannels(A16['ProbeTip'], A16['ProbeHead'], CustomAdaptor)
Spacing = 25
PrbFile = os.environ['SCRIPTSPATH'] + '/Python3/Klusta/A16-'+str(Spacing)+'.prb'

Folders = glob('*/**/*UnitRec/**/*.kwd', recursive=True)
Folders = ['/'.join(F.split('/')[:-1]) for F in Folders]
Folders = DataAnalysis.UniqueStr(Folders)

Done = []; Errors = []; ErrorsLog = []; Skipped = []
Done = [
'CaMKIIahM3Dn01-20150903-UnitRec/CaMKIIa-hM3D-LeftEar_2015-09-03_15-39-51_1800/KlustaFiles',
'CaMKIIahM3Dn01-20150903-UnitRec/CaMKIIa-hM3D-LeftEar_2015-09-03_19-23-00_4500-CNO06-Kwik/KlustaFiles',
'CaMKIIahM4Dn04-20151012-UnitRec/Test_2015-10-12_07-22-14_Kwik/KlustaFiles',
'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160703-UnitRec/KwikFiles/2016-07-03_19-03-56_NaCl/KlustaFiles',
'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160703-UnitRec/KwikFiles/2016-07-03_19-08-21_NaCl/KlustaFiles',
'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160703-UnitRec/KwikFiles/2016-07-03_19-12-50_NaCl/KlustaFiles',
'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160703-UnitRec/KwikFiles/2016-07-03_19-18-49_NaCl/KlustaFiles',
'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160703-UnitRec/KwikFiles/2016-07-03_19-23-16_NaCl/KlustaFiles',
'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160703-UnitRec/KwikFiles/2016-07-03_19-56-02_CNO/KlustaFiles',
'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160703-UnitRec/KwikFiles/2016-07-03_20-00-26_CNO/KlustaFiles',
'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160703-UnitRec/KwikFiles/2016-07-03_20-04-59_CNO/KlustaFiles',
'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160703-UnitRec/KwikFiles/2016-07-03_20-09-31_CNO/KlustaFiles',
'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160703-UnitRec/KwikFiles/2016-07-03_20-13-54_CNO/KlustaFiles',
'CaMKIIahM4Dn09/CaMKIIahM4Dn09-20160703-UnitRec/KwikFiles/2016-07-04_10-14-32_NaCl/KlustaFiles',
'CaMKIIahM4Dn09/CaMKIIahM4Dn09-20160703-UnitRec/KwikFiles/2016-07-04_10-19-01_NaCl/KlustaFiles',
'CaMKIIahM4Dn09/CaMKIIahM4Dn09-20160703-UnitRec/KwikFiles/2016-07-04_10-23-27_NaCl/KlustaFiles',
'CaMKIIahM4Dn09/CaMKIIahM4Dn09-20160703-UnitRec/KwikFiles/2016-07-04_10-27-52_NaCl/KlustaFiles',
'CaMKIIahM4Dn09/CaMKIIahM4Dn09-20160703-UnitRec/KwikFiles/2016-07-04_10-32-51_NaCl/KlustaFiles',
'CaMKIIahM4Dn09/CaMKIIahM4Dn09-20160703-UnitRec/KwikFiles/2016-07-04_11-37-17_CNO/KlustaFiles',
'EarBarTest/EarBarTest_02-20170217-UnitRec/2017-02-17_13-59-06_Ear/KlustaFiles'
]
Done = ['/'.join(F.split('/')[:-1]) for F in Done]

for Folder in Folders:
    try:
        if Folder in Done+Errors+Skipped: continue
        
        Data = Hdf5.LoadOEKwik(Folder, True, 'Bits', Map)[0]
        Proc = Hdf5.GetProc(Data, Board)
        
        Keys = list(Data[Proc]['data'].keys())
        if Data[Proc]['data'][Keys[0]].shape[1] < 16: Skipped.append(Folder); continue
        
        ExpFolder = Folder + '/KlustaFiles'; os.makedirs(ExpFolder, exist_ok=True)
        ExpName = Folder.split('/')[-1]
        for R, Rec in Data[Proc]['data'].items():
            DataInfo = {'Rate': int(Data[Proc]['info'][R]['sample_rate'])}
            DataFile = ExpName + '_Rec' + "{0:02d}".format(int(R))
            Bin.Write(DataFile+'.dat', ExpFolder, Rec, DataInfo)
        
        DataInfo = DictRead(ExpFolder+'/'+DataFile+'-Info.dict')
        raw_data_files = glob(os.getcwd()+'/'+ExpFolder+'/'+ExpName+'*.dat')
        raw_data_files.sort()
        
        Klusta.PrmWrite(ExpFolder+'/'+ExpName+'.prm', ExpName, PrbFile, raw_data_files, 
                        DataInfo['Rate'],DataInfo['Shape'][1], DataInfo['DType'])
        
        Klusta.Run(ExpName+'.prm', os.getcwd()+'/'+ExpFolder, Overwrite=True)
        Done.append(Folder)
    
    except Exception as e:
        Errors.append(Folder)
        ErrorsLog.append(e)
        print(e)


#%% Batch
Animal = 'CaMKIIahM4Dn09'
Exp = 'CaMKIIahM4Dn09-20160703-UnitRec'

Board = 'OE'
AnalogTTLs = True
StimTTLCh = 17
TimeBeforeTTL = 0
TimeAfterTTL = 300

CustomAdaptor = [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4]
A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

#==========#==========#==========#==========#

ChannelMap = DataAnalysis.RemapChannels(A16['ProbeTip'], A16['ProbeHead'], CustomAdaptor)
RecFolders = glob(Animal + '/' + Exp + '/KwikFiles/*'); RecFolders.sort()
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'
Done, Errors, ErrorsLog = [], [], []

for RecFolder in RecFolders:
    if RecFolder in Done or RecFolder in Errors: 
        print('!!!==========')
        print(RecFolder, 'already done. Skipping...')
        print('!!!==========')
        continue
    
    try:
        ExpFolder = RecFolder.split('/')[-1]
        ClusterPath = os.getcwd() + '/' + RecFolder + '/' + 'ClusterFiles'
        AnalysisKey = '/'.join(ClusterPath.split('/')[-4:-1:2])
        
        Data = Hdf5.OEKwikLoad(RecFolder, AnalogTTLs, 'uV', ChannelMap)[0]
        Proc = Hdf5.GetProc(Data, Board)
        
        for Rec in Data[Proc]['data'].keys():
            Rate = Data[Proc]['info'][Rec]['sample_rate']
            Clusters = Units.WaveClus.ClusterizeSpks(Data[Proc]['data'][Rec], Rate, 
                                                   ChannelMap, ClusterPath, AnalysisFile, 
                                                   AnalysisKey, Rec, Return=True)
            
            TTLCh = Data[Proc]['data'][Rec][:, StimTTLCh-1]
            
            Rec = "{0:02d}".format(int(Rec))
            Units.WaveClus.Spks(Clusters, AnalysisFile, AnalysisKey, Rec)
            Units.WaveClus.PSTH(Clusters, TTLCh, Rate, AnalysisFile, AnalysisKey, 
                                   Rec, TimeBeforeTTL, TimeAfterTTL, AnalogTTLs)
            
            Override = {'Rec': Rec}
            XValues = np.arange(TimeBeforeTTL, TimeAfterTTL, 1)
            UnitsRec = Hdf5.LoadUnits(AnalysisFile, AnalysisKey, Override)[0]
            FigBase = RecFolder + '/Figs/' + Exp + '_Rec' + Rec
            
            Plot.Units.SpksPSTH(UnitsRec, XValues, FigBase, Ext='svg')
            
            del(Clusters)
        
        del(Data); Done.append(RecFolder)
    
    except Exception as e:
        Errors.append(RecFolder); ErrorsLog.append(e)
        print('!!!==========')
        print(e)
        print('!!!==========')


#%% Clustering
Animal = 'CaMKIIahM4Dn08'
Exp = 'CaMKIIahM4Dn08-20160703-UnitRec'
RecFolder = 'KwikFiles/2016-07-03_19-03-56_NaCl'

AnalogTTLs = True
Board = 'OE'

CustomAdaptor = [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4]
A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

#==========#==========#==========#==========#

DataFolder = Animal + '/' + Exp + '/' + RecFolder
ClusterPath = os.getcwd() + '/' + DataFolder + '/' + 'ClusterFiles'
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'
AnalysisKey = '/'.join(ClusterPath.split('/')[-4:-1:2])

ChannelMap = DataAnalysis.RemapChannels(A16['ProbeTip'], A16['ProbeHead'], CustomAdaptor)

Data = Hdf5.LoadOEKwik(DataFolder, AnalogTTLs, 'uV', ChannelMap)[0]
Proc = Hdf5.GetProc(Data, Board)

for Rec in Data[Proc]['data'].keys():
    Rate = Data[Proc]['info'][Rec]['sample_rate']
    Clusters = Units.WaveClus.ClusterizeSpks(Data[Proc]['data'][Rec], Rate, 
                                           ChannelMap, ClusterPath, 
                                           AnalysisFile, AnalysisKey, Rec, 
                                           Return=True)


#%% Units+Spks
Animal = 'CaMKIIahM4Dn08'
Exp = 'CaMKIIahM4Dn08-20160703-UnitRec'
RecFolder = 'KwikFiles/2016-07-03_19-03-56_NaCl'

Board = 'OE'
AnalogTTLs = True
StimTTLCh = 17
TimeBeforeTTL = 0
TimeAfterTTL = 300
StimType0 = ['Sound_NaCl']
StimType1 = ['Sound_CNO']

CustomAdaptor = [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4]
A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

#==========#==========#==========#==========#

DataFolder = Animal + '/' + Exp + '/' + RecFolder
ClusterPath = os.getcwd() + '/' + DataFolder + '/' + 'ClusterFiles'
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'
AnalysisKey = '/'.join(ClusterPath.split('/')[-4:-1:2])

ChannelMap = DataAnalysis.RemapChannels(A16['ProbeTip'], A16['ProbeHead'], 
                                        CustomAdaptor)

Data = Hdf5.LoadOEKwik(DataFolder, AnalogTTLs, 'uV', ChannelMap)[0]
Proc = Hdf5.GetProc(Data, Board)
Clusters = Hdf5.LoadClusters(AnalysisFile, AnalysisKey)

for Rec in Clusters.keys():
    ORec = "{0:01d}".format(int(Rec))
    Rate = Data[Proc]['info'][ORec]['sample_rate']
    TTLCh = Data[Proc]['data'][ORec][:, StimTTLCh-1]
     
    Units.WaveClus.Spks(Clusters, AnalysisFile, AnalysisKey, Rec)
    Units.WaveClus.PSTH(Clusters, TTLCh, Rate, AnalysisFile, AnalysisKey, 
                           Rec, TimeBeforeTTL, TimeAfterTTL, AnalogTTLs)
    
    Override = {'Rec': Rec}
    XValues = np.arange(TimeBeforeTTL, TimeAfterTTL, 1)
    UnitsRec = Hdf5.LoadUnits(AnalysisFile, AnalysisKey, Override)[0]
    FigBase = DataFolder + '/Figs/' + Exp + '_Rec' + Rec
            
    Plot.Units.SpksPSTH(UnitsRec, XValues, FigBase, Ext='svg')


#%% Old Arch recs w/o TTLs
Animal = 'CaMKIIaArch3_02'
Here = os.getcwd(); Path = 'SepCh/'
ClusterPath = Here + '/' + Path

RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
A16 = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

ChannelMap = DataAnalysis.RemapChannels(A16['Tip'], A16['Head'], RHAHeadstage)

Rate = np.array([25000]); AnalogTTLs = True; Override={}

## Arch3_01
#Files = [['Intan/Arch45S_153244.int', [1652313]]]
### Arch3_02
Files = [#['IntanSound/Arch3n2-43S_184905.int', [154934, 4121425]],
         ['IntanSound/Arch3n2-35S_175902.int', [1702747]],
#         'IntanSound/Arch3n2-40S_184417.int',
         ]#['IntanSound/Arch3n2-45S_191253.int', [1283124]]]
## Arch3_03
#Files = [['IntanSound/Arch3n3-40S_123259.int', [226763]],
#         ['IntanSound/Arch3n3-43S_125706.int', [604540]],
#         ['IntanSound/Arch3n3-45S_130436.int', [327500]]]

AnalysisFile = './'+Animal+'.hdf5'
TimeBeforeTTL = 0; TimeAfterTTL = 300

for File in Files:
    Data = {}; Rec = File[0].split('-')[-1].split('_')[0] # Arch3_0{2,3}
#    Data = {}; Rec = File[0].split('_')[0][10:] # Arch3_01
    Data[Rec] = Intan.IntLoad(File[0])[0]
    Override['RecS'] = Rec
    
    Clusters = Units.WaveClus.ClusterizeSpks(Data, Rate, ClusterPath, AnalysisFile, 
                              ChannelMap, Override, Return=True)
    
#    Clusters = Hdf5.LoadClusters(AnalysisFile)
    Units.WaveClus.Spks(Clusters, AnalysisFile, Override)
#    
    Override['TTLs'] = DataAnalysis.GenerateFakeTTLsRising(File[1], int(0.003*Rate), 200, int(0.090*Rate))
#    
    Units.WaveClus.PSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
                           AnalogTTLs, Override)
    
    XValues = np.arange(TimeBeforeTTL, TimeAfterTTL, 1)
    UnitsRec = Hdf5.LoadUnits(AnalysisFile, AnalysisKey, Override)[0]
    FigBase = './Figs/' + Animal + '_Rec' + Rec
    Plot.Units.SpksPSTH(UnitsRec, XValues, FigBase, Ext='svg')
    
    XValues = [np.arange(TimeBeforeTTL, TimeAfterTTL, _) for _ in [0.3, 1, 2, 3]]
    Plot.Units.SpksPSTH(UnitsRec, XValues, FigBase, Mode='BinSizeTest', Ext='svg')


#%% Old arch SoundLight
Animal = 'CaMKIIaArch3_01'
Here = os.getcwd(); Path = 'SepCh/'; ClusterPath = Here + '/' + Path

RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
A16 = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

ChannelMap = DataAnalysis.RemapChannels(A16['Tip'], A16['Head'], RHAHeadstage)

Rate = np.array([25000]); AnalogTTLs = True; Override={}

Files = glob('Intan/*LS*') + glob('Intan/*SL*')

AnalysisFile = './'+Animal+'.hdf5'
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
##    Clusters = Hdf5.LoadClusters(AnalysisFile)
#    DataAnalysis.UnitsSpks(Clusters, AnalysisFile, Override)
##    
    Starts = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Data[Rec][:, 16])
    Override['TTLs'] = DataAnalysis.GenerateFakeTTLsRising(Starts, int(0.003*Rate), 200, int(0.090*Rate))
    FigBase = './Figs/' + Animal + '_Rec' + Rec
#    
#    DataAnalysis.UnitsPSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
#              AnalogTTLs, Override)
#    XValues = np.arange(TimeBeforeTTL, TimeAfterTTL, 1)
#    UnitsRec = Hdf5.LoadUnits(AnalysisFile, AnalysisKey, Override)[0]
#    DataAnalysis.Plot.UnitsSpksPSTH(UnitsRec, XValues, FigBase, Ext='svg')
    
    XValues = [np.arange(TimeBeforeTTL, TimeAfterTTL, _) for _ in [0.3, 1, 2, 3]]
    Plot.Units.SpksPSTH(UnitsRec, XValues, FigBase, Mode='BinSizeTest', Ext='svg')


#%% CaMKIIaArch3_04
Animal = 'CaMKIIaArch3_04'
Here = os.getcwd(); Path = 'SepCh/'; ClusterPath = Here + '/' + Path

RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
A16 = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

ChannelMap = DataAnalysis.RemapChannels(A16['Tip'], A16['Head'], RHAHeadstage)

Rate = np.array([25000]); AnalogTTLs = True; Override={}

Files = glob('Intan/*.int'); Files.sort()

AnalysisFile = './'+Animal+'.hdf5'
TimeBeforeTTL = 0; TimeAfterTTL = 300

for File in Files:
    Data = {}; Rec = File.split('-')[-1].split('_')[0]
    if Rec[-1] == 'L': continue
    Data[Rec] = Intan.IntLoad(File)[0]
    Override['RecS'] = Rec
    
#    Clusters = DataAnalysis.ClusterizeSpks(Data, Rate, ClusterPath, AnalysisFile, 
#                              ChannelMap, Override, Return=True)
    
    Clusters = Hdf5.LoadClusters(AnalysisFile)
    Units.WaveClus.Spks(Clusters, AnalysisFile, Override)
    
    Threshold = 1.5
    Override['TTLs'] = [Ind for Ind, El in enumerate(Data[Rec][:, 16]) 
                            if El > Threshold 
                            if Data[Rec][:, 16][Ind-1] < Threshold]
    
    Units.WaveClus.PSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
              AnalogTTLs, Override)
    
    XValues = np.arange(TimeBeforeTTL, TimeAfterTTL, 1)
    UnitsRec = Hdf5.LoadUnits(AnalysisFile, AnalysisKey, Override)[0]
    FigBase = './Figs/' + Animal + '_Rec' + Rec
    Plot.Units.SpksPSTH(UnitsRec, XValues, FigBase, Ext='svg')
    
    XValues = [np.arange(TimeBeforeTTL, TimeAfterTTL, _) for _ in [0.3, 1, 2, 3]]
    Plot.Units.SpksPSTH(UnitsRec, XValues, FigBase, Mode='BinSizeTest', Ext='svg')


#%% Plot raw channels
Animal = 'CaMKIIaArch3_02'
Ch = 13

RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
A16 = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}

ChannelMap = DataAnalysis.RemapChannels(A16['Tip'], A16['Head'], RHAHeadstage)
## Sound
#File = ['IntanSound/Arch3n2-35S_175902.int', [1702747]]
#Data = {}; Rec = File[0].split('-')[-1].split('_')[0]
#Data[Rec], XValues = IntanBin.IntLoad(File[0], ChannelMap)
#Rate = np.array([25000])
#TTLs = DataAnalysis.GenerateFakeTTLsRising(File[1], int(0.003*Rate), 200, int(0.090*Rate))
## SoundLight
File = 'IntanSound/Arch3n2-35SL_180751.int'
Data = {}; Rec = File.split('-')[-1].split('_')[0]
Data[Rec], XValues = Intan.IntLoad(File, ChannelMap)
Rate = np.array([25000]); AnalogTTLs = True; Override={}
Starts = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Data[Rec][:, 16])
TTLs = DataAnalysis.GenerateFakeTTLsRising([Starts[0]], int(0.003*Rate), 200, int(0.090*Rate))

Data[Rec][:,Ch-1] = DataAnalysis.FilterSignal(Data[Rec][:,Ch-1], Rate, [300, 3000])
TTLVec = DataAnalysis.GenerateTTLVector(TTLs, int(3*Rate/1000), len(Data[Rec]))

## Full block
#Slice = [int(TTLs[0]-(4*Rate)), int(TTLs[-1]+(4*Rate))]
#FigName = ''.join(['Figures/', Animal, '_Rec', Rec, '_RawCh', str(Ch), '.svg'])
## Inset
Slice = [int(TTLs[0]-(4*Rate))+125000, int(TTLs[0]-(4*Rate))+137500]
FigName = ''.join(['Figures/', Animal, '_Rec', Rec, '_RawCh', str(Ch), '_Inset.svg'])

Leg, Colors = ['Sound pulses', 'Channel at 4.2mm'], ['r', 'k']
DataAnalysis.Plot.RawCh([TTLVec[:,0], Data[Rec][:,Ch-1]], 2, 1, XValues, Slice, Leg, 
                   Colors, FigName=FigName)

