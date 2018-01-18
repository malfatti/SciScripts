#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 16:40:55 2017

@author: malfatti
"""
#%% ABRs
from DataAnalysis import ABRs
from DataAnalysis.Plot import ABRs as ABRPlot
from copy import deepcopy
from glob import glob
from IO import Hdf5, Txt

# cd '/home/cerebro/Malfatti/Data'
#%% 
Group = 'RecoveryControl'
AnalysisFile = Group + '/' + Group + '-Analysis.hdf5'


Exps = sorted(glob(Group+'/*2018*ABRs'))
# del(Exps[5], Exps[3], Exps[2], Exps[1]) # Temp override

StimType = ['Sound']

for Exp in Exps:
    Folders = glob(Exp + '/*' + Exp.split('/')[-1][:4] + '-*'); Folders.sort()
    InfoFile = glob(Exp + '/' + Exp.split('/')[-1][:4] + '*dict')[0]
    
    AnalysisPath = '-'.join(InfoFile.split('/')[-1].split('-')[:2]) + '-ABRs'
    
    ABRs.Analysis(Exp, Folders, InfoFile, AnalysisFile, StimType=StimType)
    ABRPlot.Traces(AnalysisPath, AnalysisFile, InfoFile, Group+'/Figs', Save=True, Visible=True)


#%% Convert hdf5 info to dict
from DataAnalysis.Plot import Plot
from IO import OpenEphys
Files = glob('Prevention/*ABRs/2*.hdf5'); Files.sort()

for File in Files:
    print(File)
    SoundAmpF, DataInfo = Hdf5.DataLoad('/DataInfo', File)
#    SoundAmpF = Hdf5.DataLoad('/DataInfo', File)[0]
    DataInfo.update(SoundAmpF)
    DataInfo['ExpInfo'] = Hdf5.DataLoad('/ExpInfo', File)[1]
    
    Keys = list(DataInfo['ExpInfo'].keys())
    for E in Keys:
        if E == "{0:02d}".format(int(E)): continue 
        DataInfo['ExpInfo']["{0:02d}".format(int(E))] = deepcopy(DataInfo['ExpInfo'][E])
        del(DataInfo['ExpInfo'][E])
    
    for E in DataInfo['ExpInfo'].keys():
        DataInfo['ExpInfo'][E]['StimType'] = [_.decode() for _ in DataInfo['ExpInfo'][E]['StimType']]
        
        for S, Stim in enumerate(DataInfo['ExpInfo'][E]['StimType']):
            if '_' in Stim:  DataInfo['ExpInfo'][E]['StimType'] = Stim.split('_')
    
#    if not glob('/'.join(File.split('/')[:-1]+[File.split('/')[-1][:-4]+'dict'])):
    Folder = sorted(glob('/'.join(File.split('/')[:-1])+'/2???-*'))[-1]
    Data, Rate = OpenEphys.DataLoader(Folder, AnalogTTLs=True, Unit='Bits')
    if len(Data.keys()) == 1: Proc = list(Data.keys())[0]
    
    Chs = [Data[Proc]['3'][:int(Rate[Proc]/2),Ch] for Ch in range(Data[Proc]['3'].shape[1])]
    Plot.RawCh(Chs, Lines=len(Chs), Cols=1, Save=False)
    
    ABRCh = input('ABRCh: '); ABRCh = [int(ABRCh)]
    TTLCh = input('TTLCh: '); TTLCh = int(TTLCh)
#    else:
#        DataDict = Txt.DictRead(File[:-4]+'dict')
#        ABRCh, TTLCh = DataDict['ABRCh'], DataDict['TTLCh']
    
    DataInfo.update({'ABRCh': ABRCh, 'TTLCh': TTLCh, 'FileName': File[:-4]+'dict'})
    Txt.DictWrite(File[:-4]+'dict', DataInfo)
