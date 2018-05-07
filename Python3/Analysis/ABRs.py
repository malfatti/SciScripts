#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 16:40:55 2017

@author: malfatti
"""
#%% ABRs
from DataAnalysis import ABRs, DataAnalysis
from DataAnalysis.Plot import ABRs as ABRPlot
from copy import deepcopy
from glob import glob
from IO import Hdf5, IO, Txt

# cd '/home/cerebro/Malfatti/Data'
#%% 
Group = 'ToDelete'
AnalysisFile = Group + '/' + Group + '-Analysis.hdf5'
AnalysisFile = 'Test.hdf5'


Exps = sorted(glob(Group+'/*ABRs'))
# del(Exps[5], Exps[3], Exps[2], Exps[1]) # Temp override

StimType = ['Sound']

for Exp in Exps:
    Folders = glob(Exp + '/*' + Exp.split('/')[-1][:4] + '-*'); Folders.sort()
    InfoFile = glob(Exp + '/' + Exp.split('/')[-1][:4] + '*dict')[0]
    
    AnalysisPath = '-'.join(InfoFile.split('/')[-1].split('-')[:2]) + '-ABRs'
    
    ABRs.Analysis(Folders, InfoFile, AnalysisFile, StimType=StimType)
    ABRPlot.Traces(AnalysisPath, AnalysisFile, InfoFile, Group+'/Figs', Save=False, Show=True)



#%% Single folder
Folder = '/home/malfatti/Documents/PhD/Data/Recovery/20160703-CaMKIIahM4Dn08-UnitRec/KwikFiles/2016-07-03_19-03-56_NaCl'
InfoFile = '/home/malfatti/Documents/PhD/Data/Recovery/20160703-CaMKIIahM4Dn08-UnitRec/20160703185210-CaMKIIahM4Dn08-SoundStim.dict'
AnalysisFile = 'Test.hdf5'
ABRCh = [1]
TTLCh = 0

Data, Rate = IO.DataLoader(Folder, Unit='uV', ChannelMap=[])
TTLs = DataAnalysis.QuantifyTTLsPerRec(True, Data['100']['0'][:,-1])
ABR = DataAnalysis.SliceData(Data['100']['0'][:,ABRCh[0]], TTLs, 0, 
                             int(0.01*Rate['100']), AnalogTTLs=True)

if 'plt' not in globals():
    Params = {'backend': 'Qt5Agg'}
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
for T in range(len(TTLs)): plt.plot(DataAnalysis.FilterSignal(ABR[:,T], Rate['100'], [300,3000]))
plt.plot(DataAnalysis.FilterSignal(ABR.mean(axis=1), Rate['100'], [300,3000]), 'k', lw=3)
plt.show()

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
    Data, Rate = IO.DataLoader(Folder, AnalogTTLs=True, Unit='Bits')
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

#%% Convert old dict to new dict
DataInfo['InfoFile'] = DataInfo.pop('FileName')
for K in ['Animal', 'DAqs', 'Audio', 'Laser']: DataInfo[K] = {}

for K in ['AnimalName', 'StimType']:
    if K in DataInfo: DataInfo['Animal'][K] = DataInfo.pop(K)

for K in ['SoundCh', 'TTLCh', 'ABRCh', 'BaudRate', 'AnalogTTLs']:
    if K in DataInfo: DataInfo['DAqs'][K] = DataInfo.pop(K)

for K in ['Rate', 'Intensities', 'NoiseFrequency', 'SoundPulseNo', 
          'SoundPauseBeforePulseDur', 'SoundPulseDur', 
          'SoundPauseAfterPulseDur', 'PauseBetweenIntensities', 'System', 
          'Setup', 'CalibrationFile', 'SoundAmpF']: 
    if K in DataInfo: DataInfo['Audio'][K] = DataInfo.pop(K)

for K in ['LaserStimBlockNo', 'LaserPulseNo', 'LaserPauseBeforePulseDur', 
          'LaserPulseDur', 'LaserPauseAfterPulseDur', 
          'LaserPauseBetweenStimBlocksDur', 'LaserType', 'LaserDur', 
          'LaserFreq']:
    if K in DataInfo: DataInfo['Laser'][K] = DataInfo.pop(K)

RemainingKeys = list(DataInfo.keys())
for K in RemainingKeys:
    if K == 'LaserPrePauseDur': DataInfo['Laser']['LaserPauseBeforePulseDur'] = DataInfo.pop(K)
    elif K == 'LaserPostPauseDur': DataInfo['Laser']['LaserPauseAfterPulseDur'] = DataInfo.pop(K)
    elif K == 'SoundPauseBetweenStimBlocksDur': DataInfo['Audio']['PauseBetweenIntensities'] = DataInfo.pop(K)
    elif K == 'SoundStimBlockNo': DataInfo['Audio'][K] = DataInfo.pop(K)
    elif K == 'SoundPrePauseDur': DataInfo['Audio']['SoundPauseBeforePulseDur'] = DataInfo.pop(K)
    elif K == 'SoundPostPauseDur': DataInfo['Audio']['SoundPauseAfterPulseDur'] = DataInfo.pop(K)

for K in DataInfo['ExpInfo']: DataInfo['ExpInfo'][K]['DV'] = DataInfo['ExpInfo'][K].pop('DVCoord')
