#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 16:40:55 2017

@author: malfatti
"""
#%% ABRs
import numpy as np

from DataAnalysis import ABRs, DataAnalysis
from DataAnalysis.Plot import ABRs as ABRPlot
from DataAnalysis.Plot import Plot
from copy import deepcopy
from glob import glob
from IO import Hdf5, IO, Txt

#%%
def ABRSingle(Data, Rate, ABRCh, TTLCh, SecondsBefore, SecondsAfter, FilterFreq):
    TTLs = DataAnalysis.QuantifyTTLsPerRec(True, Data[:,TTLCh-1])
    ABR = DataAnalysis.SliceData(Data[:,ABRCh[0]-1], TTLs, 
                                 int(SecondsBefore*Rate), 
                                 int(SecondsAfter*Rate), 
                                 AnalogTTLs=True)
    X = np.arange(-int(SecondsBefore*Rate), 
                  int(SecondsAfter*Rate))*1000/Rate
    
    return(ABR, X)


def PlotCheck():
    if 'plt' not in globals():
        Params = Plot.Set(Params=True)
        from matplotlib import rcParams; rcParams.update(Params)
        from matplotlib import pyplot as plt
    
        return(plt)
    else:
        return(globals()['plt'])


def ABRPSDPlot(ABR):
    plt = PlotCheck()
    
    ABRPSD = DataAnalysis.FilterSignal(ABR.mean(axis=1), Rate['100'], [200,10000])
    ABRPSD = DataAnalysis.PSD(ABRPSD, Rate['100'], WindowSize=ABRPSD.shape[0])
    plt.plot(ABRPSD[0], ABRPSD[1]); plt.show()
    
    return(None)


def ABRSinglePlot(ABR, X, Ax=None, Show=True):
    ReturnAx = True
    if not Ax:
        ReturnAx = False
        plt = PlotCheck()
        Fig, Ax = plt.subplots()
    
    for T in range(ABR.shape[1]): 
        Ax.plot(X, DataAnalysis.FilterSignal(ABR[:,T], Rate['100'], FilterFreq), 
                'r', alpha=0.05)
    Ax.plot(X, DataAnalysis.FilterSignal(ABR.mean(axis=1), Rate['100'], FilterFreq), 
            'k', lw=2)
    
    if ReturnAx: 
        return(Ax)
    else:
        Plot.Set(Ax=Ax, Fig=Fig)
        if Show: plt.show()
        else: plt.close()
        return(None)


def ABRMultiplePlot(ABRs, X, AxArgs={}, Show=True):
    Recs = sorted(list(ABRs.keys()))
    YAmp = 0; YLim = None
    for Rec in ABRs.values():
        if max(Rec)-min(Rec) > YAmp:
            YAmp = max(Rec)-min(Rec)
            YLim = [min(Rec), max(Rec)]
    AxArgs['ylim'] = YLim
    AxArgs['xtickspacing'] = 3
    
    plt = PlotCheck()
    Fig, Axes = plt.subplots(len(ABRs), 1, sharex=True)
    for R,Rec in enumerate(Recs):
        Axes[R].plot(X, ABRs[Rec], 'r', lw=2)
        Plot.Set(Ax=Axes[R], Fig=Fig, AxArgs=AxArgs)
        if R != len(Recs)-1:
            Axes[R].tick_params(bottom='off')
            Axes[R].spines['bottom'].set_visible(False)
        
        if R != len(Recs)//2:
            Axes[R].tick_params(left='off')
            Axes[R].spines['left'].set_visible(False)
            Axes[R].set_yticklabels([])
    
    if Show: plt.show()
    else: plt.close()
    
    return(None)


#%% All sessions from a group
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


#%% Multiple folders - Single ABR session
Folder = '/home/cerebro/Malfatti/Data/ToBeAssigned/A1'
InfoFile = '/home/cerebro/Malfatti/Data/ToBeAssigned/A1/20180510095719-A1-Sound.dict'
AnalysisFile = 'Test.hdf5'
ABRCh = [1]
TTLCh = 0
StimType = ['Sound']

Folders = glob(Folder + '/????-*'); Folders.sort()
AnalysisPath = '-'.join(InfoFile.split('/')[-1].split('-')[:2]) + '-ABRs'
ABRs.Analysis(Folders, InfoFile, AnalysisFile, StimType=StimType)
ABRPlot.Traces(AnalysisPath, AnalysisFile, InfoFile, Folder+'/Figs', Save=False, Show=True)


#%% Multiple recs - Single ABR frequency
Folder = '/home/cerebro/Malfatti/Data/ToBeAssigned/A1/2018-05-10_09-57-27_A1-1012'
Recs = 'all'
Proc = '100'

ABRCh = [1]
TTLCh = 0
SecondsBefore = 0.003; SecondsAfter = 0.012
FilterFreq = [600, 1500]

Data, Rate = IO.DataLoader(Folder, Unit='uV', ChannelMap=[])
if type(Recs) == str:
    if Recs.lower() == 'all': Recs = sorted(list(Data[Proc].keys()))

ABRs = {}
for Rec in Recs:
    ABR, X = ABRSingle(Data[Proc][Rec], Rate[Proc], ABRCh, TTLCh, SecondsBefore, SecondsAfter, FilterFreq)
    ABRs[Rec] = DataAnalysis.FilterSignal(ABR.mean(axis=1), Rate[Proc], FilterFreq)
    ABRs[Rec] *= 1000 # Convert to microvolt

ABRMultiplePlot(ABRs, X, AxArgs={'xlim': [(-SecondsBefore*1000), SecondsAfter*1000]})


#%% Single rec - Single ABR intensity
Folder = '/home/cerebro/Malfatti/Data/ToBeAssigned/A1/2018-05-10_09-57-27_A1-1012'
Rec = '0'
Proc = '100'

ABRCh = [1]
TTLCh = 0
SecondsBefore = 0.003; SecondsAfter = 0.012
FilterFreq = [600, 1500]

Data, Rate = IO.DataLoader(Folder, Unit='uV', ChannelMap=[])
ABR, X = ABRSingle(Data[Proc][Rec], Rate[Proc], ABRCh, TTLCh, SecondsBefore, SecondsAfter, FilterFreq)
ABRSinglePlot(ABR*1000, X)
# ABRPSDPlot(ABR)


#%% Convert hdf5 info to dict
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
