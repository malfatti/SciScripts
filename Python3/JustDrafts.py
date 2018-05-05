# -*- coding: utf-8 -*-
"""
Just drafts :)

@author: T. Malfatti <malfatti@disroot.org>
@year: 2017
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
#%% RT plots
from DataAnalysis.Plot import Plot
from glob import glob

Files = sorted(glob('Tests/RTTest/Hist-NonRT-*'))
Jitter = Plot.RTTest(Files, Return=True)


#%% TTLs figure

import ControlSoundBoard
import KwikAnalysis
import numpy as np

AnimalName = 'CaMKIIahM4Dn04'
Rate = 128000
BaudRate = 38400

CalibrationFile = '20160419093139-SoundMeasurement.hdf5'
SoundBoard = 'USBPre2_oAux-iAux'

TTLAmpF = 1
SoundTTLVal = 0.6
LaserTTLVal = 4.1

SoundPrePauseDur = 0.008
SoundPulseDur = 0.003
SoundPostPauseDur = 0.039
SoundPulseNo = 2
SoundStimBlockNo = 1
SoundPauseBetweenStimBlocksDur = 1
Intensities = [80]
NoiseFrequency = [[8000, 10000]]
SBOutAmpF = 1

LaserPrePauseDur = 0.004
LaserPulseDur = 0.01
LaserPostPauseDur = 0.036
LaserPulseNo = 2
LaserStimBlockNo = 1
LaserPauseBetweenStimBlocksDur = 1

SoundAmpF = ControlSoundBoard.dBToAmpF(Intensities, CalibrationFile)
SoundPulse = ControlSoundBoard.GenNoise(Rate, SoundPulseDur)
SoundPulseFiltered = ControlSoundBoard.BandpassFilterSound(SoundPulse, Rate, NoiseFrequency)
SoundUnit = ControlSoundBoard.ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                               NoiseFrequency, SBOutAmpF, SoundPrePauseDur, 
                               SoundPostPauseDur)
SoundTTLUnit = ControlSoundBoard.GenTTL(Rate, SoundPulseDur, TTLAmpF, 
                                        SoundTTLVal, SoundBoard, SBOutAmpF, 
                                        SoundPrePauseDur, SoundPostPauseDur)

LaserUnit = ControlSoundBoard.GenTTL(Rate, LaserPulseDur, TTLAmpF, LaserTTLVal, SoundBoard, SBOutAmpF, 
                       LaserPrePauseDur, LaserPostPauseDur)

SoundTTLLaserUnit = [LaserUnit[_]+SoundTTLUnit[_] \
                     for _ in range(len(SoundTTLUnit))]

SoundTTL = []
for El in SoundTTLUnit:
    if El > 0: SoundTTL.append(El)
    else: SoundTTL.append(0)

LaserTTL = []
for El in LaserUnit:
    if El > 0: LaserTTL.append(El)
    else: LaserTTL.append(0)

SoundLaserTTL = []
for El in SoundTTLLaserUnit:
    if El > 0: SoundLaserTTL.append(El)
    else: SoundLaserTTL.append(0)

XValues = (range(round(0.05*Rate))/np.array([Rate]))*10**3
    
Signals = {}
Signals['SoundPulse'] = {'X': XValues,
                         'Y': np.array(SoundUnit[0][0])*1000,
                         'Label': 'Sound pulse',
                         'XLabel': 'Time [ms]',
                         'YLabel': 'Voltage [mV]',
                         'Name': 'SoundPulse'}

Signals['SoundSqWave'] = {'X': XValues,
                          'Y': SoundTTLUnit,
                          'Label': 'Sound square wave',
                          'XLabel': 'Time [ms]',
                          'YLabel': 'Voltage [V]',
                          'Name': 'SoundSqWave'}

Signals['SoundTTL'] = {'X': XValues,
                       'Y': SoundTTL,
                       'Label': 'Sound TTL',
                       'XLabel': 'Time [ms]',
                       'YLabel': 'Voltage [V]',
                       'Name': 'SoundTTL'}

Signals['LaserSqWave'] = {'X': XValues,
                          'Y': LaserUnit,
                          'Label': 'Laser square wave',
                          'XLabel': 'Time [ms]',
                          'YLabel': 'Voltage [V]',
                          'Name': 'LaserSqWave'}

Signals['LaserTTL'] = {'X': XValues,
                       'Y': LaserTTL,
                       'Label': 'LaserTTL',
                       'XLabel': 'Time [ms]',
                       'YLabel': 'Voltage [V]',
                       'Name': 'LaserTTL'}

Signals['SoundLaserSqWave'] = {'X': XValues,
                               'Y': SoundTTLLaserUnit,
                               'Label': 'Sound + laser square wave',
                               'XLabel': 'Time [ms]',
                               'YLabel': 'Voltage [V]',
                               'Name': 'SoundLaserSqWave'}

Signals['SoundLaserTTL'] = {'X': XValues,
                            'Y': SoundLaserTTL,
                            'Label': 'Sound + laser TTL',
                            'XLabel': 'Time [ms]',
                            'YLabel': 'Voltage [V]',
                            'Name': 'SoundLaserTTL'}



def GeneralPlot(Data, Visible=True):
    Params = KwikAnalysis.SetPlot(Backend='TkAgg', Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(4,2))
    plt.plot(Data['X'], Data['Y'], label=Data['Label'])
    plt.ylabel(Data['YLabel']); plt.xlabel(Data['XLabel'])
    plt.legend(loc='lower right')
    plt.locator_params(tight=True)
    plt.tick_params(direction='out')
    plt.axes().spines['right'].set_visible(False)
    plt.axes().spines['top'].set_visible(False)
    plt.axes().yaxis.set_ticks_position('left')
    plt.axes().xaxis.set_ticks_position('bottom')
    plt.savefig(Data['Name']+'.svg', format='svg')
    if Visible: plt.show()
    
    return(None)


for Signal in Signals:
    GeneralPlot(Signals[Signal])


#%% GPIASIndex
import h5py
from glob import glob

Params = KwikAnalysis.SetPlot(Backend='TkAgg', Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt

Animals = ['CaMKIIahM4Dn06', 'CaMKIIahM4Dn07', 'CaMKIIahM4Dn08', 
           'CaMKIIahM4Dn09']

GPIASIndex = {Animal: {} for Animal in Animals}

for Animal in Animals:
    AnalysisFile = glob(Animal+'/*.hdf5')[0]
    
    with h5py.File(AnalysisFile, 'r') as F:
        Keys = [Key for Key in F.keys() if 'GPIAS_' in Key]; Keys.sort()
        
        for Key in Keys:
            Freqs = list(F[Key]['GPIAS']); ToPrepend = []
            for FF in ['8000-10000', '9000-11000']:
                if FF in Freqs:
                    del(Freqs[Freqs.index(FF)]); ToPrepend.append(FF)
            
            ToPrepend.sort(); Freqs = ToPrepend + Freqs
            GPIASIndex[Animal][Key] = [F[Key]['GPIAS'][Freq]['GPIASIndex'][()]
                                       for Freq in Freqs]

# Override
del(GPIASIndex['CaMKIIahM4Dn06']['CaMKIIahM4Dn06-20160623-GPIAS_2016-06-23_14-27-40-GPIAS_20160831160855'])
del(GPIASIndex['CaMKIIahM4Dn06']['CaMKIIahM4Dn06-20160630-GPIAS_2016-06-30_10-32-42-GPIAS_20160831160947'])
del(GPIASIndex['CaMKIIahM4Dn07']['CaMKIIahM4Dn07-20160217-GPIASn01Screening_2016-02-17-GPIAS_20160831161040'])
del(GPIASIndex['CaMKIIahM4Dn08']['CaMKIIahM4Dn08-20160629-GPIAS_2016-06-29_14-27-21-GPIAS_20160831161134'])
del(GPIASIndex['CaMKIIahM4Dn08']['CaMKIIahM4Dn08-20160630-GPIAS_2016-06-30_11-13-25-GPIAS_20160831161200'])
del(GPIASIndex['CaMKIIahM4Dn09']['CaMKIIahM4Dn09-20160629-GPIAS_2016-06-29_15-17-25-GPIAS_20160831161347'])
del(GPIASIndex['CaMKIIahM4Dn09']['CaMKIIahM4Dn09-20160630-GPIAS_2016-06-30_11-53-13-GPIAS_20160831161413'])

for Animal in Animals:
    for Exp in GPIASIndex[Animal]:
        if len(GPIASIndex[Animal][Exp]) != 5:
            GPIASIndex[Animal][Exp].insert(1, float('NaN'))

YLabel = 'GPIAS index'
Colors = ['-ro', '-bx', '-g^', '-ms']
Freqs = ['8$\sim$10', '9$\sim$11', '10$\sim$12', '12$\sim$14', 
         '14$\sim$16 [KHz]']
When = ['Before ANT', '2d after ANT', '1w after ANT (NaCl)', 
             '1w after ANT (CNO)']

Fig, Axes = plt.subplots(1,len(When), sharey=True)

for AInd, Animal in enumerate(Animals):
    Exps = list(GPIASIndex[Animal]); Exps.sort()
    
    for EInd, Exp in enumerate(Exps):
        Axes[EInd].plot(GPIASIndex[Animal][Exp], Colors[AInd], label=Animal)
        KwikAnalysis.SetPlot(AxesObj=Axes[EInd], Axes=True)
        Axes[EInd].xaxis.set_ticks(range(len(Freqs)))
        Axes[EInd].set_xticklabels(Freqs)
        Axes[EInd].set_xlabel(When[EInd])
        Axes[EInd].legend(loc='upper left')
        Axes[EInd].spines['bottom'].set_position(('outward', 5))
        if EInd == 0: Axes[EInd].spines['left'].set_position(('outward', 5))
        else: 
            Axes[EInd].spines['left'].set_visible(False)
            Axes[EInd].yaxis.set_ticks_position('none')
    
    Axes[0].set_ylabel(YLabel)
plt.show()


#%% Plot ABRThresholds
import numpy as np

from DataAnalysis import Stats
from DataAnalysis.DataAnalysis import StrRange
from DataAnalysis.Plot import Plot
from glob import glob
from IO import Hdf5
from itertools import combinations
#%%
Animals = ['CaMKIIahM4Dn06', 'CaMKIIahM4Dn07', 'CaMKIIahM4Dn08', 'CaMKIIahM4Dn09']

Params = Plot.Set(Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt

ABRThresholds = {Animal: {} for Animal in Animals}
for Animal in Animals:
    AnalysisFile = glob(Animal + '*.hdf5')[0]
    Exps = Hdf5.GetGroupKeys('/ABRThresholds', AnalysisFile)
    
    for Exp in Exps:
        ABRThresholds[Animal][Exp] = Hdf5.DictLoad('/ABRThresholds/'+Exp, 
                                                    AnalysisFile, Attrs=False)

Freqs = [Freq for Freq in ABRThresholds[Animal][Exp]]
Freqs = sorted(Freqs)

ThresholdsPerFreq = [[[] for Freq in Freqs] for Exp in range(len(Exps))]
for Animal in Animals:
    Exps = [Exp for Exp in ABRThresholds[Animal]]; Exps.sort()
    
    for EInd, Exp in enumerate(Exps):
        for FInd, Freq in enumerate(Freqs):
            ThresholdsPerFreq[EInd][FInd].append(ABRThresholds[Animal][Exp][Freq])

TPFExpanded = ThresholdsPerFreq[:]
for EInd in range(2,4):
    TPFExpanded[EInd] = [TPFExpanded[EInd][FInd] + TPFExpanded[EInd][FInd]
                         for FInd in range(len(Freqs))]

# For t-tests
Pairs = {}; Combinations = list(combinations(''.join(StrRange('0', str(len(Exps)))), 2))
for Pair in Combinations:
    PKey = min(Pair)+max(Pair)
    if PKey in Pairs: continue
    
    Pairs[PKey] = {}
    
    for FInd, Freq in enumerate(Freqs):
#        DataA = ThresholdsPerFreq[int(Pair[0])][FInd]
#        DataB = ThresholdsPerFreq[int(Pair[1])][FInd]
        DataA = TPFExpanded[int(Pair[0])][FInd]
        DataB = TPFExpanded[int(Pair[1])][FInd]
        CL = 1 - 0.05
        
        if np.mean(DataA) > np.mean(DataB): DataA, DataB = DataB, DataA
        Pairs[PKey][Freq] = Stats.RTTest(DataA, DataB, Confidence=CL)
        Pairs[PKey][Freq]['p.value'] *= len(Combinations)
        print('Pair', str(Pair), 'Freq', Freq + ':', 
              str(Pairs[PKey][Freq]['p.value']))

## For Anova
#Data = {'Freqs': [], 'Values':[]}
#for Freq in range(len(Freqs)):
#    for Exp in range(len(Exps)):
#        for Value in TPFExpanded[Exp][Freq]:
#            Data['Freqs'].append(Freq); Data['Values'] = Data['Values'] + [Value]

ToDelete = []
for Pair in Pairs:
    for Freq in Freqs:
        if Pairs[Pair][Freq]['p.value'] > 0.05:
            ToDelete.append([Pair, Freq])
        
        if Pairs[Pair][Freq]['p.value'] != Pairs[Pair][Freq]['p.value']:
            ToDelete.append([Pair, Freq])
        

for KeyPair in ToDelete:
    del(Pairs[KeyPair[0]][KeyPair[1]])

EmptyPairs = []
for Pair in Pairs:
    if len(Pairs[Pair]) == 0:
        EmptyPairs.append(Pair)

for Pair in EmptyPairs:
    del(Pairs[Pair])

# Plots
def DiffLine(XInd, Y, Text, Axes, lw=1):
    X = XInd[0] + ((max(XInd) - min(XInd))/2)
    props = {'connectionstyle':'bar','arrowstyle':'-',
             'shrinkA':20,'shrinkB':20,'lw':lw}
    
    Axes.annotate(Text, xy=(X, Y+7), ha='center')#, zorder=10)
    Axes.annotate('', xy=(XInd[0], Y), xytext=(XInd[1], Y), arrowprops=props)
    
    return(None)


#Colors = ['ro', 'bx', 'g^', 'ms']
Visible = True
Colors = ['r', 'g', 'b', 'm', 'k', '#ffa500']; ErrorCfg = {'ecolor':'k'}
XLabel = 'ABR measurements'; YLabel = 'Intensity [dBSPL]'
When = ['Before ANT', '48h after ANT', '2w after ANT (NaCl)', '2w after ANT (CNO)']

Thresholds = {}
Fig, Axes = plt.subplots(2, 1, gridspec_kw={'height_ratios':[3, 1]}, sharex=True, figsize=(12, 9))
for FInd, Freq in enumerate(Freqs):
    Thresholds[Freq] = {}
    Thresholds[Freq]['Means'] = [np.mean(ThresholdsPerFreq[EInd][FInd]) 
                                 for EInd in range(len(Exps))]
    Thresholds[Freq]['SEMs'] = [np.std(ThresholdsPerFreq[EInd][FInd]) / \
                                (len(ThresholdsPerFreq[EInd][FInd]))**0.5
                                for EInd in range(len(Exps))]
    
    
    for EInd, Exp in enumerate(Exps):
        Bars = range(len(Freqs)) + np.array([7*EInd])
        if EInd == 0:
            Axes[0].bar(Bars[FInd], Thresholds[Freq]['Means'][EInd],
                     color=Colors[FInd], alpha=0.4, error_kw=ErrorCfg,
                     yerr=Thresholds[Freq]['SEMs'][EInd], label=Freq+' Hz')
        else: 
            Axes[0].bar(Bars[FInd], Thresholds[Freq]['Means'][EInd],
                     color=Colors[FInd], alpha=0.4, error_kw=ErrorCfg,
                     yerr=Thresholds[Freq]['SEMs'][EInd])

ExpWidth = len(Freqs)+2
Ticks = [x/10 for x in range(int(len(Freqs)*10/2), 
                             int(ExpWidth*10*len(Exps)), 
                             ExpWidth*10)]

Counter = 1
for Pair in Pairs:
    for FInd, Freq in enumerate(Freqs):
        if Freq not in Pairs[Pair]: continue
        
        E0Ind = int(Pair[0]); E1Ind = int(Pair[1])
        XInd = [(E0Ind*ExpWidth)+FInd+0.5, (E1Ind*ExpWidth)+FInd+0.5]
        Y = 1.1* max(Thresholds[Freq]['Means'][E0Ind], 
                        Thresholds[Freq]['Means'][E1Ind])
        
        Text = 'p = ' + str(round(Pairs[Pair][Freq]['p.value'], 4))
#        DiffLine(XInd, Y, Text, Axes[1])
#        print(Counter)
        Plot.SignificanceBar(XInd, [Y]*2, Text, Axes[1], TicksDir='up')
        Counter += 1

for Ind in [0,1]:
    Axes[Ind].spines['top'].set_visible(False)
    Axes[Ind].spines['bottom'].set_visible(False)
    Axes[Ind].spines['right'].set_visible(False)
    Axes[Ind].yaxis.set_ticks_position('left')
    Axes[Ind].xaxis.set_ticks_position('none')
    Axes[Ind].spines['left'].set_position(('outward', 5))


Axes[0].set_ylabel(YLabel)
Axes[0].set_ylim(bottom=30)
Axes[0].legend(loc='best', frameon=False)

Axes[1].xaxis.set_ticks(Ticks)
Axes[1].set_xticklabels(When)
Axes[1].set_xlabel(XLabel)
#Axes[1].set_ylim(-2, 7)
Axes[1].spines['left'].set_visible(False)
Axes[1].yaxis.set_ticks_position('none')

Fig.savefig('ABRThresholds-' + Animals[0][:-2] + 
            '_'.join([_[-2:] for _ in Animals]) + '.svg', format='svg')
if Visible: plt.show()



#Fig, Axes = plt.subplots(1,4, sharey=True, figsize=(12, 6))
#for EInd, Exp in enumerate(Exps):
#    Maxs, Mins, Means = [], [], []
#    for FInd, Freq in enumerate(Freqs):
#        Mean = np.mean(ThresholdsPerFreq[EInd][FInd])
#        SEM = np.std(ThresholdsPerFreq[EInd][FInd]) / \
#              (len(ThresholdsPerFreq[EInd][FInd]))**0.5
#        
#        Maxs.append(Mean+SEM); Means.append(Mean); Mins.append(Mean-SEM)
#    
#    Axes[EInd].plot(Means, label='Mean')
#    Axes[EInd].fill_between(range(len(Maxs)), Maxs, Mins, color='0.8', label='SEM')
#    Axes[EInd].xaxis.set_ticks(range(len(Freqs)))
#    Axes[EInd].set_xticklabels(Freqs)
#    Axes[EInd].set_xlabel(When[EInd])
#    Axes[EInd].legend(loc='upper left')
#    Axes[EInd].spines['bottom'].set_position(('outward', 5))
#    if EInd == 0: Axes[EInd].spines['left'].set_position(('outward', 5))
#    else: 
#        Axes[EInd].spines['left'].set_visible(False)
#        Axes[EInd].yaxis.set_ticks_position('none')
#
#Axes[0].set_ylabel(YLabel)
#plt.show()
#

## Individual data
#Fig, Axes = plt.subplots(1,4, sharey=True)
#for AInd, Animal in enumerate(Animals):
#    Exps = [Exp for Exp in ABRThresholds[Animal]]; Exps.sort()
#    
#    for EInd, Exp in enumerate(Exps):        
#        Thresholds = [ABRThresholds[Animal][Exp][Freq] for Freq in Freqs]
#        
#        Axes[EInd].plot(Thresholds, Colors[AInd], label=Animal)
#        KwikAnalysis.SetPlot(AxesObj=Axes[EInd], Axes=True)
#        Axes[EInd].xaxis.set_ticks(range(len(ABRThresholds[Animal][Exp])))
#        Axes[EInd].set_xticklabels(Freqs)
#        Axes[EInd].set_xlabel(When[EInd])
#        Axes[EInd].legend(loc='upper left')
#        Axes[EInd].spines['bottom'].set_position(('outward', 5))
#        if EInd == 0: Axes[EInd].spines['left'].set_position(('outward', 5))
#        else: 
#            Axes[EInd].spines['left'].set_visible(False)
#            Axes[EInd].yaxis.set_ticks_position('none')
#    
#    Axes[0].set_ylabel(YLabel)
#plt.show()


#%% Tar+Wave revolution
import os, tarfile, wave

def WriteABRTar(ABRs, XValues, Path, FileName):
    print('Writing data to', FileName+'... ', end='')    
    FileList = glob(Path + '/**/*.*', recursive=True); FileList.sort()
    
    with tarfile.open(FileName, 'a') as F:
        for File in FileList:
            F.add(File)
    
    print('Done.')
    return(None)


def WriteWav(Data, ChNo, DataSize, Rate, FileName):
    print('Writing data to', FileName+'... ', end='')
    with wave.open(FileName, 'w') as F:
        F.setparams((ChNo, 4, Rate, DataSize*ChNo, 'NONE', 'uncompressed'))
        F.writeframes(Data)
    
    print('Done.')
    return(None)

def WriteABRWave(ABRs, XValues, Rate, Group, Path):
    Path = Group + '/' + Path
    Trial = len(glob(Path+'/*')); Trial = "{0:02d}".format(Trial)
    Path = Path + '/' + Trial
    os.makedirs(Path, exist_ok=True)
    
#    Keys = list(ABRs.keys()); Keys.sort()
#    ChNo = len(ABRs); 
    
#    Data = [ABRs[dB][Sample] for Sample in range(DataSize) for dB in Keys]
    
    for Key in ABRs:
#        Data = ABRs[Key].tobytes(); DataSize = len(Data)
        FileName = Path.split(sep='/'); del(FileName[0])
        FileName = Path + '/' + 'ABRs-' + '_'.join(FileName) + '_' + Key + '.wav'
#        wavfile.write(FileName, Rate, ABRs[Key])
#        WriteWav(Data, 1, DataSize, Rate, FileName)
    
    XValues = XValues.tobytes()
    FileName = Path.split(sep='/'); del(FileName[0])
    FileName = 'XValues-' + '_'.join(FileName)
    
    WriteWav(XValues, 1, len(XValues), Rate, FileName)


#%% hdf5 info to dict
Folders = sorted(glob('**/KlustaFiles', recursive=True))
#Folders = [_ for _ in Folders if 'Prevention' not in _]
Folders = [_ for _ in Folders if '3D' not in _]
Folders = [_ for _ in Folders if 'st_02' not in _]
Folders = [_ for _ in Folders if 'st_03' not in _]
#Folders = [_ for _ in Folders if '4Dn04' not in _]

for Folder in Folders:
    DataFolder = '/'.join(Folder.split('/')[:-1]) + '/KwikFiles'
    InfoFile = glob('/'.join(Folder.split('/')[:-1]) + '/*.dict')[0]
    
    if InfoFile[-4:] == 'hdf5':  DataInfo = Hdf5.DictLoad('/DataInfo', InfoFile)
    else: DataInfo = Txt.DictRead(InfoFile)
    
    ExpFolder = sorted(glob(DataFolder+'/*'))[0]
    Data, Rate = IO.DataLoader(ExpFolder, True, 'Bits')
    Proc = list(Data.keys())[0]
    Rec = Data[Proc]['0']; Rate = Rate[Proc]
    
    Chs = [Rec[:int(Rate/2),Ch] for Ch in range(Rec.shape[1])]
    Plot.RawCh(Chs, Lines=len(Chs), Cols=1, Save=False)
#    a, b = Hdf5.DataLoad('/DataInfo', InfoFile)
#    DataInfo = {**a, **b}
#    a, b = Hdf5.DataLoad('/ExpInfo', InfoFile)
#    DataInfo['ExpInfo'] = {**a, **b}
##    DataInfo['ExpInfo'] = Hdf5.DictLoad('/ExpInfo', InfoFile)
##    print(Txt.DictPrint(DataInfo))
#    print(InfoFile[:-4]+'dict')
#    Txt.DictWrite(InfoFile[:-4]+'dict', DataInfo)

