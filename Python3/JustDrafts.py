# -*- coding: utf-8 -*-
"""
Just drafts
"""
#%% Treadmill
import DataAnalysis, Hdf5F
import numpy as np

from glob import glob

Params = DataAnalysis.Plot.SetPlot(Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt

SensorCh = 17
Diameter = 60; PeaksPerCycle = 12

Paths = glob('Treadmilltina/2017-03-08*')
DataPath = Paths[0]
Data = Hdf5F.LoadOEKwik(DataPath)
Proc = Hdf5F.GetProc(Data, 'OE')

PeakDist = (np.pi*Diameter)/PeaksPerCycle

Recs = Data[Proc]['data'].keys()
for Rec in Recs:
    RecData = Data[Proc]['data'][Rec]
    Rate = Data[Proc]['info'][Rec]['sample_rate']
    SensorData = RecData[:,SensorCh]*-1
    
    Peaks = DataAnalysis.QuantifyTTLsPerRec(True, SensorData)
    
    V = [float('NaN')]*Peaks[0]
    for P in range(1, len(Peaks)):
        Samples = Peaks[P] - Peaks[P-1]; Time = Samples/Rate
        Speed = PeakDist/Time; V = V + [Speed]*Samples
    
    Ch = 7
    Pxx, F, B, I = plt.specgram(RecData[:Ch], 200, Rate, noverlap=100)
    
    
    


#%% RT plots
import KwikAnalysis

Params = KwikAnalysis.SetPlot(Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt

Files = ['Hist-NonRT-T' + str(_) for _ in range (8)]
Jitter = {}
         
for File in Files:
    with open(File, 'r') as FHist:
        for LineNo, Line in enumerate(FHist, 1):
            pass
    
    Jitter[File] = [[] for _ in range(LineNo)]
    with open(File, 'r') as FHist:
        for Ind, Line in enumerate(FHist):
            if Line == 'Avg:\n': Line = '51'
            Jitter[File][Ind] = int(Line)


plt.plot(Jitter['Hist-NonRT-T0'], 'bo'); plt.show()

#%% Nernst potential
Temperature = 25    # in Celsius
Ions = ['Sodium', 'Potassium', 'Chloride']

Sodium = {'In':8,
          'Out': 130 + 20,
          'z': +1}

Potassium = {'In': 108,
             'Out': 3 + 1.2,
             'z': +1}

Chloride = {'In': 9,
            'Out': 130 + 3 + 4.8,
            'z':-1}


import numpy as np

e = 1.6021766 * 10**-19
Na = 6.022141 * 10**23

F = e*Na 
R = 8.314472
T = Temperature + 274.150

Nernst = {}
for Ion in Ions:
    IonDict = eval(Ion)
    Nernst[Ion] = (((R*T)/(IonDict['z']*F)) * 
                   np.log(IonDict['Out']/IonDict['In']) * 1000) # in mV


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


#%% GPIASDrafts
import Hdf5F
Hdf5F.DeleteGroup('GPIAS', 'CaMKIIahM4Dn06/CaMKIIahM4Dn06-Analysis.hdf5', All=True)
Hdf5F.DeleteGroup('GPIAS', 'CaMKIIahM4Dn07/CaMKIIahM4Dn07-Analysis.hdf5', All=True)
Hdf5F.DeleteGroup('GPIAS', 'CaMKIIahM4Dn08/CaMKIIahM4Dn08-Analysis.hdf5', All=True)
Hdf5F.DeleteGroup('GPIAS', 'CaMKIIahM4Dn09/CaMKIIahM4Dn09-Analysis.hdf5', All=True)


#%% GPIASIndex
import h5py
import KwikAnalysis
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
import DataAnalysis, Hdf5F, KwikAnalysis
import numpy as np
from glob import glob
from itertools import permutations

Animals = ['CaMKIIahM4Dn06', 'CaMKIIahM4Dn07', 'CaMKIIahM4Dn08', 'CaMKIIahM4Dn09']

Params = {'backend': 'TkAgg'}
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt

ABRThresholds = {Animal: {} for Animal in Animals}
for Animal in Animals:
    AnalysisFile = glob(Animal + '/*.hdf5')[0]
    Exps = Hdf5F.GetGroupKeys('/ABRThresholds', AnalysisFile)
    
    for Exp in Exps:
        ABRThresholds[Animal][Exp] = Hdf5F.LoadDict('/ABRThresholds/'+Exp, 
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
Pairs = {}
for Pair in permutations(''.join(DataAnalysis.StrRange('0', str(len(Exps)))), 2):
    PKey = min(Pair)+max(Pair)
    if PKey in Pairs: continue
    
    Pairs[PKey] = {}
    
    for FInd, Freq in enumerate(Freqs):
#        DataA = ThresholdsPerFreq[int(Pair[0])][FInd]
#        DataB = ThresholdsPerFreq[int(Pair[1])][FInd]
        DataA = TPFExpanded[int(Pair[0])][FInd]
        DataB = TPFExpanded[int(Pair[1])][FInd]
        CL = 1 - (0.05/6)
        
        Pairs[PKey][Freq] = DataAnalysis.Stats.RTTest(DataA, DataB, Confidence=CL)
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
        DataAnalysis.Plot.SignificanceBar(XInd[0], XInd[1], Counter, Text, Axes[1], TicksDir='up')
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
Axes[1].set_ylim(-2, 7)
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

ABRCh = [10]         # [RightChannel, LeftChannel], if order matters
ABRTTLCh = 17            # TTL ch for ABR
ABRTimeBeforeTTL = 3    # in ms
ABRTimeAfterTTL = 9    # in ms
FilterFreq = [300, 3000]         # frequeutter order
AnalogTTLs = True
FilterOrder = 5         # b
StimType = ['Sound_NaCl', 'Sound_CNO']
#StimType = ['Sound']
Board = 'OE'
Override = {}

import Hdf5F
import numpy as np
import os
import struct
import tarfile
import wave
from datetime import datetime
from glob import glob
from scipy import signal
from scipy.io import wavfile

FileName = glob('*.hdf5')[0]

def GetProc(Raw, Board):
    print('Get proc no. for', Board, 'board... ', end='')
    ProcChs = {Proc: len(Raw[Proc]['data']['0'][1,:]) 
               for Proc in Raw.keys()}
    
    for Proc, Chs in ProcChs.items():
        if Chs == max(ProcChs.values()): OEProc = Proc
        else: RHAProc = Proc
    
    if 'RHAProc' not in locals(): RHAProc = OEProc
    
    if Board == 'OE': Proc = OEProc
    elif Board == 'RHA': Proc = RHAProc
    else: print("Choose Board as 'OE' or 'RHA'."); return(None)
    
    print('Done.')
    return(OEProc, RHAProc, Proc)


def FilterSignal(Signal, Rate, Frequency, FilterOrder=4, Type='bandpass'):
    if Type not in ['bandpass', 'lowpass', 'highpass']:
        print("Choose 'bandpass', 'lowpass' or 'highpass'.")
    
    elif len(Frequency) not in [1, 2]:
        print('Frequency must have 2 elements for bandpass; or 1 element for \
        lowpass or highpass.')
    
    else:
        passband = [_/(Rate/2) for _ in Frequency]
        f2, f1 = signal.butter(FilterOrder, passband, Type)
        Signal = signal.filtfilt(f2, f1, Signal, padtype='odd', padlen=0)
        
        return(Signal)


def GetRecKeys(Raw, Events, AnalogTTLs):
    for Processor in Raw.keys():
        if '0' not in list(Raw[Processor]['data'].keys()):
            print('Rec numbers are wrong. Fixing...')
            for iKey in Raw[Processor].keys():
                Recs = list(Raw[Processor][iKey].keys())
                Recs = [int(_) for _ in Recs]; Min = min(Recs)
                
                for Key in Recs:
                    Raw[Processor][iKey][str(Key-Min)] = \
                        Raw[Processor][iKey].pop(str(Key))
                
                if AnalogTTLs:
                    return(Raw)
                else:
                    EventRec = Events['TTLs']['recording'][:]
                    for _ in range(len(EventRec)): 
                        EventRec[_] = EventRec[_] - Min
                    return(Raw, EventRec)
            
            print('Fixed.')
    
    else:
        if AnalogTTLs: return(Raw)
        else: EventRec = Events['TTLs']['recording']; return(Raw, EventRec)


def GetTTLInfo(Events, EventRec, TTLCh):
    print('Get TTL data...')
    EventID = Events['TTLs']['user_data']['eventID']
    EventCh = Events['TTLs']['user_data']['event_channels']
    EventSample = Events['TTLs']['time_samples']

#    TTLChs = np.nonzero(np.bincount(EventCh))[0]
    TTLRecs = np.nonzero(np.bincount(EventRec))[0]
    TTLRecs = ["{0:02d}".format(_) for _ in TTLRecs]
    TTLsPerRec = {Rec: [EventSample[_] for _ in range(len(EventRec)) 
                         if EventRec[_] == int(Rec)
                         and EventCh[_] == TTLCh-1 
                         and EventID[_] == 1]
                  for Rec in TTLRecs}
#    TTLRising = Kwik.get_rising_edge_times(Files['kwe'], TTLCh-1)
    
    return(TTLsPerRec)


def QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, ChTTL=-1, Proc='', TTLsPerRec=[], 
                       Rate=[]):
    print('Get TTL timestamps... ', end='')
    if AnalogTTLs:
        TTLCh = Raw[Proc]['data'][Rec][:, ChTTL-1]
        Threshold = max(TTLCh)/2
        TTLs = []
        for _ in range(1, len(TTLCh)):
            if TTLCh[_] > Threshold:
                if TTLCh[_-1] < Threshold: TTLs.append(_)
        
        print('Done.')
        return(TTLs)
    else:
#        TTLNo = [0]
#        for _ in range(1, len(TTLsPerRec)+1):
#            TTLNo = TTLNo + [len(TTLsPerRec[_-1]) + TTLNo[-1]]
#        TTLNo = [0] + [TTLNo[_]-1 for _ in range(1, len(TTLNo))]
#        
#        if Rec == 0:
#            sTTLNo = 0
#        else:
#            sTTLNo = TTLNo[Rec] + 1
#        
#        return(TTLNo, sTTLNo)
        RawTime = [_*Rate for _ in Raw['timestamps'][int(Rec)]]
        TTLs = TTLsPerRec[Rec]
        
        print('Done.')
        return(RawTime, TTLs)


def FixTTLs(Array, TTLsToFix):
    for TTL in TTLsToFix:
        nInd = np.random.randint(1, 100)
        while nInd == TTL: nInd = np.random.randint(0, 100)
        
        print('TTL', str(TTL), 'was replaced by', str(nInd))
        Array[TTL] = Array[nInd]
    
    return(Array)


def SliceData(Data, Proc, Rec, TTLs, DataCh, NoOfSamplesBefore, 
              NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime=[]):
    print('Slicing data around TTL...')
    Array = [[0 for _ in range(NoOfSamples)] for _ in range(len(TTLs))]
    TTLsToFix = []
    
    for TTL in range(len(TTLs)):
        if AnalogTTLs: TTLLoc = int(TTLs[TTL])
        else: TTLLoc = int(RawTime.index(TTLs[TTL]))#)/Rate)
        
        Start = TTLLoc-NoOfSamplesBefore
        End = TTLLoc+NoOfSamplesAfter
        
        if Start < 0: Start = 0; End = End+(Start*-1); TTLsToFix.append(TTL)
        
        Array[TTL] = Data[Proc]['data'][Rec][Start:End, DataCh[0]-1] * \
                     Data[Proc]['channel_bit_volts'][Rec][DataCh[0]-1] # in mV
        
        if len(Array[TTL]) != End-Start: TTLsToFix.append(TTL)
            
    Array = FixTTLs(Array, TTLsToFix)
    
    print('Done.')
    return(Array)


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
        wavfile.write(FileName, Rate, ABRs[Key])
#        WriteWav(Data, 1, DataSize, Rate, FileName)
    
    XValues = XValues.tobytes()
    FileName = Path.split(sep='/'); del(FileName[0])
    FileName = 'XValues-' + '_'.join(FileName)
    
    WriteWav(XValues, 1, len(XValues), Rate, FileName)


print('Load DataInfo...')
DirList = glob('KwikFiles/*'); DirList.sort()
DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)

AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
Now = datetime.now().strftime("%Y%m%d%H%M%S")
Here = os.getcwd().split(sep='/')[-1]
Group = Here + '-ABRs_' + Now
    
for Stim in StimType:
    if Override != {}: 
        if 'Stim' in Override.keys(): Stim = Override['Stim']
    
    Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
    
    for RecFolder in Exps:
        ABRs = {}; Info = {}
        ExpInfo = Hdf5F.ExpExpInfo(RecFolder, DirList, FileName)
        
        if AnalogTTLs: 
            Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
        else: 
            Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
        
        OEProc, RHAProc, ABRProc = GetProc(Raw, Board)
        
        if AnalogTTLs: Raw = GetRecKeys(Raw, [0], AnalogTTLs)
        else:
            Raw, EventRec = GetRecKeys(Raw, Events, AnalogTTLs)
            TTLsPerRec = GetTTLInfo(Events, EventRec, ABRTTLCh)
        
        Rate = Raw[OEProc]['info']['0']['sample_rate']
        NoOfSamplesBefore = ABRTimeBeforeTTL*int(Rate*10**-3)
        NoOfSamplesAfter = ABRTimeAfterTTL*int(Rate*10**-3)
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        
        SFreq = ''.join([str(DataInfo['NoiseFrequency'][ExpInfo['Hz']][0]),
                         '-',
                         str(DataInfo['NoiseFrequency'][ExpInfo['Hz']][1])])
        
        XValues = (range(-NoOfSamplesBefore, 
                                 NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
        
        for Rec in Raw[OEProc]['data'].keys():
            print('Slicing and filtering ABRs Rec ', str(Rec), '...')
            
            if AnalogTTLs:
                TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, ABRTTLCh, 
                                          OEProc)
                ABR = SliceData(Raw, ABRProc, Rec, TTLs, ABRCh, 
                                NoOfSamplesBefore, NoOfSamplesAfter, 
                                NoOfSamples, AnalogTTLs)
            else:
                RawTime, TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, 
                                                   TTLsPerRec=TTLsPerRec)
                ABR = SliceData(Raw, ABRProc, Rec, TTLs, ABRCh, 
                                NoOfSamplesBefore, NoOfSamplesAfter, 
                                AnalogTTLs, RawTime)
            
            for TTL in range(len(TTLs)):
                ABR[TTL] = FilterSignal(ABR[TTL], Rate, [min(FilterFreq)], 
                                        FilterOrder, 'highpass')
            
            ABR = np.mean(ABR, axis=0)
            ABR = FilterSignal(ABR, Rate, [max(FilterFreq)], FilterOrder, 
                               'lowpass')
            
            dB = str(DataInfo['Intensities'][int(Rec)]) + 'dB'
            ABRs[dB] = ABR[:]; del(ABR)
        
        Path = Stim + '/' + ExpInfo['DVCoord'] + '/' + SFreq
        Hdf5F.WriteABR(ABRs, Info['XValues'], Group, Path, AnalysisFile)


#%% Fix old ExpInfo files
import h5py

for ind, key in enumerate(F['ExpInfo'].keys()):
    if ind <5: F['ExpInfo'][key].attrs['StimType'] = np.string_(['Sound_NaCl'])
    else: F['ExpInfo'][key].attrs['StimType'] = np.string_(['Sound_CNO'])




#%% SavGol
        # SavGol smooth
#        GPIAS[Freq]['NoGap'] = signal.savgol_filter(GPIAS[Freq]['NoGap'], 5, 2, 
#                                                    mode='nearest')
#        GPIAS[Freq]['Gap'] = signal.savgol_filter(GPIAS[Freq]['Gap'], 5, 2, 
#                                                  mode='nearest')
        

#%% RenameFiles
import os
from glob import glob

ModifyFreq = ['20160217-GPIASn01Screening']

for Folder in ModifyFreq:
    Files = glob(Folder + '/*.mat')
    
    for File in Files:
        NewFile = File.replace('  ', '_')
        os.rename(File, NewFile)


## 
ModifyName = glob('*GPIAS*'); ModifyName.sort()

for Folder in ModifyName:
    Files = glob(Folder + '/*.mat')
    Start = len(Folder) + 1
        
    for File in Files:
        NewFile = File[Start:Start+14] + '_' + File[Start+15:]
        NewFile = NewFile.split('_')
        Date = NewFile[0]; Animal = NewFile[1]
        DateIndex = NewFile.index(Date); AnimalIndex = NewFile.index(Animal)
        NewFile[0], NewFile[1] = NewFile[AnimalIndex], NewFile[DateIndex]
        NewFile = File[:Start] + '_'.join(NewFile)
        os.rename(File, NewFile)
#        print(File + '\n', NewFile)

##
ModifyFreq = ['20160217-GPIASn01Screening']

for Folder in ModifyFreq:
    Files = glob(Folder + '/*.mat')
    Start = len(Folder) + 1
    
    for File in Files:
        if File[Start:Start+4] == '1OD-': NewFile = File[:Start] + '1R_' + File[Start+4:]
        elif File[Start:Start+4] == '1OE-': NewFile = File[:Start] + '1L_' + File[Start+4:]
        elif File[Start:Start+4] == '2OD-': NewFile = File[:Start] + '2R_' + File[Start+4:]
        elif File[Start:Start+4] == '2OE-': NewFile = File[:Start] + '2L_' + File[Start+4:]
        elif File[Start:Start+5] == '1ODE-': NewFile = File[:Start] + '1LR_' + File[Start+5:]
        
        os.rename(File, NewFile)

#%% Freq sort

FKeys = list(Sound.keys()); ToPrepend = []
for FF in ['8000-10000', '9000-11000']:
    if FF in FKeys:
        del(FKeys[FKeys.index(FF)]); ToPrepend.append(FF)

ToPrepend.sort(); FKeys = ToPrepend + FKeys
