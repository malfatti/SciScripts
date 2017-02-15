# -*- coding: utf-8 -*-
"""
    Copyright (C) 2015  T. Malfatti
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

This is a script to generate a sound white noise at several frequencies and
intensities, play and record them at the same time, for testing purposes.
Next, it calculates the intensity in RMS and dB for each frequency at each
amplification factor. In our setup, we use this to calibrate the audio
equipment.
"""

#%% Set parameters of the experiment

Rate = 192000
# Use one that was used in SoundBoardCalibration.py
SBAmpFsFile = '/home/cerebro/Malfatti/Test/20170213142143-SBAmpFs.hdf5'
#SBAmpFsFile = '/home/malfatti/Documents/PhD/Tests/20170214093602-SBAmpFs.hdf5'
SoundBoard = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn'
Setup = 'UnitRec'

## Fill all durations in SECONDS! If

## Sound
# Pulse duration. Avoid using more than 0.5. If you need long pulses, use
# SoundPulseDur = 0.5 and SoundPulseNo = DesiredDurationInSec/0.5
SoundPulseDur = 2
# Amount of pulses per block
SoundPulseNo = 1
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
#NoiseFrequency = [[8000, 10000], [12000, 14000]]
NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], [12000, 14000], 
                  [14000, 16000], [16000, 18000], [8000, 18000]]
# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0
# Mic sensitivity, from mic datasheet, in dB re V/Pa
MicSens_dB = -47.46

#==========#==========#==========#==========#

#import array
import ControlSoundBoard
import h5py
import Hdf5F
import math
import numpy as np
import os
import pandas
import sounddevice as SD
import time
from datetime import datetime
from scipy import signal

Params = {'backend': 'Qt5Agg'}
from matplotlib import rcParams; rcParams.update(Params)
import matplotlib.pyplot as plt

SBOutAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard, 'SBOutAmpF')
SBInAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard, 'SBInAmpF')

#SoundAmpF = [2.15, 1.5, 1, 0.5, 0.1, 0]
SoundAmpF = np.hstack((
                np.arange(2.15, 1.05, -0.1), np.arange(1.0, 0.4, -0.05),
                np.arange(0.4, 0.15, -0.01), np.arange(0.15, 0.03, -0.005),
                np.arange(0.03, 0.01, -0.0005), np.arange(0.01, 0.001, -0.0001),
                np.arange(0.001, 0, -0.00002), np.array(0.0)
                ))
#SoundAmpF = FRange(2.0, 1.0, 0.1) + FRange(1.0, 0.4, 0.05) + \
#            FRange(0.4, 0.15, 0.01) + FRange(0.15, 0.03, 0.005) + \
#            FRange(0.03, 0.01, 0.0005) + FRange(0.01, 0.001, 0.0001) + \
#            FRange(0.001, 0, 0.00002) + [0.0]

## Prepare dict w/ experimental setup
Date = datetime.now().strftime("%Y%m%d%H%M%S")
Folder = ''.join([Setup, '-', Date, '-SoundMeasurement'])
FileName = Folder + '/' + Folder + '.hdf5'

DataInfo = dict((Name, eval(Name)) 
                for Name in ['Rate', 'SoundPulseDur', 'SoundPulseNo', 
                             'SoundAmpF', 'NoiseFrequency', 'TTLAmpF', 
                             'SoundBoard', 'Setup', 'MicSens_dB', 'Folder'])

## Preallocate dict for recordings
SoundRec = {}
#for FKey in Sound:
#    SoundRec[FKey] = {}

# Prepare sound pulses
Freqs = [str(Freq[0]) + '-' + str(Freq[1]) for Freq in NoiseFrequency]
SoundAmpF = {str(Freq): SoundAmpF for Freq in Freqs}

# Prepare audio objects
SD.default.device = 'system'
SD.default.samplerate = Rate
SD.default.channels = 2


for Freq in NoiseFrequency:
    SoundRec[str(Freq[0]) + '-' + str(Freq[1])] = {}
    Sound = ControlSoundBoard.SoundStim(Rate, SoundPulseDur, SoundAmpF, 
                                        [Freq], TTLAmpF, SoundBoard, 
                                        TTLs=False)[0]
    
    ## Run!
    FullTime = (len(SoundAmpF['8000-10000'])*len(NoiseFrequency)*(SoundPulseDur*SoundPulseNo))/60
#    input('Press enter to start sound measurement.')
    print('Full test will take', str(round(FullTime, 2)), 'min to run.')
    print('Current time: ', datetime.now().strftime("%H:%M:%S"))
    print('Cover your ears!!')
    print('5 ', end='')
    time.sleep(1)
    print('4 ', end='')
    time.sleep(1)
    print('3 ', end='')
    time.sleep(1)
    print('2 ', end='')
    time.sleep(1)
    print('1 ')
    time.sleep(1)
    print('Sound measurement running...')    
    for FKey in Sound:
        for AKey in Sound[FKey]:
            print(FKey, AKey)
            for Pulse in range(SoundPulseNo):
                SoundRec[FKey][AKey] = SD.playrec(Sound[FKey][AKey].T, blocking=True)
    
    print('Done playing/recording. Saving data...')
    
    ## Save!!!
    os.makedirs(Folder, exist_ok=True)
    Group = Folder
    with h5py.File(FileName) as F:
        if Group not in F: F.create_group(Group)
        if 'SoundRec' not in F[Group]: F[Group].create_group('SoundRec')
        for FKey in SoundRec:
            F[Group]['SoundRec'].create_group(FKey)
            for AKey, AVal in SoundRec[FKey].items():
                F[Group]['SoundRec'][FKey][AKey] = AVal
        
        for Key, Value in DataInfo.items():
            F[Group]['SoundRec'].attrs[Key] = Value
    
    del(Sound, SoundRec[str(Freq[0]) + '-' + str(Freq[1])])
    print('Data saved.')

print('Done \O/!')


#%% Analysis
# If needed:
#import array
#import datetime
#import glob
#import h5py
#import Hdf5F
#import KwikAnalysis
#import math
#import numpy as np
#import pandas
#from scipy import signal
#DataInfo = Hdf5F.LoadSoundMeasurement(FileName, 'DataInfo')
#SBInAmpF = Hdf5F.SoundCalibration(DataInfo['SBAmpFsFile'], SoundBoard, 'SBInAmpF')
SoundRec = Hdf5F.LoadSoundMeasurement(FileName, 'SoundRec')

DataInfo['MicSens_VPa'] = 10**(DataInfo['MicSens_dB']/20)

print('Calculating PSD, RMS and dBSLP...')
RecordingData = {}; Intensity = {}

for FKey in SoundRec:
    Intensity[FKey] = {}
    
    for AKey in SoundRec[FKey]:
        Intensity[FKey][AKey] = {}
        
        print('Saving data for', FKey, 'at', AKey)        
#        SliceStart = int(DataInfo['Rate']*0.25)-1
#        SliceEnd = SliceStart + int(DataInfo['Rate']*1)
#        RecordingData[Freq][AmpF] = RecordingData[Freq][AmpF][
#                                                        SliceStart:SliceEnd
#                                                        ]
        SoundRec[FKey][AKey] = SoundRec[FKey][AKey] * SBInAmpF
#        RecordingData[Freq][AmpF] = [_/DataInfo['SBInAmpF']
#                                     for _ in RecordingData[Freq][AmpF]]
        
        if SoundRec[FKey][AKey].shape[-1] == 2:
            Window = signal.hanning(len(SoundRec[FKey][AKey][:,0])//
                                    (DataInfo['Rate']/1000))
            F, PxxSp = signal.welch(SoundRec[FKey][AKey][:,0], 
                                    DataInfo['Rate'], Window, 
                                    nperseg=len(Window), 
                                    noverlap=0, 
                                    scaling='density')
        else:
            Window = signal.hanning(len(SoundRec[FKey][AKey])//
                                    (DataInfo['Rate']/1000))
            F, PxxSp = signal.welch(SoundRec[FKey][AKey], DataInfo['Rate'], 
                                    Window, noverlap=0, 
                                    scaling='density')
        
        FreqBand = [int(_) for _ in FKey.split('-')]
        
        Start = np.where(F > FreqBand[0])[0][0]-1
        End = np.where(F > FreqBand[1])[0][0]-1
        BinSize = F[1] - F[0]
        RMS = (sum(PxxSp[Start:End]))**0.5
#        RMS = (sum(PxxSp) * BinSize)**0.5
        
        Intensity[FKey][AKey]['PSD'] = [F, PxxSp]
        Intensity[FKey][AKey]['RMS'] = RMS
        Intensity[FKey][AKey]['dB'] = 20*(math.log(RMS/DataInfo['MicSens_VPa'], 10)) + 94
        
        del(F, PxxSp, BinSize, RMS)

del(SoundRec)

SoundIntensity = {}
for Freq in Intensity:
    SoundIntensity[Freq] = {}
    
    for AmpF in Intensity[Freq]:
        SoundIntensity[Freq][AmpF] = Intensity[Freq][AmpF]['dB']

## Save analyzed data
print('Saving analyzed data...')
os.makedirs(Folder, exist_ok=True)
Group = Folder
with h5py.File(FileName) as h5:
    if Group not in h5: h5.create_group(Group)
    h5[Group].create_group('SoundIntensity')
    
    for Freq in SoundIntensity:
        h5[Group]['SoundIntensity'].create_group(Freq)
        
        for AKey, AVal in SoundIntensity[Freq].items():
            h5[Group]['SoundIntensity'][Freq][AKey] = AVal


AmpFList = sorted(list(SoundIntensity['8000-10000'].keys()), reverse=True)
ToAppend = AmpFList[:4] + [AmpFList[-1]]
AmpFList = AmpFList[4:-1] + ToAppend
TexTable = pandas.DataFrame([[float(AmpF)] + [Intensity[Freq][AmpF]['dB'] 
                                              for Freq in Intensity]
                             for AmpF in AmpFList])

File = open(Folder+'/'+'IntensityTable.tex', 'w')
File.write(r"""
%% Configs =====
\documentclass[12pt,a4paper]{report}

\usepackage[left=0.5cm,right=0.5cm,top=0.5cm,bottom=0.5cm]{geometry}
\usepackage{longtable}
\usepackage{setspace}
\usepackage{titlesec}
\titleformat{\chapter}{\normalfont\LARGE\bfseries}{\thechapter.}{1em}{}

%% Document ======
\begin{document}

\cleardoublepage
\chapter{Sound measurements}
\input{IntensityTable-Contents.tex}

\end{document}
""")
File.close()
del(File)

File = open(Folder+'/'+'IntensityTable-Contents.tex', 'w')
File.write(TexTable.to_latex(longtable=True))
File.close()
del(File)

Colors = ['r', 'g', 'b', 'm', 'k', '#ffa500', '#00b2b2']
FreqList = list(SoundIntensity.keys()); FreqList.sort(); ToPrepend = []
for FF in ['8000-10000', '9000-11000']:
    if FF in FreqList:
        del(FreqList[FreqList.index(FF)]); ToPrepend.append(FF)

ToPrepend.sort(); FreqList = ToPrepend + FreqList

plt.figure(figsize=(6, 2))
for FKey in SoundIntensity:
    FInd = FreqList.index(FKey)
   
    FigTitle = 'Intensity curve per Freq per AmpF'
    YLabel = 'Intensity [dBSPL]'
    XLabel = 'Sound amplification factor'
    LineLabel = FKey + ' Hz'
    
    plt.semilogx(AmpFList, [SoundIntensity[FKey][AmpF] for AmpF in AmpFList], 
                 label=LineLabel, color=Colors[FInd])

plt.ylabel(YLabel); plt.xlabel(XLabel)
plt.legend(loc='lower right')
plt.locator_params(tight=True)
plt.tick_params(direction='out')
plt.axes().spines['right'].set_visible(False)
plt.axes().spines['top'].set_visible(False)
plt.axes().yaxis.set_ticks_position('left')
plt.axes().xaxis.set_ticks_position('bottom')
plt.title(FigTitle)
plt.savefig(Folder+'/'+'SoundMeasurement.svg', format='svg')


Colormaps = [plt.get_cmap(Map) for Map in ['Reds', 'Greens', 'Blues', 
                                           'Purples', 'Greys', 'Oranges', 'GnBu']]
Fig, Axes = plt.subplots(len(DataInfo['NoiseFrequency']), 
                         sharex=True, figsize=(6, 12))

for FKey in FreqList:
    Freq = FreqList.index(FKey)
    
    FigTitle = 'Power spectral density'
    AxTitle = FKey + ' Hz'
    YLabel = 'PSD [V²/Hz]'
    XLabel = 'Frequency [Hz]'
    
#    AmpFList = list(SoundIntensity[FKey].keys()); AmpFList.sort()
    Intensities = [80, 70, 60, 50, 40, 0]
    AmpFList = ControlSoundBoard.dBToAmpF(Intensities, FileName)
    for AFKey in AmpFList:
        AmpFList[AFKey] = [str(_) for _ in AmpFList[AFKey]]
    
    for AKey in AmpFList[FKey]:
        AmpF = AmpFList[FKey].index(AKey)
        
        Colors = [Colormaps[Freq](255-(AmpF*255//len(AmpFList)))]
        
        LineLabel = str(round(SoundIntensity[FKey][AKey])) + ' dB'
        
        for IAmpF in Intensity[FKey]:
            IdB = Intensity[FKey][IAmpF]['dB']
            
            if SoundIntensity[FKey][AKey] == IdB:
#                print(IAmpF)
                Axes[Freq].semilogy(Intensity[FKey][IAmpF]['PSD'][0], 
                                    Intensity[FKey][IAmpF]['PSD'][1], 
                                    label=LineLabel, color=Colors[0])
            else: continue
        
    Axes[Freq].set_xlim(left=5000, right=20000)
    Axes[Freq].set_ylabel(YLabel)
    Axes[Freq].spines['right'].set_visible(False)
    Axes[Freq].spines['top'].set_visible(False)
    Axes[Freq].yaxis.set_ticks_position('left')
    Axes[Freq].xaxis.set_ticks_position('bottom')
    Axes[Freq].set_title(AxTitle)
    if Freq < 2:
        Axes[Freq].legend(loc='lower right')
    else:
        Axes[Freq].legend(loc='lower left')

Fig.suptitle(FigTitle)
Fig.tight_layout()
Fig.subplots_adjust(top=0.93)
Fig.savefig(Folder+'/'+'PowerSpectrumDensity.svg', format='svg')
plt.show()
print('Done.')