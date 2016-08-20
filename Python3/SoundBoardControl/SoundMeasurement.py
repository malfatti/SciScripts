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

Rate = 128000
# Use one that was used in SoundBoardCalibration.py
SoundBoard = 'USBPre2_oAux-iAux'
Setup = 'UnitRec'

## Fill all durations in SECONDS! If

## Sound
# Pulse duration. Avoid using more than 0.5. If you need long pulses, use
# SoundPulseDur = 0.5 and SoundPulseNo = DesiredDurationInSec/0.5
SoundPulseDur = 0.5
# Amount of pulses per block
SoundPulseNo = 4
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], [12000, 14000], 
                  [14000, 16000], [16000, 18000]]
# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0
# Mic sensitivity, from mic datasheet, in dB re V/Pa
MicSens_dB = -47.46

# Path to file saved after Python3/SoundBoardControl/SoundBoardCalibration.py
SBAmpFsFile = '/home/cerebro/Malfatti/Data/Test/20160418173048-SBAmpFs.hdf5'
#==========#==========#==========#==========#

#import array
import ControlSoundBoard
import h5py
import math
import numpy as np
import os
import pandas
import pyaudio
import time
from datetime import datetime
from scipy import signal

with h5py.File(SBAmpFsFile) as h5:
    SBOutAmpF = h5[SoundBoard]['SBOutAmpF'][0]
    SBInAmpF = h5[SoundBoard]['SBInAmpF'][0]

def FRange(Start, End, Step):
    Range = [round(x/(1/Step), 5) 
             for x in range(round(Start/Step), round(End/Step), -1)]
    return(Range)

SoundAmpF = FRange(2.0, 1.0, 0.1) + FRange(1.0, 0.4, 0.05) + \
            FRange(0.4, 0.15, 0.01) + FRange(0.15, 0.03, 0.005) + \
            FRange(0.03, 0.01, 0.0005) + FRange(0.01, 0.001, 0.0001) + \
            FRange(0.001, 0, 0.00002) + [0.0]

#SoundAmpF = [1, 0.5, 0.25, 0]

MicSens_VPa = 10**(MicSens_dB/20)

Date = datetime.now().strftime("%Y%m%d%H%M%S")
Folder = ''.join([Setup, '-', Date, '-SoundMeasurement'])

## Prepare dict w/ experimental setup
FileName = Folder + '/' + Folder + '.hdf5'
DataInfo = dict((Name, eval(Name)) 
                for Name in ['Rate', 'SoundPulseDur', 'SoundPulseNo', 
                             'SoundAmpF', 'NoiseFrequency', 'TTLAmpF', 
                             'MicSens_dB', 'MicSens_VPa', 'SBOutAmpF', 
                             'SBInAmpF', 'Folder'])

## Prepare sound objects
SoundRec = {}
for Freq in range(len(NoiseFrequency)):
    FKey = str(NoiseFrequency[Freq][0]) + '-' + str(NoiseFrequency[Freq][1])
    SoundRec[FKey] = {}
    for AmpF in SoundAmpF:
        SoundRec[FKey][str(AmpF)] = []
        
Sound, Stimulation = ControlSoundBoard.SoundMeasurementOut(
                                   Rate, SoundPulseDur, SoundPulseNo, 
                                   SoundAmpF, NoiseFrequency, TTLAmpF, 
                                   SBOutAmpF
                                                           )

# Define input objects
q = pyaudio.PyAudio()

InOn = False
RecStop = False
def InCallBack(in_data, frame_count, time_info, status):
    if InOn:
        global SoundRec
        SoundRec[Freq][AmpF].append(in_data)
        
    if RecStop:
        InFlag = pyaudio.paComplete
    else:
        InFlag = pyaudio.paContinue
    return(None, InFlag)

Reading = q.open(format=pyaudio.paFloat32,
                     channels=1,
                     rate=Rate,
                     #frames_per_buffer = len(Sound[0][0])//2,
                     input=True,
                     output=False,
                     stream_callback=InCallBack)


## Run!
fTime = (len(SoundAmpF)*len(NoiseFrequency)*(SoundPulseDur*SoundPulseNo))/60
input('Press enter to start sound measurement.')
print('Full test will take', str(round(fTime, 2)), 'min to run.')
print('Current time: ', datetime.now().strftime("%H:%M:%S"))
print('Cover your ears!!')
print('5...')
time.sleep(1)
print('4...')
time.sleep(1)
print('3...')
time.sleep(1)
print('2...')
time.sleep(1)
print('1...')
time.sleep(1)
print('Sound measurement running...')
time.sleep(1)

Reading.start_stream()
for Freq in SoundRec:
    for AmpF in SoundRec[Freq]:
        AmpF = str(AmpF)
        
        RecStop = False; InOn = True        
        for Pulse in range(SoundPulseNo):
            Stimulation.write(Sound[Freq][AmpF])        
        InOn = False; RecStop = True
Reading.stop_stream()

print('Done playing/recording. Saving data...')

## Save!!!
os.makedirs(Folder, exist_ok=True)
with h5py.File(FileName) as F:
    F.create_group('SoundRec')
    for Freq in SoundRec:
        F['SoundRec'].create_group(Freq)
        for AKey, AVal in SoundRec[Freq].items():
            F['SoundRec'][Freq][AKey]= AVal
    
    for Key, Value in DataInfo.items():
        F['SoundRec'].attrs[Key] = Value

del(Sound, Stimulation, q, Reading)
print('Data saved.')


#%% Analysis
# If needed:
#import array
import datetime
#import glob
import h5py
import Hdf5F
import KwikAnalysis
import math
import numpy as np
import pandas
from scipy import signal
DataInfo = Hdf5F.LoadSoundMeasurement(FileName, 'DataInfo')
SoundRec = Hdf5F.LoadSoundMeasurement(FileName, 'SoundRec')

print('Calculating PSD, RMS and dBSLP...')
RecordingData = {}; Intensity = {}

for Freq in SoundRec:
    RecordingData[Freq] = {}
    Intensity[Freq] = {}
    
    for AmpF in SoundRec[Freq]:
        Intensity[Freq][AmpF] = {}
        
        print('Saving data for', Freq, 'at', AmpF)
        
        RecordingData[Freq][AmpF] = np.fromstring(b''.join(SoundRec[Freq][AmpF]), 
                                                  dtype='f')
#        RecordingData[Freq][AmpF] = array.array('f', 
#                                                b''.join(SoundRec[Freq][AmpF]))
#        RecordingData[Freq][AmpF] = array.array('f', SoundRec[Freq][AmpF])
        
        SliceStart = int(DataInfo['Rate']*0.25)-1
        SliceEnd = SliceStart + int(DataInfo['Rate']*1)
        RecordingData[Freq][AmpF] = RecordingData[Freq][AmpF][
                                                        SliceStart:SliceEnd
                                                        ]
        
        RecordingData[Freq][AmpF] = RecordingData[Freq][AmpF]/DataInfo['SBInAmpF']
#        RecordingData[Freq][AmpF] = [_/DataInfo['SBInAmpF']
#                                     for _ in RecordingData[Freq][AmpF]]
        
        Window = signal.hanning(len(RecordingData[Freq][AmpF])//
                                    (DataInfo['Rate']/1000))
        F, PxxSp = signal.welch(RecordingData[Freq][AmpF], DataInfo['Rate'], 
                                Window, nperseg=len(Window), noverlap=0, 
                                scaling='density')
        
        FreqBand = [int(_) for _ in Freq.split('-')]
        
        Start = np.where(F > FreqBand[0])[0][0]-1
        End = np.where(F > FreqBand[1])[0][0]-1
        BinSize = F[1] - F[0]
        RMS = sum(PxxSp[Start:End] * BinSize)**0.5
#        RMS = sum(PxxSp * BinSize)**0.5
        
        Intensity[Freq][AmpF]['PSD'] = [F, PxxSp]
        Intensity[Freq][AmpF]['RMS'] = RMS
        Intensity[Freq][AmpF]['dB'] = 20*(math.log(RMS/DataInfo['MicSens_VPa'], 10)) + 94
        
        del(F, PxxSp, BinSize, RMS)

del(SoundRec)

SoundIntensity = {}
for Freq in Intensity:
    SoundIntensity[Freq] = {}
    
    for AmpF in Intensity[Freq]:
        SoundIntensity[Freq][AmpF] = Intensity[Freq][AmpF]['dB']

## Save analyzed data
print('Saving analyzed data...')
Group = Folder
with h5py.File(FileName) as h5:
    h5.create_group(Group)
    
    for Freq in SoundIntensity:
        h5[Group].create_group(Freq)
        
        for AKey, AVal in SoundIntensity[Freq].items():
            h5[Group][Freq][AKey] = AVal


TexTable = pandas.DataFrame([[DataInfo['SoundAmpF'][AmpF]] + 
                             [Intensity[Freq][AmpF]['dB'] 
                             for Freq in Intensity] 
                             for AmpF in Intensity[Freq]])

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

Params = KwikAnalysis.SetPlot(Backend='TkAgg', Params=True)
from matplotlib import rcParams; rcParams.update(Params)
import matplotlib.pyplot as plt

Colors = ['r', 'g', 'b', 'm', 'k', '#ffa500']
FreqList = list(SoundIntensity.keys()); FreqList.sort(); ToPrepend = []
for FF in ['8000-10000', '9000-11000']:
    if FF in FreqList:
        del(FreqList[FreqList.index(FF)]); ToPrepend.append(FF)

ToPrepend.sort(); FreqList = ToPrepend + FreqList

plt.figure(figsize=(6, 2))
for FKey in SoundIntensity:
    FInd = FreqList.index(FKey)
    AmpFList = sorted(list(SoundIntensity[FKey].keys()), reverse=True)
    ToAppend = AmpFList[:4] + [AmpFList[-1]]
    AmpFList = AmpFList[4:-1] + ToAppend
    
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
plt.show()


Colormaps = [plt.get_cmap(Map) for Map in ['Reds', 'Greens', 'Blues', 
                                           'Purples', 'Greys', 'Oranges']]
Fig, Axes = plt.subplots(len(DataInfo['NoiseFrequency']), 
                         sharex=True, figsize=(6, 12))

for FKey in FreqList:
    Freq = FreqList.index(FKey)
    
    FigTitle = 'Power spectral density'
    AxTitle = FKey + ' Hz'
    YLabel = 'PSD [V$^{2}$/Hz]'
    XLabel = 'Frequency [Hz]'
    
#    AmpFList = list(SoundIntensity[FKey].keys()); AmpFList.sort()
    Intensities = [80, 70, 60, 50, 40, 0]
    AmpFList = ControlSoundBoard.dBToAmpF(Intensities, FileName)
    for AFKey in AmpFList:
        AmpFList[AFKey] = [str(_) for _ in AmpFList[AFKey]]
    
    for AKey in AmpFList[FKey]:
        AmpF = AmpFList[FKey].index(AKey)
        if AKey == '0.0': AKey = '0'
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