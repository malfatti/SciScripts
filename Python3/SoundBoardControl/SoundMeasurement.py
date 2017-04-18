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

## Paths
# Use one that was used in SoundBoardCalibration.py
SBAmpFsFile = '/home/cerebro/Malfatti/Test/20170403123604-SBAmpFs.hdf5'
#SBAmpFsFile = '/home/malfatti/Documents/PhD/Tests/20170214093602-SBAmpFs.hdf5'
SoundSystem = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn'
Setup = 'GPIAS'

Folder = 'SoundMeasurements'
FileName = Folder + '/' + Folder + '.hdf5'
Group = '/'.join([Setup, SoundSystem])

## Sound (Durations in sec)
Rate = 192000
SoundPulseDur = 2
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
#NoiseFrequency = [[8000, 10000], [12000, 14000]]
NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], [12000, 14000], 
                  [14000, 16000], [8000, 16000]]
# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0
# Mic sensitivity, from mic datasheet, in dB re V/Pa or in V/Pa
MicSens_dB = -47.46
#MicSens_VPa = 0.0042364

#==========#==========#==========#==========#

import ControlSoundBoard, DataAnalysis, Hdf5F
import os, pandas
import numpy as np
import sounddevice as SD
from time import sleep
from datetime import datetime

Params = {'backend': 'Qt5Agg'}
from matplotlib import rcParams; rcParams.update(Params)
import matplotlib.pyplot as plt

os.makedirs(Folder, exist_ok=True)

#SoundAmpF = [1, 0.5, 0]
SoundAmpF = np.hstack((
                np.arange(2.15, 1.05, -0.1), np.arange(1.0, 0.4, -0.05),
                np.arange(0.4, 0.15, -0.01), np.arange(0.15, 0.03, -0.005),
                np.arange(0.03, 0.01, -0.0005), np.arange(0.01, 0.001, -0.0001),
                np.arange(0.001, 0, -0.00002), np.array(0.0)
                ))

## Prepare dict w/ experimental setup
Date = datetime.now().strftime("%Y%m%d%H%M%S")
DataInfo = dict((Name, eval(Name)) 
                for Name in ['Rate', 'SoundPulseDur', 'Date','SoundAmpF', 
                             'NoiseFrequency', 'TTLAmpF', 'SoundSystem', 
                             'Setup', 'MicSens_dB', 'Folder'])

Freqs = [str(Freq[0]) + '-' + str(Freq[1]) for Freq in NoiseFrequency]
SoundAmpF = {Freq: SoundAmpF for Freq in Freqs}

# Prepare audio objects
SD.default.device = 'system'
SD.default.samplerate = Rate
SD.default.channels = 2

# Warn user
FullTime = (len(SoundAmpF[Freqs[0]])*len(NoiseFrequency)*(SoundPulseDur))/60
FullTime = str(round(FullTime, 2))
print('Full test will take', FullTime, 'min to run.')
print('Current time: ', datetime.now().strftime("%H:%M:%S"))
print('This can be loud - cover your ears!!')

for i in range(5, 0, -1): print(i, end=' '); sleep(1)

print('Sound measurement running...')
SoundRec = {}
for Freq in NoiseFrequency:
    SoundRec[str(Freq[0]) + '-' + str(Freq[1])] = {}
    Sound = ControlSoundBoard.SoundStim(Rate, SoundPulseDur, SoundAmpF, 
                                        [Freq], TTLAmpF, SoundSystem, 
                                        TTLs=False, Map=[2,1])
    
    for FKey in Sound:
        for AKey in Sound[FKey]:
            print(FKey, AKey)
            SoundRec[FKey][AKey] = SD.playrec(Sound[FKey][AKey], blocking=True)
    
    print('Done playing/recording', FKey + '.')
    Hdf5F.WriteSoundMeasurement(SoundRec, DataInfo, Group, FileName)
    del(Sound, SoundRec)

print('Finished recording \O/')


#%% Analysis
#DataInfo = Hdf5F.LoadSoundMeasurement(FileName, 'DataInfo')
SoundRec = Hdf5F.LoadSoundMeasurement(FileName, Group, 'SoundRec')
SBInAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundSystem, 'SBInAmpF')

if 'MicSens_dB' in DataInfo:
    DataInfo['MicSens_VPa'] = 10**(DataInfo['MicSens_dB']/20)
    dB = True
#elif 'MicSens_VPa' not in DataInfo:
#    print('No mic sensitivity info. Impossible to calculate dBSPL values.')
#    dB = False

print('Calculating PSD, RMS and dBSLP...')
SoundIntensity = {}

for FKey in SoundRec:
    SoundIntensity[FKey] = {}
    
    for AKey in SoundRec[FKey]:
        SoundRec[FKey][AKey] = SoundRec[FKey][AKey] * SBInAmpF
        FreqBand = [int(_) for _ in FKey.split('-')]
        
        Intensity = DataAnalysis.SignalIntensity(SoundRec[FKey][AKey], 
                                                 DataInfo['Rate'], FreqBand, 
                                                 DataInfo['MicSens_VPa'])
        
        SoundIntensity[FKey][AKey] = Intensity['dB']

del(SoundRec)


## Save analyzed data
print('Saving analyzed data...')
os.makedirs(Folder, exist_ok=True)
Hdf5F.WriteSoundIntensity(SoundIntensity, Group, FileName)


AmpFList = sorted(list(SoundIntensity['8000-10000'].keys()), reverse=True)
#ToAppend = AmpFList[:4] + [AmpFList[-1]]
#AmpFList = AmpFList[4:-1] + ToAppend

TexTable = pandas.DataFrame([[float(AmpF)] + [Intensity[Freq][AmpF]['dB'] 
                                              for Freq in Intensity]
                             for AmpF in AmpFList])

File = open(Folder+'/'+'IntensityTable_'+DataInfo['Date']+'.tex', 'w')
File.write(r"""%% Configs =====
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
\input{IntensityTable-Contents""" +'_'+DataInfo['Date']+'.tex}')
File.write(r"""
\end{document}
""")
File.close()
del(File)

File = open(Folder+'/'+'IntensityTable-Contents_'+DataInfo['Date']+'.tex', 'w')
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
plt.savefig(Folder+'/'+'SoundMeasurement_'+DataInfo['Date']+'.svg', format='svg')


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
    AmpFList = ControlSoundBoard.dBToAmpF(Intensities, FileName, Group)
    for AFKey in AmpFList:
        AmpFList[AFKey] = [str(_) for _ in AmpFList[AFKey]]
    
    for AKey in AmpFList[FKey]:
        AmpF = AmpFList[FKey].index(AKey)
        
        Colors = [Colormaps[Freq](255-(AmpF*255//len(AmpFList)))]
        
        try:
            LineLabel = str(round(SoundIntensity[FKey][AKey])) + ' dB'
        except KeyError:
            AKey = '0'
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
Fig.savefig(Folder+'/'+'PowerSpectrumDensity_'+DataInfo['Date']+'.svg', format='svg')
plt.show()
print('Done.')