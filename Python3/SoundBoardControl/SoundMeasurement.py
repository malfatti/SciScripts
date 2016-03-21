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

## Fill all durations in SECONDS! If

## Sound
# Pulse duration. Avoid using more than 0.5. If you need long pulses, use
# SoundPulseDur = 0.5 and SoundPulseNo = DesiredDurationInSec/0.5
SoundPulseDur = 0.5
# Amount of pulses per block
SoundPulseNo = 4
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
NoiseFrequency = [[8000, 10000], [10000, 12000], [12000, 14000], 
                  [14000, 16000], [16000, 18000]]
# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0
# Mic sensitivity, from mic datasheet, in dB re V/Pa
MicSens_dB = -47.46
# Sound board amp factor. See Python3/SoundBoardControl/SoundBoardCalibration.py
SBOutAmpF = 1.7
SBInAmpF = 0.4852
#==========#==========#==========#==========#

import array
import ControlSoundBoard
import datetime
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas
import pyaudio
import shelve
import time
from scipy import signal

def FRange(Start, End, Step):
    Range = [round(x/(1/Step), 5) 
             for x in range(round(Start/Step), round(End/Step), -1)]
    return(Range)

AmpFList = FRange(2, 1, 0.1) + FRange(1, 0.4, 0.05) + \
           FRange(0.4, 0.15, 0.01) + FRange(0.15, 0.03, 0.005) + \
           FRange(0.03, 0.01, 0.0005) + FRange(0.0001, 0, 0.0001) + \
           FRange(0.0001, 0, 0.00002) + [0]

#AmpFList = FRange(0.025, 0, 0.0005) + [0]

SoundAmpF = [AmpFList for _ in range(len(NoiseFrequency))]
MicSens_VPa = 10**(MicSens_dB/20)

Date = datetime.datetime.now()
Folder = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-SoundMeasurement'])

## Prepare dict w/ experimental setup
DataInfo = dict((Name, eval(Name)) for Name in ['Rate', 'SoundPulseDur', 
                                                'SoundPulseNo', 'SoundAmpF', 
                                                'NoiseFrequency', 'TTLAmpF', 
                                                'MicSens_dB', 'MicSens_VPa',
                                                'SBOutAmpF', 'SBInAmpF',
                                                'Folder'])

## Prepare sound objects
Sound, SoundRec, Stimulation = ControlSoundBoard.SoundMeasurementOut(
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
input('Press enter to start sound measurement.')
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
for Freq in range(len(NoiseFrequency)):
    for AmpF in range(len(SoundAmpF[Freq])):
            RecStop = False
            InOn = True
            for _ in range(SoundPulseNo):
                Stimulation.write(Sound[Freq][AmpF])
            InOn = False
            RecStop = True
Reading.stop_stream()
print('Done playing/recording. Saving data...')

## Save!!!
os.makedirs(Folder)
with shelve.open(Folder + '/' + Folder) as Shelve:
    Shelve['SoundRec'] = SoundRec
    Shelve['DataInfo'] = DataInfo

print('Data saved.')

## Analysis

# If needed:
#
#
#import array
#import math
#import matplotlib.pyplot as plt
#import numpy as np
#import pandas
#import shelve
#from scipy import signal
#Folder = '20160315153450-SoundMeasurement'
#with shelve.open(Folder + '/' + Folder) as Shelve:
#    SoundRec = Shelve['SoundRec']
#    DataInfo = Shelve['DataInfo']
#MicSens_dB = -47.46
#MicSens_VPa = 10**(MicSens_dB/20)
#SBOutAmpF = 1.7
#SBInAmpF = 0.4852

print('Calculating LSD, RMS and dBSLP...')
RecordingData = [0]*len(SoundRec)
Intensity = [0]*len(SoundRec)

for Freq in range(len(SoundRec)):   
    RecordingData[Freq] = [0]*len(SoundRec[Freq])
    Intensity[Freq] = [0]*len(SoundRec[Freq])
    
    for AmpF in range(len(SoundRec[Freq])):
        Intensity[Freq][AmpF] = {}
        
        print('Saving data for ', DataInfo['NoiseFrequency'][Freq], 
              ' at ', DataInfo['SoundAmpF'][Freq][AmpF])
        
        RecordingData[Freq][AmpF] = array.array('f', 
                                                b''.join(SoundRec[Freq][AmpF]))
        
        SliceStart = round(DataInfo['Rate']*0.25)-1
        SliceEnd = SliceStart + round(DataInfo['Rate']*1)
        RecordingData[Freq][AmpF] = RecordingData[Freq][AmpF][
                                                        SliceStart:SliceEnd
                                                        ]
        
        RecordingData[Freq][AmpF] = [_/SBInAmpF
                                     for _ in RecordingData[Freq][AmpF]]
        
        Window = signal.hanning(len(RecordingData[Freq][AmpF])//
                                    (DataInfo['Rate']/1000))
        F, PxxSp = signal.welch(RecordingData[Freq][AmpF], DataInfo['Rate'], 
                                Window, nperseg=len(Window), noverlap=0, 
                                scaling='density')
        
        FreqBand = [DataInfo['NoiseFrequency'][Freq][0], 
                    DataInfo['NoiseFrequency'][Freq][1]]
        
        Start = np.where(F > FreqBand[0])[0][0]-1
        End = np.where(F > FreqBand[1])[0][0]-1
        BinSize = F[1] - F[0]
        RMS = sum(PxxSp[Start:End] * BinSize)**0.5
        
        Intensity[Freq][AmpF]['LSD'] = [F, PxxSp]
        Intensity[Freq][AmpF]['RMS'] = RMS
        Intensity[Freq][AmpF]['dB'] = 20*(math.log(RMS/MicSens_VPa, 10)) + 94
        
        del(F, PxxSp, BinSize, RMS)

#SoundIntensity = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
#for Freq in range(len(DataInfo['NoiseFrequency'])):
#    SoundIntensity[Freq] = {str(DataInfo['SoundAmpF'][Freq][_]): 
#                                Intensity[Freq][_]['dB'] 
#                            for _ in range(len(DataInfo['SoundAmpF'][Freq]))}

SoundIntensity = {}
for Freq in range(len(DataInfo['NoiseFrequency'])):
    Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' + \
          str(DataInfo['NoiseFrequency'][Freq][1])
    SoundIntensity[Key] = {str(DataInfo['SoundAmpF'][Freq][_]): 
                                Intensity[Freq][_]['dB'] 
                            for _ in range(len(DataInfo['SoundAmpF'][Freq]))}
                                
## Save analyzed data
print('Saving analyzed data...')
with shelve.open(Folder + '/SoundIntensity') as Shelve:
    Shelve['SoundIntensity'] = SoundIntensity
    Shelve['SBOutAmpF'] = SBOutAmpF
    Shelve['SBInAmpF'] = SBInAmpF
    Shelve['DataInfo'] = DataInfo

TexTable = pandas.DataFrame([[DataInfo['SoundAmpF'][Freq][AmpF]] + 
                             [Intensity[Freq][AmpF]['dB'] 
                             for Freq in range(
                                             len(DataInfo['NoiseFrequency']))] 
                             for AmpF in range(
                                             len(DataInfo['SoundAmpF'][Freq]))]
                             )

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

Colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y']

for Freq in range(len(DataInfo['NoiseFrequency'])):
    FigTitle = 'Intensity\ curve\ per\ Freq\ per\ AmpF'
    YLabel = 'Intensity\ [dBSPL]'
    XLabel = 'Sound\ amplification\ factor'
    LineLabel = str(DataInfo['NoiseFrequency'][Freq]) + '\ Hz'
    
    plt.plot(DataInfo['SoundAmpF'][Freq], 
             [Intensity[Freq][_]['dB'] 
              for _ in range(len(DataInfo['SoundAmpF'][Freq]))], 
             label='$'+LineLabel +'$', color=Colors[Freq])

plt.ylabel('$'+YLabel+'$'); plt.xlabel('$'+XLabel+'$')
plt.legend(loc='lower right')
plt.locator_params(tight=True)
plt.tick_params(direction='out')
plt.axes().spines['right'].set_visible(False)
plt.axes().spines['top'].set_visible(False)
plt.axes().yaxis.set_ticks_position('left')
plt.axes().xaxis.set_ticks_position('bottom')
plt.title('$'+FigTitle+'$')
plt.savefig(Folder+'/'+'SoundMeasurement.svg', format='svg')


Fig, Axes = plt.subplots(len(DataInfo['NoiseFrequency']), 
                         sharex=True, figsize=(6, 12))
for Freq in range(len(DataInfo['NoiseFrequency'])):
    Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' + \
          str(DataInfo['NoiseFrequency'][Freq][1])
    FigTitle = 'Linear\ spectral\ density'
    AxTitle = str(DataInfo['NoiseFrequency'][Freq]) + '\ Hz'
    YLabel = 'LSD\ [V*RMS/\sqrt{Hz}]'
    XLabel = 'Frequency\ [Hz]'
    
    for AmpF in range(len(DataInfo['SoundAmpF'][Freq])):
        LineLabel = str(SoundIntensity[Key][
                        str(DataInfo['SoundAmpF'][Freq][AmpF])
                        ])[:5] + ' dB'
        
        Axes[Freq].semilogy(Intensity[Freq][AmpF]['LSD'][0], 
                     Intensity[Freq][AmpF]['LSD'][1], 
                     label='$'+LineLabel+'$')
    
    Axes[Freq].set_xlim(left=5000, right=20000)
    Axes[Freq].set_ylabel('$'+YLabel+'$')
    Axes[Freq].set_xlabel('$'+XLabel+'$')
    Axes[Freq].spines['right'].set_visible(False)
    Axes[Freq].spines['top'].set_visible(False)
    Axes[Freq].yaxis.set_ticks_position('left')
    Axes[Freq].xaxis.set_ticks_position('bottom')
    Axes[Freq].set_title('$'+AxTitle+'$')
#    if Freq < 2:
#        Axes[Freq].legend(loc='lower right')
#    else:
#        Axes[Freq].legend(loc='lower left')

Fig.suptitle('$'+FigTitle+'$')
Fig.tight_layout()
Fig.subplots_adjust(top=0.93)
Fig.savefig(Folder+'/'+'LinearSpectrum.svg', format='svg')