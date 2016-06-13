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

## Fill all durations in SECONDS! If

## Sound
# Pulse duration. Avoid using more than 0.5. If you need long pulses, use
# SoundPulseDur = 0.5 and SoundPulseNo = DesiredDurationInSec/0.5
SoundPulseDur = 0.5
# Amount of pulses per block
SoundPulseNo = 4
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
#NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], [12000, 14000], 
#                  [14000, 16000], [16000, 18000]]
NoiseFrequency = [[9000, 11000]]
# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0
# Mic sensitivity, from mic datasheet, in dB re V/Pa
MicSens_dB = -47.46

# Path to file saved after Python3/SoundBoardControl/SoundBoardCalibration.py
SBAmpFsFile = '/home/cerebro/Malfatti/Data/Test/20160418173048-SBAmpFs.hdf5'
#==========#==========#==========#==========#

import array
import ControlSoundBoard
import datetime
import h5py
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas
import pyaudio
import time
from scipy import signal

with h5py.File(SBAmpFsFile) as h5:
    SBOutAmpF = h5[SoundBoard]['SBOutAmpF'][0]
    SBInAmpF = h5[SoundBoard]['SBInAmpF'][0]

def FRange(Start, End, Step):
    Range = [round(x/(1/Step), 5) 
             for x in range(round(Start/Step), round(End/Step), -1)]
    return(Range)

SoundAmpF = FRange(2, 1, 0.1) + FRange(1, 0.4, 0.05) + \
            FRange(0.4, 0.15, 0.01) + FRange(0.15, 0.03, 0.005) + \
            FRange(0.03, 0.01, 0.0005) + FRange(0.01, 0.001, 0.0001) + \
            FRange(0.001, 0, 0.00002) + [0]

#SoundAmpF = [1, 0.5, 0.25, 0]

MicSens_VPa = 10**(MicSens_dB/20)

#Date = datetime.datetime.now()
#Folder = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-SoundMeasurement'])

## Prepare dict w/ experimental setup
#FileName = Folder + '/' + Folder + '.hdf5'
FileName = '20160419093139-SoundMeasurement/20160419093139-SoundMeasurement.hdf5'
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
fTime = (len(SoundAmpF)*len(NoiseFrequency)*(SoundPulseDur*SoundPulseNo))/60
input('Press enter to start sound measurement.')
print('Full test will take', str(round(fTime, 2)), 'min to run.')
print('Current time: ', datetime.datetime.now().strftime("%H:%M:%S"))
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
    for AmpF in range(len(SoundAmpF)):
            RecStop = False
            InOn = True
            for _ in range(SoundPulseNo):
                Stimulation.write(Sound[Freq][AmpF])
            InOn = False
            RecStop = True
Reading.stop_stream()
print('Done playing/recording. Saving data...')

## Save!!!
#os.makedirs(Folder, exist_ok=True)
with h5py.File(FileName) as h5:
#    h5.create_group('SoundRec')
    for Freq in range(len(SoundRec)):
        Key = str(NoiseFrequency[Freq][0]) + '-' + str(NoiseFrequency[Freq][1])
        h5['SoundRec'].create_group(Key)
        for AmpF in range(len(SoundRec[Freq])):
            h5['SoundRec'][Key].create_dataset(str(SoundAmpF[AmpF]), 
                                                data=SoundRec[Freq][AmpF])
    
    for Key, Value in DataInfo.items():
        h5['SoundRec'].attrs[Key] = Value

del(Sound, Stimulation, q, Reading)
print('Data saved.')


## Analysis
# If needed:
#import array
#import datetime
#import glob
#import h5py
#import LoadHdf5Files
#import math
#import matplotlib.pyplot as plt
#import numpy as np
#import pandas
#from scipy import signal
#FileName = glob.glob('*.hdf5'); FileName = FileName[0]
#DataInfo = LoadHdf5Files.SoundMeasurement(FileName, 'DataInfo')
#SoundRec = LoadHdf5Files.SoundMeasurement(FileName, 'SoundRec')

print('Calculating LSD, RMS and dBSLP...')
RecordingData = [0]*len(DataInfo['NoiseFrequency'])
Intensity = [0]*len(DataInfo['NoiseFrequency'])

for Freq in range(len(DataInfo['NoiseFrequency'])):   
    RecordingData[Freq] = [0]*len(DataInfo['SoundAmpF'])
    Intensity[Freq] = [0]*len(DataInfo['SoundAmpF'])
    
    for AmpF in range(len(DataInfo['SoundAmpF'])):
        Intensity[Freq][AmpF] = {}
        
        print('Saving data for ', DataInfo['NoiseFrequency'][Freq], 
              ' at ', DataInfo['SoundAmpF'][AmpF])
        
        RecordingData[Freq][AmpF] = array.array('f', 
                                                b''.join(SoundRec[Freq][AmpF]))
        
        SliceStart = int(DataInfo['Rate']*0.25)-1
        SliceEnd = SliceStart + int(DataInfo['Rate']*1)
        RecordingData[Freq][AmpF] = RecordingData[Freq][AmpF][
                                                        SliceStart:SliceEnd
                                                        ]
        
        RecordingData[Freq][AmpF] = [_/DataInfo['SBInAmpF']
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
#        RMS = sum(PxxSp * BinSize)**0.5
        
        Intensity[Freq][AmpF]['LSD'] = [F, PxxSp]
        Intensity[Freq][AmpF]['RMS'] = RMS
        Intensity[Freq][AmpF]['dB'] = 20*(math.log(RMS/DataInfo['MicSens_VPa'], 10)) + 94
        
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
    SoundIntensity[Key] = {str(DataInfo['SoundAmpF'][_]): 
                                Intensity[Freq][_]['dB'] 
                            for _ in range(len(DataInfo['SoundAmpF']))}
                                
## Save analyzed data
print('Saving analyzed data...')
#GroupName = 'SoundIntensity-' + datetime.datetime.now().strftime("%Y%m%d%H%M%S")
with h5py.File(FileName) as h5:
#    if 'SoundIntensity' in h5.keys():
#            h5[GroupName] = h5['SoundIntensity']; del(h5['SoundIntensity'])
#    
#    h5.create_group('SoundIntensity')
    for Freq in range(len(DataInfo['NoiseFrequency'])):
        Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' + \
              str(DataInfo['NoiseFrequency'][Freq][1])
        h5['SoundIntensity'].create_group(Key)
        for AmpF in range(len(DataInfo['SoundAmpF'])):
            h5['SoundIntensity'][Key].create_dataset(
                str(DataInfo['SoundAmpF'][AmpF]), 
                data=[Intensity[Freq][AmpF]['dB']]
                                                     )

TexTable = pandas.DataFrame([[DataInfo['SoundAmpF'][AmpF]] + 
                             [Intensity[Freq][AmpF]['dB'] 
                             for Freq in range(
                                             len(DataInfo['NoiseFrequency']))] 
                             for AmpF in range(
                                             len(DataInfo['SoundAmpF']))]
                             )

#File = open(Folder+'/'+'IntensityTable.tex', 'w')
File = open('IntensityTable.tex', 'w')
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

#File = open(Folder+'/'+'IntensityTable-Contents.tex', 'w')
File = open('IntensityTable-Contents.tex', 'w')
File.write(TexTable.to_latex(longtable=True))
File.close()
del(File)

Colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y']

for Freq in range(len(DataInfo['NoiseFrequency'])):
    Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' + \
          str(DataInfo['NoiseFrequency'][Freq][1])
    FigTitle = 'Intensity\ curve\ per\ Freq\ per\ AmpF'
    YLabel = 'Intensity\ [\si{\dB SPL}]'
    XLabel = 'Sound\ amplification\ factor'
    LineLabel = '$'+ Key + '\ \si{\hertz}' +'$'
    
    plt.plot(DataInfo['SoundAmpF'], 
             [Intensity[Freq][_]['dB'] 
              for _ in range(len(DataInfo['SoundAmpF']))], 
             label=LineLabel, color=Colors[Freq])

plt.ylabel('$'+YLabel+'$'); plt.xlabel('$'+XLabel+'$')
plt.legend(loc='lower right')
plt.locator_params(tight=True)
plt.tick_params(direction='out')
plt.axes().spines['right'].set_visible(False)
plt.axes().spines['top'].set_visible(False)
plt.axes().yaxis.set_ticks_position('left')
plt.axes().xaxis.set_ticks_position('bottom')
plt.title(FigTitle)
#plt.savefig(Folder+'/'+'SoundMeasurement.svg', format='svg')
plt.savefig('SoundMeasurement.svg', format='svg')


Fig, Axes = plt.subplots(len(DataInfo['NoiseFrequency']), 
                         sharex=True, figsize=(6, 12))
for Freq in range(len(DataInfo['NoiseFrequency'])):
    Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' + \
          str(DataInfo['NoiseFrequency'][Freq][1])
    FigTitle = 'Linear\ spectral\ density'
    AxTitle = Key + '\ \si{\hertz}'
    YLabel = 'LSD\ [\si{V *RMS}/\sqrt{\si{\hertz}}]'
    XLabel = 'Frequency\ [\si{\hertz}]'
    
    for AmpF in range(len(DataInfo['SoundAmpF'])):
        LineLabel = '$'+str(SoundIntensity[Key][
                        str(DataInfo['SoundAmpF'][AmpF])
                        ])[:5] + ' \si{\dB}' +'$'
        
        Axes[Freq].semilogy(Intensity[Freq][AmpF]['LSD'][0], 
                     Intensity[Freq][AmpF]['LSD'][1], 
                     label=LineLabel)
    
    Axes[Freq].set_xlim(left=5000, right=20000)
    Axes[Freq].set_ylabel('$'+YLabel+'$')
    Axes[Freq].set_xlabel('$'+XLabel+'$')
    Axes[Freq].spines['right'].set_visible(False)
    Axes[Freq].spines['top'].set_visible(False)
    Axes[Freq].yaxis.set_ticks_position('left')
    Axes[Freq].xaxis.set_ticks_position('bottom')
    Axes[Freq].set_title(AxTitle)
#    if Freq < 2:
#        Axes[Freq].legend(loc='lower right')
#    else:
#        Axes[Freq].legend(loc='lower left')

Fig.suptitle(FigTitle)
Fig.tight_layout()
Fig.subplots_adjust(top=0.93)
#Fig.savefig(Folder+'/'+'LinearSpectrum.svg', format='svg')
Fig.savefig('LinearSpectrum.svg', format='svg')

print('Done.')