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

This is a script to record sound from sound board input in sinchrony with a 
DAQ. Next, it calculates the intensity in RMS and dB for each frequency at each
amplification factor. In our setup, we use this to calibrate the audio 
equipment.
"""

#%% Set parameters of the experiment

Rate = 128000

# Amplification factors. Set the same in the DAQ.
SoundAmpF =  [2, 1.75, 1.5, 1.25, 1, 0.9, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 
              0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.08, 0.06, 0.04, 
              0.03, 0.02, 0.015, 0.01, 0.008, 0.006, 0.005, 0.004, 0.003, 
              0.002, 0.001, 0]
# Noise frequency. Set the same in the DAQ.
NoiseFrequency = [[8000, 10000], [10000, 12000], [12000, 14000], 
                  [14000, 16000], [16000, 18000]]
# Time to record in seconds.
TimeRec = 2
# Mic sensitivity, from mic datasheet, in dB re V/Pa.
MicSens_dB = -47.46
#==========#==========#==========#==========#

import array
#import ControlSoundBoard
import datetime
import math
import matplotlib.pyplot as plt
import os
import pandas
import pickle
import pyaudio
import serial
import serial.tools.list_ports
from scipy import signal

Date = datetime.datetime.now()
Folder = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-SoundMeasurementDAQ'])
os.makedirs(Folder)

DataInfo = [Rate, SoundAmpF, NoiseFrequency, TimeRec, MicSens_dB, Folder]

# Define in/out objects
q = pyaudio.PyAudio()
Reading = q.open(format=pyaudio.paFloat32,
                     channels=1,
                     rate=Rate,
                     input=True,
                     output=False)

Port = serial.tools.list_ports.comports()
Arduino = serial.Serial(Port[-1][0], 19200)

# Preallocate memory
Time =  Rate * (TimeRec+1)

SoundRec = [[]*_ for _ in range(len(NoiseFrequency))]
for Freq in range(len(NoiseFrequency)):
    SoundRec[Freq] = [[]*_ for _ in range(len(SoundAmpF))]

#%% Run!
print('Sound measurement running...')
for Freq in range(len(NoiseFrequency)):
    for AmpF in range(len(SoundAmpF)):
        Arduino.write(b'A')
        SoundRec[Freq][AmpF] = Reading.read(Time)
print('Done playing/recording. Saving data...')

## Save!!!

File = open(Folder+'/'+'SoundRec.pckl', 'wb')
pickle.dump(SoundRec, File)
File.close()
del(File)

File = open(Folder+'/'+'DataInfo.pckl', 'wb')
pickle.dump(DataInfo, File)
File.close()
del(File)
print('Data saved.')

#%% Analysis

## If needed:
#File = open(Folder+'/'+'SoundRec.pckl', 'rb')
#SoundRec = pickle.load(File)
#File.close()
#del(File)
#
#File = open(Folder+'/'+'DataInfo.pckl', 'rb')
#DataInfo = pickle.load(File)
#File.close()
#del(File)
#
#Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, \
#NoiseFrequency, TTLAmpF, MicSens_dB, Folder = DataInfo

print('Calculating PSD, RMS and dBSLP...')
RecordingData = [0]*len(SoundRec)
Intensity = [0]*len(SoundRec)
MicSens_VPa = 10**(MicSens_dB/20)

for Freq in range(len(SoundRec)):   
    RecordingData[Freq] = [0]*len(SoundRec[Freq])
    Intensity[Freq] = [0]*len(SoundRec[Freq])
    
    for AmpF in range(len(SoundRec[Freq])):
        Intensity[Freq][AmpF] = {}
        print('Saving data for ', NoiseFrequency[Freq], 
              ' at ', SoundAmpF[AmpF])
        
        RecordingData[Freq][AmpF] = array.array('f', SoundRec[Freq][AmpF])
        RecordingData[Freq][AmpF] = RecordingData[Freq][AmpF][
                                    round(Rate*0.1)-1:round(Rate*0.1)+Rate]
        
        F, PxxSp = signal.welch(RecordingData[Freq][AmpF], Rate, 
                                nperseg=1024, scaling='spectrum')
        
        Intensity[Freq][AmpF]['PSD'] = [F, PxxSp]
        
        Intensity[Freq][AmpF]['RMS'] = (sum(Intensity[Freq][AmpF]['PSD'][1])*
                                        (Intensity[Freq][AmpF]['PSD'][0][1] -
                                         Intensity[Freq][AmpF]['PSD'][0][0])
                                         )**0.5
        
        Intensity[Freq][AmpF]['dB'] = 20*(math.log(
                                    Intensity[Freq][AmpF]['RMS']/0.00002, 10))
        
        del(F, PxxSp)


## Save analyzed data
print('Saving analyzed data...')
TexTable = pandas.DataFrame([[SoundAmpF[AmpF]] + 
                             [Intensity[Freq][AmpF]['dB'] 
                             for Freq in range(len(NoiseFrequency))] 
                             for AmpF in range(len(SoundAmpF))])

File = open(Folder+'/'+'IntensityTable.tex', 'w')
File.write(r"""
%% Configs =====
\documentclass[12pt,a4paper]{report}

\usepackage{lmodern}
\usepackage[utf8]{inputenc}
\usepackage[english]{babel}
\usepackage{graphicx}
\usepackage[font=small,labelfont=bf,justification=justified,singlelinecheck=false]{caption}

\usepackage{indentfirst}

\usepackage{siunitx}

\usepackage[left=0.5cm,right=0.5cm,top=0.5cm,bottom=0.5cm]{geometry}
\usepackage{setspace}
\usepackage{titlesec}
\titleformat{\chapter}{\normalfont\LARGE\bfseries}{\thechapter.}{1em}{}

\renewcommand{\rmdefault}{phv}
\renewcommand{\sfdefault}{phv}

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
File.write(TexTable.to_latex())
File.close()
del(File)

Colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y']

plt.figure(1)
for Freq in range(len(NoiseFrequency)):
    plt.plot(SoundAmpF, 
             [Intensity[Freq][_]['dB'] for _ in range(len(SoundAmpF))], 
             label=str(NoiseFrequency[Freq]), color=Colors[Freq])
plt.ylabel('Intensity [dBSPL]'); plt.xlabel('Sound amplification factor')
plt.legend(loc='best', frameon=False)
plt.locator_params(tight=True)
plt.tick_params(direction='out')
plt.axes().spines['right'].set_visible(False)
plt.axes().spines['top'].set_visible(False)
plt.axes().yaxis.set_ticks_position('left')
plt.axes().xaxis.set_ticks_position('bottom')
plt.show()
#input('Press enter to save figure 1.')
plt.savefig(Folder+'/'+'SoundMeasurement.pdf', transparent=True)

#%%
F, ((A, B), (C, D)) = plt.subplots(2, 2)
Axes = [A, B, C, D]
for Freq in range(len(NoiseFrequency)):
    for AmpF in range(len(SoundAmpF)):
        Axes[Freq].semilogy(Intensity[Freq][AmpF]['PSD'][0], 
                     Intensity[Freq][AmpF]['PSD'][1], 
                     label=str(SoundAmpF[AmpF]), color=Colors[Freq])
    Axes[Freq].set_xlim(left=5000, right=20000)
    Axes[Freq].set_ylabel('Linear spectrum [V RMS]')
    Axes[Freq].set_xlabel('Frequency [Hz]')
    Axes[Freq].spines['right'].set_visible(False)
    Axes[Freq].spines['top'].set_visible(False)
    Axes[Freq].yaxis.set_ticks_position('left')
    Axes[Freq].xaxis.set_ticks_position('bottom')
    Axes[Freq].set_title(str(NoiseFrequency[Freq]))
F.show()
#input('Press enter to save figure 2.')
F.savefig(Folder+'/'+'LinearSpectrum.pdf', transparent=True)
