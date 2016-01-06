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

"""==========#==========#==========#=========="""
Rate = 128000

## Fill all durations in SECONDS! If

## Sound
# Silence before pulse
SoundPrePauseDur = 0
# Pulse duration
SoundPulseDur = 0.01
# Silence after pulse
SoundPostPauseDur = 0
# Amount of pulses per block
SoundPulseNo = 300
# Number of blocks
SoundStimBlockNo = 1
# Duration of pause between blocks
SoundPauseBetweenStimBlocksDur = 0
# Amplification factor (consider using values <= 1). If using one value, keep it in a list.
SoundAmpF = [0.3, 0.2, 0.1]
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
NoiseFrequency = [[8000, 10000], [10000, 12000], [12000, 14000], [14000, 16000]]
# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 1
# Mic sensitivity, from mic datasheet, in dB re V/Pa
MicSens_dB = -47.46
"""==========#==========#==========#=========="""

import array
import ControlSoundBoard
import datetime
import math
import os
import pickle
import pyaudio


## Prepare dict w/ experimental setup
DataInfo = [Rate, SoundPrePauseDur, SoundPulseDur, SoundPostPauseDur, 
            SoundPulseNo, SoundStimBlockNo, SoundPauseBetweenStimBlocksDur, 
            SoundAmpF, NoiseFrequency, TTLAmpF]


## Prepare sound objects
Sound, SoundPauseBetweenStimBlocks, SoundRec, Stimulation = \
    ControlSoundBoard.GenSoundRec(Rate, SoundPrePauseDur, SoundPulseDur, 
                                  SoundPostPauseDur, SoundPulseNo, 
                                  SoundStimBlockNo, 
                                  SoundPauseBetweenStimBlocksDur, SoundAmpF, 
                                  NoiseFrequency, TTLAmpF)

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
print('Sound measurement running...')
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
print('Done. Saving data...')
#PlayRec(Sound, SoundPulseNo, SoundAmpF, NoiseFrequency, Stimulation, Reading)


## Save!!!
Date = datetime.datetime.now()
Folder = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-SoundMeasurement'])
os.makedirs(Folder)

File = open(Folder+'/'+'SoundRec.pckl', 'wb')
pickle.dump(SoundRec, File)
File.close()
del(File)

File = open(Folder+'/'+'DataInfo.pckl', 'wb')
pickle.dump(DataInfo, File)
File.close()
del(File)
print('Done saving data.')

#Continue = input(''.join(['\n',
#'You can stop the program now and run analysis later. Do you want to run the ','\n',
#'analysis now?','\n',
#'[y/N]: ']))
#
#if Continue in ('y', 'Y', 'yes', 'Yes', 'YES'):
#    pass
#else:
#    How to stop the script?
#
#print('Test.')
#%% Analysis

## If needed:
#import datetime
#import os
#import pickle
#import array
#import math
#import matplotlib.pyplot as plt
#
#File = open('SoundRec.pckl', 'rb')
#SoundRec = pickle.load(File)
#File.close()
#del(File)
#
#File = open('DataInfo.pckl', 'rb')
#DataInfo = pickle.load(File)
#File.close()
#del(File)
#
#Rate, SoundPrePauseDur, SoundPulseDur, SoundPostPauseDur, SoundPulseNo, \
#    SoundStimBlockNo, SoundPauseBetweenStimBlocksDur, SoundAmpF, NoiseFrequency, \
#    TTLAmpF = DataInfo

RecordingData = [[0]*len(SoundRec[Freq])]*len(SoundRec)
Intensity = [0]*len(SoundRec)
MicSens_VPa = 10**(MicSens_dB/20)

for Freq in range(len(SoundRec)):   
    RecordingData[Freq] = [0]*len(SoundRec[Freq])
    Intensity[Freq] = [0]*len(SoundRec[Freq])
    
    for AmpF in range(len(SoundRec[Freq])):
        Intensity[Freq][AmpF] = {}
        print('Saving data for ', NoiseFrequency[Freq], ' at ', SoundAmpF[AmpF])
        
        RecordingData[Freq][AmpF] = array.array('f', b''.join(SoundRec[Freq][AmpF]))
        DataSquare = [RecordingData[Freq][AmpF][_]**2 for _ in range(len(RecordingData[Freq][AmpF]))]
        
        Intensity[Freq][AmpF]['RMS'] = (sum(DataSquare)/len(DataSquare))**.5
        P = Intensity[Freq][AmpF]['RMS']/MicSens_VPa
        Intensity[Freq][AmpF]['dB'] = 20*(math.log(P/0.00002, 10))
        
        del(DataSquare, P)


#%% Save again!
File = open(Folder+'/'+'RecordingData.pckl', 'wb')
pickle.dump(RecordingData, File)
File.close()
del(File)

File = open(Folder+'/'+'Intensity.pckl', 'wb')
pickle.dump(Intensity, File)
File.close()
del(File)