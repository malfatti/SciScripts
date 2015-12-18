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
intensities, play them and record them at the same time, for testing purposes. 
Next, it calculates the RMS of each frequency at each intensity. In our setup, 
we use this to calibrate the audio equipment.
"""

#%% Set parameters of the experiment

import array
import datetime
import pyaudio
import pickle
import os
import random
import scipy.signal
import time

import matplotlib.pyplot as plt

#==========#==========#==========#==========#
#==========#==========#==========#==========#
Rate = 128000

## Fill all durations in SECONDS!

## Sound
SoundPulseDur = 0.01
#SoundPostPauseDur = 0.01
# Amount of pulses per block
SoundPulseNo = 200
# Amplification factor (consider using values <= 1)
SoundAmpF = [0.3, 0.2, 0.1]
NoiseFrequency = [[8000, 10000], [10000, 12000], [12000, 14000], [14000, 16000]]

#==========#==========#==========#==========#
#==========#==========#==========#==========#

## Generate sound stimulus
SoundNoise = [random.random() for _ in range(0, round(Rate*SoundPulseDur))]
SoundPulse = [SoundNoise[ElI]*2-1 for ElI,ElV in enumerate(SoundNoise)]
SoundPulseFiltered = [0]*len(NoiseFrequency)
SoundUnit = [0]*len(NoiseFrequency)
SoundList = [0]*len(NoiseFrequency)
Sound = [0]*len(NoiseFrequency)
SoundRec = [0]*len(NoiseFrequency)

for Freq in range(len(NoiseFrequency)):
    passband = [NoiseFrequency[Freq][i]/(Rate/2) for i,j in enumerate(NoiseFrequency[Freq])]
    f2, f1 = scipy.signal.butter(4, passband, 'bandpass')
    SoundPulseFiltered[Freq] = scipy.signal.filtfilt(f2, f1, SoundPulse, padtype='odd', padlen=0)
    SoundPulseFiltered[Freq] = SoundPulseFiltered[Freq].tolist()
    SoundPulseFiltered[Freq][-1] = 0
    
    SoundUnit[Freq] = [0]*len(SoundAmpF)
    SoundList[Freq] = [0]*len(SoundAmpF)
    Sound[Freq] = [0]*len(SoundAmpF)
    SoundRec[Freq] = [[] for _ in range(len(SoundAmpF))]
    
    for AmpF in range(len(SoundAmpF)):
        SoundUnit[Freq][AmpF] = [SoundEl*SoundAmpF[AmpF] for SoundEl in SoundPulseFiltered[Freq]]
        
        SoundList[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))
        Sound[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))
        
        for _ in range(len(SoundUnit[Freq][AmpF])):
            SoundList[Freq][AmpF][_ *2] = SoundUnit[Freq][AmpF][_]
            SoundList[Freq][AmpF][_ *2+1] = 0
        
        Sound[Freq][AmpF] = array.array('f')
        Sound[Freq][AmpF].fromlist(SoundList[Freq][AmpF])
        Sound[Freq][AmpF] = bytes(Sound[Freq][AmpF])


## Generate sound objects and threads
p = pyaudio.PyAudio()
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


Stimulation = p.open(format=pyaudio.paFloat32,
                     channels=2,
                     rate=Rate,
                     #frames_per_buffer = len(Sound[0][0]),
                     input=False,
                     output=True)

Reading = q.open(format=pyaudio.paFloat32,
                     channels=1,
                     rate=Rate,
                     #frames_per_buffer = len(Sound[0][0])//2,
                     input=True,
                     output=False,
                     stream_callback=InCallBack)


#%% Run!
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


#%% Save!!!
Date = datetime.datetime.now()
Folder = Date.strftime("%Y%m%d%H%M%S")
os.makedirs(Folder)
File = open(Folder+'/'+'SoundRec.pckl', 'wb')
pickle.dump(SoundRec, File)
File.close()
del File


#%% Analysis

## If needed:
#import pickle
#import array
#import matplotlib.pyplot as plt
#File = open('SoundRec.pckl', 'rb')
#SoundRec = pickle.load(File)
#File.close()
#del File

RecordingData = [0]*len(SoundRec)
RMS = [0]*len(SoundRec)

for Freq in range(len(SoundRec)):   
    RecordingData[Freq] = [0]*len(SoundRec[Freq])
    RMS[Freq] = [0]*len(SoundRec[Freq])
    
    for AmpF in range(len(SoundRec[Freq])):
        print('Saving data for ', NoiseFrequency[Freq], ' at ', SoundAmpF[AmpF])
        RecordingData[Freq][AmpF] = array.array('f', b''.join(SoundRec[Freq][AmpF]))
        
        DataSquare = [RecordingData[Freq][AmpF][_]**2 for _ in range(len(RecordingData[Freq][AmpF]))]
        RMS[Freq][AmpF] = (sum(DataSquare)/len(DataSquare))**.5
        del DataSquare
    


#%% Save again!
File = open(Folder+'/'+'RecordingData.pckl', 'wb')
pickle.dump(RecordingData, File)
File.close()
del File

File = open(Folder+'/'+'RMS.pckl', 'wb')
pickle.dump(RMS, File)
File.close()
del File
