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

CalibrationFile = '/home/cerebro/Malfatti/Data/Test/' + \
                  '20160419093139-SoundMeasurement/' + \
                  '20160419093139-SoundMeasurement.hdf5'

## Fill all durations in SECONDS! If

## Sound
# Pulse duration. Avoid using more than 0.5. If you need long pulses, use
# SoundPulseDur = 0.5 and SoundPulseNo = DesiredDurationInSec/0.5
SoundPulseDur = 0.5
# Amount of pulses per block
SoundPulseNo = 4
# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0
# Mic sensitivity, from mic datasheet, in dB re V/Pa
MicSens_dB = -47.46

Frequency = '12000-14000'
Intensities = [80, 60, 40, 0]


import array
import ControlSoundBoard
import LoadHdf5Files
import math
import numpy as np
import pyaudio
from scipy import signal

SoundIntensity = LoadHdf5Files.SoundMeasurement(CalibrationFile, 'SoundIntensity')
DataInfo = LoadHdf5Files.SoundMeasurement(CalibrationFile, 'DataInfo')

NoiseFrequency = [[int(Val) for Val in Frequency.split('-')]]
SoundAmpF = [float(min(SoundIntensity[Frequency].keys(), 
              key=lambda i: abs(SoundIntensity[Frequency][i]-dB))) 
              for dB in Intensities]

MicSens_VPa = 10**(MicSens_dB/20)


## Prepare sound objects
Sound, SoundRec, Stimulation = ControlSoundBoard.SoundMeasurementOut(
                                   Rate, SoundPulseDur, SoundPulseNo, 
                                   SoundAmpF, NoiseFrequency, TTLAmpF, 
                                   DataInfo['SBOutAmpF']
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
                     input=True,
                     output=False,
                     stream_callback=InCallBack)

## Start!
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

RecordingData = [0]*len(SoundRec)
Intensity = [0]*len(SoundRec)

for Freq in range(len(SoundRec)):   
    RecordingData[Freq] = [0]*len(SoundRec[Freq])
    Intensity[Freq] = [0]*len(SoundRec[Freq])
    
    for AmpF in range(len(SoundRec[Freq])):
        Intensity[Freq][AmpF] = {}
        
        print('Saving data for ', Frequency, ' at ', SoundAmpF[AmpF])
        
        RecordingData[Freq][AmpF] = array.array('f', 
                                                b''.join(SoundRec[Freq][AmpF]))
        
        SliceStart = int(Rate*0.25)-1
        SliceEnd = SliceStart + int(Rate*1)
        RecordingData[Freq][AmpF] = RecordingData[Freq][AmpF][
                                                        SliceStart:SliceEnd
                                                        ]
        
        RecordingData[Freq][AmpF] = [_/DataInfo['SBInAmpF']
                                     for _ in RecordingData[Freq][AmpF]]
        
        Window = signal.hanning(len(RecordingData[Freq][AmpF])//
                                    (Rate/1000))
        F, PxxSp = signal.welch(RecordingData[Freq][AmpF], Rate, 
                                Window, nperseg=len(Window), noverlap=0, 
                                scaling='density')
        
        FreqBand = [NoiseFrequency[Freq][0], NoiseFrequency[Freq][1]]
        
        Start = np.where(F > FreqBand[0])[0][0]-1
        End = np.where(F > FreqBand[1])[0][0]-1
        BinSize = F[1] - F[0]
        RMS = sum(PxxSp[Start:End] * BinSize)**0.5
        
        Intensity[Freq][AmpF]['LSD'] = [F, PxxSp]
        Intensity[Freq][AmpF]['RMS'] = RMS
        Intensity[Freq][AmpF]['dB'] = 20*(math.log(RMS/MicSens_VPa, 10)) + 94
        
        del(F, PxxSp, BinSize, RMS)
        print('Sent: ', Intensities[AmpF], ' Measured: ', Intensity[Freq][AmpF]['dB'])