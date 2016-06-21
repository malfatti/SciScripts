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

"""
#%% Set parameters

AnimalName = 'TTLsLatencyTest'
Rate = 128000
BaudRate = 38400

CalibrationFile = '/home/cerebro/Malfatti/Data/Test/' + \
                  '20160419093139-SoundMeasurement/' + \
                  '20160419093139-SoundMeasurement.hdf5'

# Sound board used
SoundBoard = 'USBPre2_oAux-iAux'

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 1

## Fill all durations in SECONDS!

## Sound
# Silence before pulse
SoundPrePauseDur = 0.004
# Pulse duration
SoundPulseDur = 0.003
# Silence after pulse
SoundPostPauseDur = 0.013
# Amount of pulses per block
SoundPulseNo = 4000
# Number of blocks
SoundStimBlockNo = 1
# Duration of pause between blocks
SoundPauseBetweenStimBlocksDur = 1
# Intensities tested, in order, in dB. Supports floats :)
#Intensities = [80, 75, 70, 65, 60, 55, 50, 45]
Intensities = [80]
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
# USE ONLY FREQUENCY BANDS THAT WERE CALIBRATED. To check the calibrated freqs, 
# just run the cell once and then list(SoundIntensity).
NoiseFrequency = [[8000, 10000]]


import ControlArduino
import ControlSoundBoard
import datetime
import h5py
import LoadHdf5Files
import numpy as np
import pyaudio
import os

SoundIntensity = LoadHdf5Files.SoundMeasurement(CalibrationFile, 
                                                'SoundIntensity')

Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-SoundStim.hdf5'])

Folder = FileName[:-5]
os.makedirs(Folder, exist_ok=True)    # Figs folder

SoundAmpF = {Hz: [float(min(SoundIntensity[Hz].keys(), 
              key=lambda i: abs(SoundIntensity[Hz][i]-dB))) 
              for dB in Intensities] 
         for Hz in list(SoundIntensity)}

DataInfo = dict((Name, eval(Name)) for Name in ['AnimalName', 'Rate', 
                                       'BaudRate', 'SoundPrePauseDur', 
                                       'SoundPulseDur', 'SoundPostPauseDur', 
                                       'SoundPulseNo', 'SoundStimBlockNo', 
                                       'SoundPauseBetweenStimBlocksDur', 
                                       'TTLAmpF', 'Intensities', 
                                       'NoiseFrequency', 'CalibrationFile', 
                                       'SoundBoard', 'FileName'])

with h5py.File(Folder+'/'+FileName) as h5:
    h5.create_group('DataInfo')
    for Key, Value in DataInfo.items():
        h5['DataInfo'].attrs[Key] = Value
    
    h5['DataInfo'].create_group('SoundAmpF')
    for Key, Value in SoundAmpF.items():
        h5['DataInfo']['SoundAmpF'][Key] = Value


## Create Audio and Arduino objects
p = pyaudio.PyAudio()
Stimulation = p.open(format=pyaudio.paFloat32,
                channels=2,
                rate=Rate,
                output=True)

Arduino = ControlArduino.CreateObj(BaudRate)

Sound, PlaySound = ControlSoundBoard.GenSound(Rate, SoundPulseDur, 
                                              SoundPulseNo, SoundAmpF, 
                                              NoiseFrequency, TTLAmpF, 
                                              CalibrationFile, SoundBoard, 
                                              'AllPulses', SoundPrePauseDur, 
                                              SoundPostPauseDur, 
                                              SoundStimBlockNo, 
                                              SoundPauseBetweenStimBlocksDur)

#%% Run!

DVCoord = 'Out'
Freq = 0

print('Running...')
Key = str(NoiseFrequency[Freq][0]) + '-' + str(NoiseFrequency[Freq][1])
for AmpF in range(len(SoundAmpF[Key])):
    print('Playing', str(NoiseFrequency[Freq]), 'at', 
          str(Intensities[AmpF]), 'dB')
    
    Arduino.write(b'P')
    PlaySound(Freq, AmpF)
    Arduino.write(b'P')

print('Done. Saving info...')
with h5py.File(Folder+'/'+FileName) as h5:
    h5.create_group('ExpInfo')
    h5['ExpInfo'].create_group(str(len(list(h5['ExpInfo']))))
    Key = list(h5['ExpInfo'].keys())[-1]
    
    h5['ExpInfo'][Key].attrs['StimType'] = [np.string_('Sound')]
    h5['ExpInfo'][Key].attrs['DVCoord'] = DVCoord
    h5['ExpInfo'][Key].attrs['Hz'] = Freq

print('Saved.')
print('Played Freq', str(Freq), 'at', DVCoord, 'µm DV')
