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

This is a script to generate pulses and send to the soundboard, and then to a 
sound amplifier and an Arduino board. Basically it generates sound pulses, 
sound square pulses (TTLs), and laser square pulses. The square pulses will be 
sent to the left channel and the sound pulses will be sent to the right 
channel. In our current setup, we have an arduino uno set to read an analog 
pin and when the voltage is within a range, it turns on a digital port (Check 
ControlArduinoWithSoundBoard.ino code). This way, you can control one sound 
channel, one laser and "separate" TTLs for both.

The first cell will set the experiment parameters.

The second, third and fourth cells create audio objects for sound playback and 
set some attributes of it. Then, threads are created so .write() commands can 
run on its own thread.

All the following cells send the stimulus to the sound board, each one with its 
own settings. 
"""
#%% Set Parameters
AnimalName = 'CaMKIIahM4Dn07'
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
SoundPostPauseDur = 0.093
# Amount of pulses per block
SoundPulseNo = 529
# Number of blocks
SoundStimBlockNo = 1
# Duration of pause between blocks
SoundPauseBetweenStimBlocksDur = 10
# Intensities tested, in order, in dB. Supports floats :)
Intensities = [80, 75, 70, 65, 60, 55, 50, 45, 40, 35]
#Intensities = [100, 75, 70]
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
# USE ONLY FREQUENCY BANDS THAT WERE CALIBRATED. To check the calibrated freqs, 
# just run the cell once and then list(SoundIntensity).
NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], 
                  [12000, 14000], [14000, 16000]]
#NoiseFrequency = [[8000, 10000]]

## Laser
# Silence before pulse
LaserPrePauseDur = 0
# Pulse duration
LaserPulseDur = 0.01
# Silence after pulse
LaserPostPauseDur = 0.09
# Amount of pulses per block
LaserPulseNo = 200
# Number of blocks
LaserStimBlockNo = 1
# Duration of pause between blocks
LaserPauseBetweenStimBlocksDur = 5
#==========#==========#==========#==========#

import ControlArduino
import ControlSoundBoard
import datetime
import h5py
import LoadHdf5Files
import numpy as np
import pyaudio

SoundIntensity = LoadHdf5Files.SoundMeasurement(CalibrationFile, 
                                                'SoundIntensity')

Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-SoundStim.hdf5'])

SoundAmpF = {Hz: [float(min(SoundIntensity[Hz].keys(), 
              key=lambda i: abs(SoundIntensity[Hz][i]-dB))) 
              for dB in Intensities] 
         for Hz in list(SoundIntensity)}

DataInfo = dict((Name, eval(Name)) for Name in ['AnimalName', 'Rate', 
                                       'BaudRate', 'SoundPrePauseDur', 
                                       'SoundPulseDur', 'SoundPostPauseDur', 
                                       'SoundPulseNo', 'SoundStimBlockNo', 
                                       'SoundPauseBetweenStimBlocksDur', 
                                       'Intensities', 'NoiseFrequency', 
                                       'LaserPrePauseDur', 'LaserPulseDur', 
                                       'LaserPostPauseDur', 'LaserPulseNo', 
                                       'LaserStimBlockNo', 
                                       'LaserPauseBetweenStimBlocksDur', 
                                       'CalibrationFile', 'FileName'])

with h5py.File(FileName) as h5:
    h5.create_group('DataInfo')
    h5.create_group('ExpInfo')
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


#%% Prepare sound stimulation
Sound, SoundPauseBetweenStimBlocks, _ = \
    ControlSoundBoard.GenSound(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, 
                               NoiseFrequency, TTLAmpF, CalibrationFile, 
                               SoundBoard, SoundPrePauseDur, SoundPostPauseDur, 
                               SoundStimBlockNo, 
                               SoundPauseBetweenStimBlocksDur)


#%% Prepare laser stimulation
Laser, LaserPauseBetweenStimBlocks, _ = \
    ControlSoundBoard.GenLaser(Rate, LaserPulseDur, LaserPulseNo, TTLAmpF, 
                               CalibrationFile, SoundBoard, LaserPrePauseDur, 
                               LaserPostPauseDur, LaserStimBlockNo, 
                               LaserPauseBetweenStimBlocksDur)


#%% Prepare sound and laser simultaneous stimulation
SoundAndLaser, SoundAndLaserPauseBetweenStimBlocks, _ = \
    ControlSoundBoard.GenSoundLaser(Rate, SoundPulseDur, SoundPulseNo, 
                                    SoundAmpF, NoiseFrequency, LaserPulseDur, 
                                    LaserPulseNo, TTLAmpF, CalibrationFile, 
                                    SoundBoard, SoundPrePauseDur, 
                                    SoundPostPauseDur, SoundStimBlockNo, 
                                    SoundPauseBetweenStimBlocksDur, 
                                    LaserPrePauseDur, LaserPostPauseDur, 
                                    LaserStimBlockNo, 
                                    LaserPauseBetweenStimBlocksDur)


#%% Run sound 135-160
#Hz = input('Choose Freq index: ')
#DVCoord = input('Choose DVCoord (in µm): '); 
DVCoord = 'Out'
Hz = 3
Hz = int(Hz)

print('Running...')
Key = str(NoiseFrequency[Hz][0]) + '-' + str(NoiseFrequency[Hz][1])
for AmpF in range(len(SoundAmpF[Key])):
    Arduino.write(b'P')
    for OnePulse in range(SoundPulseNo):
#        print(str(SoundPulseNo), end='')
#        Arduino.write(b'a')
        Stimulation.write(Sound[Hz][AmpF])
#        Arduino.write(b'z')
#        print('Finished AmpF', str(AmpF))
    
    Arduino.write(b'P')
    Stimulation.write(SoundPauseBetweenStimBlocks)


print('Done. Saving info...')
#ExpFileName = AnimalName + '-SoundStim.hdf5'

with h5py.File(FileName) as h5:
    h5['ExpInfo'].create_group(str(len(list(h5['ExpInfo']))))
    Key = list(h5['ExpInfo'].keys())[-1]
    
    h5['ExpInfo'][Key].attrs['StimType'] = [np.string_('Sound')]
    h5['ExpInfo'][Key].attrs['DVCoord'] = DVCoord
    h5['ExpInfo'][Key].attrs['Hz'] = Hz

print('Saved.')
print('Played Freq ' + str(Hz) + ' at ' + DVCoord + 'µm DV')

#%% Run laser 45 85
#DVCoord = input('Choose DVCoord (in µm): '); 
DVCoord = 'Out'
#
print('Running...')
Arduino.write(b'P')
for OneBlock in range(LaserStimBlockNo):
    for OnePulse in range(LaserPulseNo):
#        Arduino.write(b'b')
        Stimulation.write(Laser)
#        Arduino.write(b'y')
    
    Stimulation.write(LaserPauseBetweenStimBlocks)
Arduino.write(b'P')

#print('Done. Saving info...')
#lHz = 1000/round((LaserPulseDur+LaserPostPauseDur)*1000)
#with h5py.File(FileName) as h5:
#    h5.create_group(str(len(list(h5)) - 1))
#    h5[list(h5.keys())[-2]].attrs['StimType'] = [np.string_('Laser')]
#    h5[list(h5.keys())[-2]].attrs['DVCoord'] = DVCoord
#    h5[list(h5.keys())[-2]].attrs['lHz'] = lHz
#
#print('Saved.')
#print('Ran laser pulses at ' + str(lHz) + ' at ' + DVCoord + 'µm DV')


#%% Run sound and laser >200
#Hz = input('Choose Freq index: ')
#DVCoord = input('Choose DVCoord (in µm): '); 
DVCoord = 'Out'
Hz = 0
Hz = int(Hz)

print('Running...')
Key = str(NoiseFrequency[Hz][0]) + '-' + str(NoiseFrequency[Hz][1])
for AmpF in range(len(SoundAmpF[Key])):
    Arduino.write(b'P')
    for OnePulse in range(SoundPulseNo):
        Stimulation.write(SoundAndLaser[Hz][AmpF])

    Arduino.write(b'P')
    Stimulation.write(SoundAndLaserPauseBetweenStimBlocks)
 
#print('Done. Saving info...')
#lHz = 1000/round((LaserPulseDur+LaserPostPauseDur)*1000)
#with h5py.File(FileName) as h5:
#    h5.create_group(str(len(list(h5)) - 1))
#    h5[list(h5.keys())[-2]].attrs['StimType'] = [np.string_('Sound'), 
#                                                 np.string_('Laser')]
#    h5[list(h5.keys())[-2]].attrs['DVCoord'] = DVCoord
#    h5[list(h5.keys())[-2]].attrs['Hz'] = Hz
#    h5[list(h5.keys())[-2]].attrs['lHz'] = lHz
#
#print('Saved.')
#print('Played Freq ' + str(Hz) + ' at ' + DVCoord + 'µm DV')
#print('Ran laser pulses at ' + str(lHz) + ' at ' + DVCoord + 'µm DV')