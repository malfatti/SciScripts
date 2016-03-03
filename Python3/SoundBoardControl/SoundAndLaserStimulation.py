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
AnimalName = 'TestSetup01'
Rate = 128000

CalibrationFile = '/home/malfatti/NotSynced/SoftwareTest/' + \
                  'SoundMeasurements/20160125114052-SoundMeasurement/' + \
                  'SoundIntensity'
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
SoundPulseNo = 200
# Number of blocks
SoundStimBlockNo = 1
# Duration of pause between blocks
SoundPauseBetweenStimBlocksDur = 5
# Amplification factor (consider using values <= 1). If using one value, keep it in a list.
SoundAmpF = [0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
NoiseFrequency = [[8000, 10000], [10000, 12000], [12000, 14000], [14000, 16000]]

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

import ControlSoundBoard
import datetime
import pyaudio
import shelve

with shelve.open(CalibrationFile) as Shelve:
    SoundIntensity = Shelve['SoundIntensity']

Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-SoundStim'])

DataInfo = dict((Name, eval(Name)) for Name in ['AnimalName', 
                                       'SoundPrePauseDur', 'SoundPulseDur', 
                                       'SoundPostPauseDur', 'SoundPulseNo', 
                                       'SoundStimBlockNo', 
                                       'SoundPauseBetweenStimBlocksDur', 
                                       'SoundAmpF', 'NoiseFrequency', 
                                       'LaserPrePauseDur', 'LaserPulseDur', 
                                       'LaserPostPauseDur', 'LaserPulseNo',
                                       'LaserStimBlockNo', 
                                       'LaserPauseBetweenStimBlocksDur', 
                                       'FileName'])

with shelve.open(FileName) as Shelve:
    Shelve['DataInfo'] = DataInfo
del(Date, FileName, DataInfo)

p = pyaudio.PyAudio()
Stimulation = p.open(format=pyaudio.paFloat32,
                channels=2,
                rate=Rate,
                output=True)


#%% Prepare sound stimulation
Sound, SoundPauseBetweenStimBlocks, _ = \
    ControlSoundBoard.GenSound(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, 
                               NoiseFrequency, TTLAmpF, SoundPrePauseDur, 
                               SoundPostPauseDur, SoundStimBlockNo, 
                               SoundPauseBetweenStimBlocksDur)


#%% Prepare laser stimulation
Laser, LaserPauseBetweenStimBlocks, _ = \
    ControlSoundBoard.GenLaser(Rate, LaserPrePauseDur, LaserPulseDur, 
                               LaserPostPauseDur, LaserPulseNo, 
                               LaserStimBlockNo, 
                               LaserPauseBetweenStimBlocksDur, TTLAmpF)


#%% Prepare sound and laser simultaneous stimulation
SoundAndLaser, SoundAndLaserPauseBetweenStimBlocks, _ = \
    ControlSoundBoard.GenSoundLaser(Rate, SoundPrePauseDur, SoundPulseDur, 
                                    SoundPostPauseDur, SoundPulseNo, 
                                    SoundStimBlockNo, 
                                    SoundPauseBetweenStimBlocksDur, SoundAmpF, 
                                    NoiseFrequency, LaserPrePauseDur, 
                                    LaserPulseDur, LaserPostPauseDur, 
                                    LaserPauseBetweenStimBlocksDur, TTLAmpF)


#%% Run sound
Date = datetime.datetime.now()

Freq = 0; AmpF = 0
NoiseFreq = NoiseFrequency[Freq]; SoundAmp = SoundAmpF[AmpF]
for Block in range(SoundStimBlockNo):
    for Pulse in range(SoundPulseNo):
        Stimulation.write(Sound[Freq][AmpF])
    
    Stimulation.write(SoundPauseBetweenStimBlocks)

print('Done. Saving info...')
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-SoundStim-', str(NoiseFreq[0]), '_', str(NoiseFreq[1]), 
                    '-', str(SoundIntensity[Freq][str(SoundAmpF[AmpF])])[:5], 
                    'dB'])

ExpInfo = dict((Name, eval(Name)) for Name in ['Freq', 'AmpF', 'SoundAmp', 
                                               'NoiseFreq', 'FileName'])

with shelve.open(FileName) as Shelve:
    Shelve['ExpInfo'] = ExpInfo
del(Date, FileName, ExpInfo)
print('Saved.')

#%% Run laser
for Block in range(LaserStimBlockNo):
    for Pulse in range(LaserPulseNo):
        Stimulation.write(Laser)
    
    Stimulation.write(LaserPauseBetweenStimBlocks)

print('Done. No need to save experiment info.')


#%% Run sound and laser
Date = datetime.datetime.now()

Freq = 0; AmpF = 0
NoiseFreq = NoiseFrequency[Freq]; SoundAmp = SoundAmpF[AmpF]
for Block in range(SoundStimBlockNo):
    for Pulse in range(SoundPulseNo):
        Stimulation.write(SoundAndLaser[Freq][AmpF])

    Stimulation.write(SoundAndLaserPauseBetweenStimBlocks)

ExpInfo = dict((Name, eval(Name)) for Name in ['Freq', 'AmpF', 'SoundAmp', 
                                               'NoiseFreq', 'FileName'])

print('Done. Saving info...')
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-SoundAndLaserStim-', str(NoiseFreq[0]), '_', 
                    str(NoiseFreq[1]), '-', 
                    str(SoundIntensity[Freq][str(SoundAmpF[AmpF])])[:5], 'dB'])

with shelve.open(FileName) as Shelve:
    Shelve['ExpInfo'] = ExpInfo
del(Date, FileName, ExpInfo)
