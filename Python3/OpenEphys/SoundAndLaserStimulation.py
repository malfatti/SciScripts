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
AnimalName = 'LongEvansTest01'
Rate = 192000
BaudRate = 38400

#CalibrationFile = '/home/cerebro/Malfatti/Test/' + \
#                  '20160419093139-SoundMeasurement/' + \
#                  '20160419093139-SoundMeasurement.hdf5'
CalibrationFile = '/home/malfatti/Documents/PhD/Tests/' + \
                  '20160419093139-SoundMeasurement/' + \
                  '20160419093139-SoundMeasurement.hdf5'

# Sound board used
#SoundBoard = 'USBPre2_oAux-iAux'
SoundBoard = 'Intel_oAnalog-iAnalog'

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 1

## Fill all durations in SECONDS!

## Sound
# Silence before pulse
SoundPrePauseDur = 0.004
# Pulse duration
SoundPulseDur = 0.005
# Silence after pulse
SoundPostPauseDur = 0.091
# Amount of pulses per block
SoundPulseNo = 529
# Number of blocks
SoundStimBlockNo = 1
# Duration of pause between blocks
SoundPauseBetweenStimBlocksDur = 10
# Intensities tested, in order, in dB. Supports floats :)
#Intensities = [80, 75, 70, 65, 60, 55, 50, 45, 40, 35]
Intensities = [80, 70, 60, 50, 40]
#Intensities = [80]
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
# USE ONLY FREQUENCY BANDS THAT WERE CALIBRATED. To check the calibrated freqs, 
# just run the cell once and then list(SoundIntensity).
NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], 
                  [12000, 14000], [14000, 16000]]
#NoiseFrequency = [[20000, 30000], [40000, 70000]]

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
import Hdf5F
import numpy as np

Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-SoundStim.hdf5'])

SoundAmpF = ControlSoundBoard.dBToAmpF(Intensities, CalibrationFile)

DataInfo = dict((Name, eval(Name)) 
                for Name in ['AnimalName', 'Rate', 'BaudRate', 
                             'SoundPrePauseDur', 'SoundPulseDur', 
                             'SoundPostPauseDur', 'SoundPulseNo', 
                             'SoundStimBlockNo', 
                             'SoundPauseBetweenStimBlocksDur', 'Intensities', 
                             'NoiseFrequency', 
#                             'LaserPrePauseDur', 'LaserPulseDur', 
#                             'LaserPostPauseDur', 'LaserPulseNo', 
#                             'LaserStimBlockNo', 
#                             'LaserPauseBetweenStimBlocksDur', 
                             'CalibrationFile', 'FileName'])

Hdf5F.WriteDict(DataInfo, '/DataInfo', FileName)
Hdf5F.WriteDict(SoundAmpF, '/DataInfo/SoundAmpF', FileName)

#Arduino = ControlArduino.CreateObj(BaudRate)


#%% Prepare sound stimulation
Sound, PlaySound = ControlSoundBoard.SoundStim(Rate, SoundPulseDur, 
                                               SoundPulseNo, SoundAmpF, 
                                               NoiseFrequency, TTLAmpF, 
                                               SoundBoard, 'AllPulses', 
                                               SoundPrePauseDur, 
                                               SoundPostPauseDur, 
                                               SoundStimBlockNo, 
                                               SoundPauseBetweenStimBlocksDur)

Stimulation = ControlSoundBoard.GenAudioObj(Rate, 'out')
SoundPauseBetweenStimBlocks = ControlSoundBoard.GenPause(
                                        Rate, SoundPauseBetweenStimBlocksDur
                                                              )

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


#%% Run sound
DVCoord = '3000'
#Freq = 4
#Freq = int(Freq)

FKeys = list(Sound.keys()); ToPrepend = []
for FF in ['8000-10000', '9000-11000']:
    if FF in FKeys:
        del(FKeys[FKeys.index(FF)]); ToPrepend.append(FF)

ToPrepend.sort(); FKeys = ToPrepend + FKeys

while True:
    print('Remember to change folder name in OE!')
    print('Choose frequency:')
    for Ind, K in enumerate(FKeys):
        print(str(Ind), '=' , K)
    
    print(str(len(FKeys)), '=', 'Cancel')
    FKey = input(': ')
    
    if FKey == str(len(FKeys)): break
    
    try:
        FKey = FKeys[int(FKey)]
    except IndexError:
        print('=== Wrong Freq index. Stopping... ===')
        print('')
        break
    
    AKeys = list(Sound[FKey].keys()); AKeys = sorted(AKeys, reverse=True)
    for AmpF, AKey in enumerate(AKeys):
        print('Playing', FKey, 'at', str(Intensities[AmpF]), 'dB')
        
#        Arduino.write(b'P')
        SD
        PlaySound(Freq, AmpF)
        Stimulation.write(SoundPauseBetweenStimBlocks)
#        Arduino.write(b'P')
    
    Hdf5F.WriteExpInfo('Sound', DVCoord, Freq, FileName)
    print('Played Freq', str(Freq), 'at', DVCoord, 'µm DV')

#%% Run laser
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


#%% Run sound and laser
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